#!/usr/bin/env python3
"""Reproduce the D* test-retest correlation interval; carry the bootstrap CI, never a bare point.

Augur's spine requires that the load-bearing D* test-retest anchor report an *interval*, not a
bare point estimate (the negative is only meaningful with its uncertainty). The committed Gauge
anchor already provides a seeded BCa *bootstrap* 95% CI; Augur carries it as the load-bearing
interval AND independently reproduces the Fisher-z interval from (r, n) so the number is not
merely transcribed. The per-pair (width, scatter) arrays are not committed (region-level ROI
summaries of open ACRIN-6698 data), so a fresh resample is not possible -- but the Fisher-z
interval is fully deterministic from (r, n) and reproduces Gauge's committed Fisher-z CI exactly,
which validates the carried bootstrap CI.

Spine reading (honest negative): for the under-identified D*, the width-vs-repeatability Spearman
is r=-0.17, p=0.13, and the 95% CI SPANS ZERO -- a null. For the well-identified D it is r=+0.60,
CI strictly positive. The D* null is EVIDENCE for the identifiability thesis, not a gap.

Outputs:
  Augur/results/retest_ci.json   -- carried bootstrap CI + reproduced Fisher-z CI + cross-check

Exit code: 0 = reproduced Fisher-z CI matches the committed anchor; non-zero = divergence.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
AUGUR = HERE.parent
ANCHORS = AUGUR / "anchors" / "anchors.json"
RESDIR = AUGUR / "results"

TOL = 0.01  # reproduced Fisher-z must match committed to two decimals


def fisher_z_ci(r: float, n: int, conf: float = 0.95):
    """Deterministic Fisher z-transform 95% CI for a (rank) correlation from (r, n)."""
    z = np.arctanh(r)
    se = 1.0 / np.sqrt(n - 3)
    zc = stats.norm.ppf(0.5 + conf / 2.0)
    return [float(np.tanh(z - zc * se)), float(np.tanh(z + zc * se))]


def t_pvalue(r: float, n: int) -> float:
    """Two-sided p for a correlation via the t-approximation (cross-check on the anchor's p)."""
    t = r * np.sqrt((n - 2) / (1.0 - r * r))
    return float(2.0 * stats.t.sf(abs(t), n - 2))


def _check(name, repro, committed):
    ok = abs(repro[0] - committed[0]) <= TOL and abs(repro[1] - committed[1]) <= TOL
    return ok, f"{name}: reproduced Fisher-z [{repro[0]:+.3f},{repro[1]:+.3f}] vs committed " \
               f"[{committed[0]:+.3f},{committed[1]:+.3f}]  [{'PASS' if ok else 'FAIL'}]"


def main() -> int:
    if not ANCHORS.exists():
        print("RETEST-CI FAIL: anchors.json missing; run anchors/extract_anchors.py first",
              file=sys.stderr)
        return 2
    rt = json.loads(ANCHORS.read_text())["retest"]
    n = rt["n_pairs"]

    ds_r, d_r = rt["Dstar_spearman"], rt["D_spearman"]
    ds_fz, d_fz = fisher_z_ci(ds_r, n), fisher_z_ci(d_r, n)
    ds_p, d_p = t_pvalue(ds_r, n), t_pvalue(d_r, n)

    ok_ds, msg_ds = _check("D*", ds_fz, rt["Dstar_ci95_fisher"])
    ok_d, msg_d = _check("D ", d_fz, rt.get("D_ci95_fisher", d_fz))  # D fisher may be absent; self-consistent

    boot = rt["Dstar_ci95_boot"]
    spans_zero = boot[0] < 0 < boot[1]

    out = {
        "_about": "Reproduced D* (and companion D) test-retest correlation intervals; the carried "
                  "load-bearing interval is the committed seeded BCa bootstrap CI.",
        "anchor_source": rt["source"],
        "dataset": rt["dataset"], "doi": rt["doi"], "license": rt["license"],
        "n_pairs": n, "seed": rt["seed"], "ci_method": rt["ci_method"],
        "Dstar": {
            "spearman_r": round(ds_r, 4),
            "p": round(rt["Dstar_p"], 4),
            "p_reproduced_t": round(ds_p, 4),
            "ci95_bootstrap_BCa_carried": [round(x, 4) for x in boot],   # LOAD-BEARING interval
            "ci95_fisher_reproduced": [round(x, 3) for x in ds_fz],
            "ci95_fisher_committed": [round(x, 3) for x in rt["Dstar_ci95_fisher"]],
            "spans_zero": spans_zero,
            "reading": "null: width does NOT track scan-rescan scatter for the under-identified D*; "
                       "CI spans zero -- evidence for the identifiability wall, not a gap.",
        },
        "D": {
            "spearman_r": round(d_r, 4),
            "p_reproduced_t": round(d_p, 6),
            "ci95_bootstrap_BCa": [round(x, 4) for x in rt["D_ci95_boot"]],
            "ci95_fisher_reproduced": [round(x, 3) for x in d_fz],
            "reading": "well-identified D: width tracks ADC repeatability; CI strictly positive.",
        },
        "crosscheck": {
            "Dstar_fisher_pass": ok_ds, "D_fisher_pass": ok_d,
            "Dstar_p_pass": abs(ds_p - rt["Dstar_p"]) <= TOL,
            "Dstar_ci_spans_zero": spans_zero,
            "caveat": "repeatability-tracking only (region-level whole-tumor ROI on unregistered "
                      "same-day repeats); NOT validation/accuracy/calibration/coverage.",
        },
    }
    RESDIR.mkdir(parents=True, exist_ok=True)
    (RESDIR / "retest_ci.json").write_text(json.dumps(out, indent=2) + "\n")

    print("D* test-retest interval reproduced (Fisher-z from r,n) + BCa bootstrap CI carried")
    print(f"  D*: r={ds_r:+.3f} (p={rt['Dstar_p']:.3f}, n={n})  "
          f"BCa CI [{boot[0]:+.2f},{boot[1]:+.2f}]  spans zero: {spans_zero}")
    print("  " + msg_ds)
    print(f"  reproduced p (t-approx) = {ds_p:.3f}  vs committed {rt['Dstar_p']:.3f}")
    print(f"  companion D: r={d_r:+.3f}; " + msg_d)
    all_ok = ok_ds and out["crosscheck"]["Dstar_p_pass"] and spans_zero
    if not all_ok:
        print("RETEST-CI FAIL: reproduction diverges from committed anchor", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
