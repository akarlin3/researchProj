#!/usr/bin/env python
"""Boundary-railing driver (CP2 primary analysis + CP3 replication).

Runs the read-only-reused railing diagnostic on each configured human-abdominal
cohort, with bootstrap CIs and a tight/wide-bounds sensitivity check, and applies
the pre-registered replication thresholds. Writes a single seeded results JSON
that the manuscript's consistency gate reads.

Cohorts (OSIPI human-abdominal, CC-BY-4.0; download with scripts/fetch_osipi.py):
  * abdomen_homogeneous -- the original curated ROI (reproduces Fashion's 54.7%).
  * abdomen_full        -- the full-abdomen ROI (internal generalisation test).
Additional independent cohorts (e.g. TCGA-LIHC) can be appended via --extra.

Run:  python scripts/run_railing.py
"""
from __future__ import annotations

import argparse
import datetime
import json
import os
import sys

# Make the sextant package importable without installation.
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "sextant-core"))

import glob                                                     # noqa: E402

import numpy as np                                              # noqa: E402

from sextant.bootstrap import DEFAULT_B, bootstrap_fraction      # noqa: E402
from sextant.cohorts import Cohort, load_array_cohort, load_nifti_cohort  # noqa: E402
from sextant.railing import analyze_cohort                       # noqa: E402
from sextant.seeding import GLOBAL_SEED                          # noqa: E402

_EXTRACT = os.path.join(_ROOT, "data", "osipi", "extracted")
_LIHC = os.path.join(_ROOT, "data", "tcga_lihc")
_RESULTS = os.path.join(_ROOT, "results", "railing_results.json")

# --------------------------------------------------------------------------- #
# Pre-registered replication thresholds (fixed BEFORE running; see VERIFICATION).
# --------------------------------------------------------------------------- #
REPL_POINT_MIN = 0.30      # railing must be a substantial minority of fits
REPL_CI_LO_MIN = 0.20      # robustly so (95% bootstrap CI lower bound)
STRONG_LO, STRONG_HI = 0.40, 0.70   # "strong" band around the original 54.7%

OSIPI_COHORTS = [
    ("abdomen_homogeneous", "abdomen.nii.gz", "mask_abdomen_homogeneous.nii.gz", "abdomen.bval"),
    ("abdomen_full",        "abdomen.nii.gz", "mask_full_temp.nii.gz",          "abdomen.bval"),
]


def verdict(point: float, ci_lo: float) -> dict:
    replicates = (point >= REPL_POINT_MIN) and (ci_lo >= REPL_CI_LO_MIN)
    strong = STRONG_LO <= point <= STRONG_HI
    if replicates:
        label = "REPLICATES-STRONG" if strong else "REPLICATES"
    else:
        label = "SCOPE-DOWN"
    return {
        "label": label, "replicates": bool(replicates), "strong": bool(strong),
        "thresholds": {"point_min": REPL_POINT_MIN, "ci_lo_min": REPL_CI_LO_MIN,
                       "strong_band": [STRONG_LO, STRONG_HI]},
    }


def _summary(name, coh, n_boot, family, extra=None, with_wide=True):
    tight = analyze_cohort(coh, bounds="tight")
    wide = analyze_cohort(coh, bounds="wide") if with_wide else None
    ci = bootstrap_fraction(tight.railed_hi, n_boot=n_boot)
    out = {
        "name": name, "family": family,
        "bvals": coh.bvals.tolist(),
        "n_voxels": coh.n_voxels, "n_high_snr": coh.n_high_snr,
        "n_analyzed": tight.n_analyzed, "subsampled": tight.subsampled,
        "tight": tight.to_dict(), "wide": (wide.to_dict() if wide else None),
        "bootstrap_ci": ci.to_dict(), "verdict": verdict(ci.point, ci.lo),
    }
    if extra:
        out.update(extra)
    return out


def run_cohort(name, img, mask, bval, n_boot):
    coh = load_nifti_cohort(name, img, mask, bval)
    return _summary(name, coh, n_boot, family="OSIPI",
                    extra={"mask": os.path.basename(mask)})


def run_lihc(n_boot):
    """Independent replication: pool TCGA-LIHC liver series per acquisition scheme."""
    files = sorted(glob.glob(os.path.join(_LIHC, "*.npz")))
    by_scheme = {}
    for f in files:
        scheme, pid = os.path.basename(f)[:-4].split("__")
        d = np.load(f)
        coh = load_array_cohort(f"{scheme}:{pid}", d["signals4d"], d["bvals"],
                                meta={"patient": pid, "scheme": scheme})
        by_scheme.setdefault(scheme, []).append((pid, coh))

    results = []
    for scheme, lst in sorted(by_scheme.items()):
        per_subject = []
        for pid, coh in lst:
            t = analyze_cohort(coh, bounds="tight", max_high_snr=8000)
            per_subject.append({"patient": pid, "n_high_snr": coh.n_high_snr,
                                "n_analyzed": t.n_analyzed, "frac_railed": t.frac_railed})
        pooled = Cohort(
            name=f"lihc_{scheme}",
            fit_signals=np.vstack([c.fit_signals for _, c in lst]),
            snrs=np.concatenate([c.snrs for _, c in lst]),
            bvals=lst[0][1].bvals, n_voxels=sum(c.n_voxels for _, c in lst),
            meta={"n_subjects": len(lst), "patients": [p for p, _ in lst]})
        results.append(_summary(
            f"lihc_{scheme}", pooled, n_boot, family="TCGA-LIHC",
            with_wide=("4b" in scheme),   # wide sensitivity on the clean 4-b cohort
            extra={"n_subjects": len(lst), "per_subject": per_subject,
                   "note": ("4-b scheme includes b=0 (clean IVIM normalisation)"
                            if "4b" in scheme else
                            "3-b scheme, lowest b=50 (no b=0; normalised by b=50 -- approximate)")}))
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-boot", type=int, default=DEFAULT_B)
    ap.add_argument("--extra", nargs=4, action="append", metavar=("NAME", "IMG", "MASK", "BVAL"),
                    default=[], help="append an independent cohort (e.g. TCGA-LIHC)")
    args = ap.parse_args()

    if not os.path.isdir(_EXTRACT):
        sys.exit("OSIPI data not found; run `python scripts/fetch_osipi.py` first.")

    cohorts = []
    for name, img, mask, bval in OSIPI_COHORTS:
        cohorts.append((name, os.path.join(_EXTRACT, img),
                        os.path.join(_EXTRACT, mask), os.path.join(_EXTRACT, bval)))
    cohorts.extend((n, i, m, b) for n, i, m, b in args.extra)

    def _report(r):
        t = r["tight"]; ci = r["bootstrap_ci"]
        wide = f"{r['wide']['frac_railed']:.4f}" if r.get("wide") else "n/a"
        sub = f" [subsampled {r['n_analyzed']}/{r['n_high_snr']}]" if r.get("subsampled") else ""
        print(f"  [{r['family']}] {r['name']}  b={r['bvals']}  n_hi={r['n_high_snr']}{sub}\n"
              f"    railed(tight)={t['frac_railed']:.4f}  CI[{ci['lo']:.4f},{ci['hi']:.4f}]  "
              f"lower={t['frac_lower']:.3f}/upper={t['frac_upper']:.3f}  "
              f"wide={wide}  -> {r['verdict']['label']}")

    results = []
    for name, img, mask, bval in cohorts:
        print(f"[railing] {name} ...", flush=True)
        r = run_cohort(name, img, mask, bval, args.n_boot)
        _report(r)
        results.append(r)

    if os.path.isdir(_LIHC) and glob.glob(os.path.join(_LIHC, "*.npz")):
        print("[railing] TCGA-LIHC independent replication ...", flush=True)
        for r in run_lihc(args.n_boot):
            _report(r)
            results.append(r)

    out = {
        "meta": {
            "seed": GLOBAL_SEED, "n_boot": args.n_boot,
            "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "fashion_reported_homogeneous": 0.547,
            "fashion_reported_n_high_snr": 1618,
            "registered_thresholds": {
                "point_min": REPL_POINT_MIN, "ci_lo_min": REPL_CI_LO_MIN,
                "strong_band": [STRONG_LO, STRONG_HI]},
        },
        "cohorts": results,
    }
    os.makedirs(os.path.dirname(_RESULTS), exist_ok=True)
    with open(_RESULTS, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"[railing] wrote {_RESULTS}")


if __name__ == "__main__":
    main()
