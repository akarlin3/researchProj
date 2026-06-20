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

from sextant.bootstrap import DEFAULT_B, bootstrap_fraction      # noqa: E402
from sextant.cohorts import load_nifti_cohort                    # noqa: E402
from sextant.railing import analyze_cohort                       # noqa: E402
from sextant.seeding import GLOBAL_SEED                          # noqa: E402

_EXTRACT = os.path.join(_ROOT, "data", "osipi", "extracted")
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


def run_cohort(name, img, mask, bval, n_boot):
    coh = load_nifti_cohort(name, img, mask, bval)
    tight = analyze_cohort(coh, bounds="tight")
    wide = analyze_cohort(coh, bounds="wide")
    ci = bootstrap_fraction(tight.railed_hi, n_boot=n_boot)
    return {
        "name": name,
        "n_voxels": coh.n_voxels,
        "n_high_snr": coh.n_high_snr,
        "tight": tight.to_dict(),
        "wide": wide.to_dict(),
        "bootstrap_ci": ci.to_dict(),
        "verdict": verdict(ci.point, ci.lo),
        "mask": os.path.basename(mask),
    }


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

    results = []
    for name, img, mask, bval in cohorts:
        print(f"[railing] {name} ...", flush=True)
        r = run_cohort(name, img, mask, bval, args.n_boot)
        t = r["tight"]; ci = r["bootstrap_ci"]
        print(f"  railed(tight)={t['frac_railed']:.4f}  "
              f"CI[{ci['lo']:.4f},{ci['hi']:.4f}]  "
              f"lower={t['frac_lower']:.3f}/upper={t['frac_upper']:.3f}  "
              f"wide={r['wide']['frac_railed']:.4f}  -> {r['verdict']['label']}")
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
