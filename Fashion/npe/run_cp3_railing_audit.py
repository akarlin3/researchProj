"""
run_cp3_railing_audit.py  (Friction remediation — Checkpoint 3, HC5)
====================================================================
NLLS-baseline survivorship audit. The reviewer concern (HC5): the in-vivo NLLS
spread comparison silently *excludes* the boundary-railed NLLS voxels, so the
reported D* spread ratio is a survivorship-filtered number.

This script recomputes the Supplementary-Figure-S4 abdominal spread comparison
BOTH ways -- railed-included (all voxels) and railed-excluded (non-railed only)
-- directly from the committed per-voxel CSV, and prints them side by side so
the inclusion policy is explicit and the OLD (railed-excluded) numbers stay
labelled.

It also documents the brain (N=500) held-out-b coverage inclusion policy, which
is verified from run_f_realdata.py (not recomputed here -- raw in-vivo signals
are not committed): all 500 sampled GM voxels are kept, railed NLLS fits are
included, and failed (NaN) fits are counted as NOT covered (they lower, not
inflate, the NLLS coverage). The brain "0.90" is therefore already
railed-included; there is no survivorship filter to undo on that number.

Outputs: cp3_railing_audit.csv + printed summary.
"""
from __future__ import annotations

import argparse
import csv
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CSV = os.path.join(HERE, "..", "figures", "manuscript", "figS4_invivo_illustration.csv")


def iqr(x):
    return float(np.subtract(*np.percentile(x, [75, 25])))


def main():
    ap = argparse.ArgumentParser(description="CP3 NLLS railing/survivorship audit (abdominal S4).")
    ap.add_argument("--csv", default=DEFAULT_CSV)
    ap.add_argument("--out", default=os.path.join(HERE, "cp3_railing_audit.csv"))
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.csv)))
    nlls = np.array([float(r["nlls_dstar"]) for r in rows])
    npe = np.array([float(r["npe_dstar"]) for r in rows])
    railed = np.array([int(r["nlls_railed"]) for r in rows]).astype(bool)
    n = len(rows)

    def block(mask, label):
        nl, np_ = nlls[mask], npe[mask]
        return {
            "subset": label, "n": int(mask.sum()),
            "nlls_sd": float(np.std(nl)), "npe_sd": float(np.std(np_)),
            "sd_ratio_npe_over_nlls": float(np.std(np_) / np.std(nl)),
            "nlls_iqr": iqr(nl), "npe_iqr": iqr(np_),
            "iqr_ratio_npe_over_nlls": float(iqr(np_) / iqr(nl)),
            "nlls_median": float(np.median(nl)), "npe_median": float(np.median(np_)),
        }

    recs = [block(np.ones(n, bool), "all_voxels_railed_INCLUDED"),
            block(~railed, "non_railed_railed_EXCLUDED_OLD")]

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(recs[0].keys()))
        w.writeheader()
        for r in recs:
            w.writerow(r)

    print(f"Abdominal in-vivo case (Supp Fig S4): n={n}, railed={railed.sum()} "
          f"({railed.mean()*100:.1f}%)\n")
    print(f"{'subset':<32} {'n':>5} {'SD ratio':>9} {'IQR ratio':>10} "
          f"{'NLLS SD':>9} {'NPE SD':>8}")
    for r in recs:
        print(f"{r['subset']:<32} {r['n']:>5} {r['sd_ratio_npe_over_nlls']:>9.3f} "
              f"{r['iqr_ratio_npe_over_nlls']:>10.3f} {r['nlls_sd']:>9.4f} {r['npe_sd']:>8.4f}")
    print(f"\n-> per-voxel audit written to {args.out}")
    print("\nInclusion-policy note (brain N=500 held-out-b coverage, from run_f_realdata.py):")
    print("  all 500 GM voxels kept; railed NLLS fits included; failed (NaN) fits count as")
    print("  NOT covered. The NLLS held-out-b coverage (0.904 @0.95) is already railed-included.")


if __name__ == "__main__":
    main()
