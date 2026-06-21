"""
make_liver_cohort_figure.py
===========================
Figure for the in-vivo COHORT validation (liver, N patients).

Panel A: per-patient held-out-b coverage at nominal 0.95, NPE vs NLLS (paired).
Panel B: cohort-mean reliability curve (empirical vs nominal) for NPE and NLLS.

Reads npe/liver_cohort_coverage.csv + npe/liver_cohort_summary.json.
Writes figures/manuscript/figS6_liver_cohort.{png,pdf} and a small CSV of the
plotted aggregates.
"""
from __future__ import annotations

import os
import csv
import json
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LEVELS = [0.50, 0.68, 0.80, 0.90, 0.95]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="npe/liver_cohort_coverage.csv")
    ap.add_argument("--json", default="npe/liver_cohort_summary.json")
    ap.add_argument("--out-dir", default="figures/manuscript")
    args = ap.parse_args()

    rows = []
    with open(args.csv) as f:
        for r in csv.DictReader(f):
            rows.append(r)
    summary = json.load(open(args.json))
    n = len(rows)

    npe95 = np.array([float(r["npe_cov_0.95"]) for r in rows])
    nlls95 = np.array([float(r["nlls_cov_0.95"]) for r in rows])

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.6))

    # Panel A: paired per-patient coverage at 0.95
    rng = np.random.default_rng(0)
    xj = lambda x0: x0 + rng.uniform(-0.06, 0.06, size=n)
    for i in range(n):
        axA.plot([1, 2], [npe95[i], nlls95[i]], color="0.8", lw=0.5, zorder=1)
    axA.scatter(xj(1), npe95, s=18, color="#c0392b", label="NPE", zorder=3)
    axA.scatter(xj(2), nlls95, s=18, color="#27ae60", label="NLLS", zorder=3)
    axA.axhline(0.95, ls="--", color="k", lw=1, label="nominal 0.95")
    axA.set_xticks([1, 2]); axA.set_xticklabels(["NPE\n(amortized)", "NLLS\n(per-voxel)"])
    axA.set_ylabel("Held-out-b coverage @ nominal 0.95")
    axA.set_ylim(0, 1.0)
    axA.set_title(f"A  Per-patient coverage (N={n} patients)\nNPE < NLLS in "
                  f"{summary['frac_patients_npe_undercovers_nlls']*100:.0f}% of patients")
    axA.legend(loc="lower left", fontsize=8, framealpha=0.9)

    # Panel B: cohort-mean reliability curve
    npe_mean = [summary["per_level_mean"][str(l)]["npe"] for l in LEVELS]
    nlls_mean = [summary["per_level_mean"][str(l)]["nlls"] for l in LEVELS]
    npe_sd = [np.std([float(r[f"npe_cov_{l}"]) for r in rows]) for l in LEVELS]
    nlls_sd = [np.std([float(r[f"nlls_cov_{l}"]) for r in rows]) for l in LEVELS]
    axB.plot([0.4, 1.0], [0.4, 1.0], ls="--", color="k", lw=1, label="perfect calibration")
    axB.errorbar(LEVELS, npe_mean, yerr=npe_sd, marker="o", color="#c0392b",
                 capsize=3, label="NPE (amortized)")
    axB.errorbar(LEVELS, nlls_mean, yerr=nlls_sd, marker="s", color="#27ae60",
                 capsize=3, label="NLLS (per-voxel)")
    axB.set_xlabel("Nominal credibility level")
    axB.set_ylabel("Cohort-mean empirical coverage")
    axB.set_xlim(0.45, 1.0); axB.set_ylim(0.0, 1.0)
    axB.set_title("B  Cohort-mean reliability (liver, GE 1.5T)")
    axB.legend(loc="upper left", fontsize=8, framealpha=0.9)

    fig.suptitle("In-vivo cohort validation: off-scheme NPE is systematically more "
                 "overconfident than per-voxel NLLS", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    os.makedirs(args.out_dir, exist_ok=True)
    png = os.path.join(args.out_dir, "figS6_liver_cohort.png")
    pdf = os.path.join(args.out_dir, "figS6_liver_cohort.pdf")
    fig.savefig(png, dpi=300); fig.savefig(pdf)
    print(f"Wrote {png} and {pdf}")

    # plotted aggregates CSV
    agg = os.path.join(args.out_dir, "figS6_liver_cohort.csv")
    with open(agg, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["nominal", "npe_mean", "npe_sd", "nlls_mean", "nlls_sd"])
        for i, l in enumerate(LEVELS):
            w.writerow([l, npe_mean[i], npe_sd[i], nlls_mean[i], nlls_sd[i]])
    print(f"Wrote {agg}")


if __name__ == "__main__":
    main()
