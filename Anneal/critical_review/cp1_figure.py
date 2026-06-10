"""CP1 supplementary figure + table: within-stratum Weibull shape k(bin) per N,
for the primary (|Delta phi_0|) and decisive (reduced predicted lifetime) binnings.
All numbers come from cp1_stratified.json (produced by cp1_stratified.py)."""
from __future__ import annotations

import csv
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "cp1_stratified")
blob = json.load(open(os.path.join(OUT, "cp1_stratified.json")))

PANELS = [("absdphi0", r"$|\Delta\phi_0|$  (collective phase gap)"),
          ("t_capture", r"reduced-model predicted lifetime  $t_{\rm capture}$")]

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
Ns = sorted((int(n) for n in blob["pooled"]))
cmap = plt.cm.viridis(np.linspace(0, 0.92, len(Ns)))

for ax, (var, xlabel) in zip(axes, PANELS):
    strat = blob["stratified"][var]["by_N"]
    for c, N in zip(cmap, Ns):
        cells = strat[str(N)]
        xs = [cc["x_mean"] for cc in cells]
        ks = [cc["k"] for cc in cells]
        lo = [cc["k"] - cc["k_lo"] for cc in cells]
        hi = [cc["k_hi"] - cc["k"] for cc in cells]
        ax.errorbar(xs, ks, yerr=[lo, hi], marker="o", ms=4, lw=1.3, capsize=2,
                    color=c, label=f"N={N}")
    ax.axhline(1.0, color="crimson", ls="--", lw=1, zorder=0)
    ax.set_xlabel(xlabel)
    ax.set_title(f"within-stratum $k$  vs  {var}", fontsize=10)
    if var == "t_capture":
        ax.set_xscale("log")
axes[0].set_ylabel(r"Weibull shape $\hat k$ (profile-likelihood 95% CI)")
axes[0].legend(fontsize=7, ncol=2, frameon=False)
fig.suptitle("CP1 (C1): aging shape $k>1$ survives conditioning on the collective initial "
             "condition at $A=0.5$\n(every stratum: $\\hat k>1$, 95% CI excludes 1; "
             "red dashes = memoryless $k=1$)", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(OUT, "cp1_k_vs_bin.png"), dpi=150)
fig.savefig(os.path.join(OUT, "cp1_k_vs_bin.pdf"))
print("[fig] cp1_k_vs_bin.png/.pdf")

# flat CSV table of every stratified cell
with open(os.path.join(OUT, "cp1_strata_table.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["binning_var", "N", "bin", "x_mean", "edge_lo", "edge_hi", "n", "n_events",
                "censor_frac", "k", "k_lo", "k_hi", "ci_excludes_1", "cv_tau", "median_tau",
                "lrt_p"])
    for var in ("absdphi0", "t_capture", "Rincoh0"):
        for N in sorted(blob["stratified"][var]["by_N"], key=int):
            for cc in blob["stratified"][var]["by_N"][N]:
                w.writerow([var, N, cc["bin"], f"{cc['x_mean']:.4f}", f"{cc['edge_lo']:.4f}",
                            f"{cc['edge_hi']:.4f}", cc["n"], cc["n_events"],
                            f"{cc['censor_frac']:.3f}", f"{cc['k']:.3f}", f"{cc['k_lo']:.3f}",
                            f"{cc['k_hi']:.3f}", int(cc["ci_excludes_1"]),
                            f"{cc['cv_tau_events']:.3f}", f"{cc['median_tau']:.2f}",
                            f"{cc['lrt_p']:.2e}"])
print("[table] cp1_strata_table.csv")

# headline summary line
for var in ("absdphi0", "t_capture", "Rincoh0"):
    v = blob["verdicts"][var]
    print(f"  {var:>10}: k in [{v['k_min']:.2f},{v['k_max']:.2f}], "
          f"CI excl 1 in {v['n_ci_excl_1']}/{v['n_cells']}")
