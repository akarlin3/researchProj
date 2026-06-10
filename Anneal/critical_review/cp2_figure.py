"""CP2 supplementary figure: the three pre-death behaviours in the N=256 cells —
canonical single-arc chimera, multi-headed chimera, and the rare degenerate decay channel.
Space-time rho_k(t,k) + a representative mature snapshot for each."""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cp2_fields import run_field, rho_profile, detect_death_ring, EPS_STD, DT_HOLD, T_MAX

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "cp2_validation")

EXAMPLES = [
    (0.130, 256, 20266309, "canonical single-arc chimera\n(struct. to death, m*=1)"),
    (0.110, 256, 20265110, "multi-headed chimera\n(struct. to death, m*=2)"),
    (0.130, 256, 20266310, "degenerate decay channel (rare, ~2-3%)\n(chimera dissolves early -> limbo)"),
]

fig, axes = plt.subplots(len(EXAMPLES), 2, figsize=(12, 3.0 * len(EXAMPLES)),
                         gridspec_kw={"width_ratios": [3, 1]})
for row, (beta, N, seed, title) in enumerate(EXAMPLES):
    t, rmean, rstd, th, P = run_field(beta, N, seed)
    tau, ev = detect_death_ring(t, rstd, EPS_STD, DT_HOLD, T_MAX)
    di = int(np.searchsorted(t, tau, side="right")) - 1
    field = np.array([rho_profile(th[i], P) for i in range(di + 1)])
    ax = axes[row, 0]
    im = ax.imshow(field.T, aspect="auto", origin="lower", cmap="magma",
                   extent=[t[0], t[di], 0, N], vmin=0, vmax=1)
    ax.axvline(tau, color="cyan", lw=1.2)
    ax.set_ylabel("ring index $k$")
    ax.set_title(f"$\\beta$={beta}, N={N}, seed {seed}, $\\tau$={tau:.0f}:  {title}", fontsize=9)
    if row == len(EXAMPLES) - 1:
        ax.set_xlabel("time $t$")
    plt.colorbar(im, ax=ax, fraction=0.025, label="$\\rho_k$")

    # representative snapshot at mid-life (mature) for the chimeras; for degenerate show the limbo
    i = int(0.55 * di)
    ax = axes[row, 1]
    ax.plot(field[i], color="crimson", lw=1.2)
    ax.set_ylim(0, 1.02); ax.set_xlim(0, N)
    ax.set_title(f"$\\rho_k$ at $t$={t[i]:.0f}\n$\\rho_{{std}}$={rstd[i]:.3f}", fontsize=8)
    if row == len(EXAMPLES) - 1:
        ax.set_xlabel("ring index $k$")

fig.suptitle("CP2 (C2): what the $\\rho_{std}<0.04$ criterion kills at N=256 — "
             "97–98% genuine chimera deaths, 2–3% degenerate channel\n"
             "filtered Weibull $k$ unchanged (0.110: 1.220$\\to$1.241; 0.130: 1.466$\\to$1.481, "
             "CIs exclude 1)", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(os.path.join(OUT, "cp2_representative.png"), dpi=150)
fig.savefig(os.path.join(OUT, "cp2_representative.pdf"))
print("[fig] cp2_representative.png/.pdf")
