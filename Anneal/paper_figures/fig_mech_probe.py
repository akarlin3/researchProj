"""Mechanism-probe figure (Appendix B addendum) — both re-injection candidates.

(a) Multiplicative breath-locked noise (sigma_eff = c_m * (1 - min(rho1,rho2)),
    N-independent amplitude by construction): prolongation factor vs c_m, one
    line per N, with the measured ~3.2x target and the physical amplitude
    (pooled median of the measured fluctuation-to-envelope ratio; shaded band =
    per-N spread of that estimate).
(b) Watanabe-Strogatz constants decomposition: per-seed absorption time of the
    constants-uniformized (Poisson-projected, same per-population Z1) partner
    vs the actual seed, with per-N medians and the reduced-flow medians.

Plot-only: reads noise_results/mech_probe_results.json and
manifold_results/ws_decomposition_N{8,16,32,64}.jsonl (both produced by
committed, seeded scripts). Writes paper_figures/fig_mech_probe.{pdf,png}.
Run: python3 paper_figures/fig_mech_probe.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools/paper-figures"))
import style  # noqa: E402

style.apply_style()
import matplotlib.pyplot as plt  # noqa: E402

P = style.PALETTE
NS = [8, 16, 32, 64]
colors = {8: P["blue"], 16: P["green"], 32: P["orange"], 64: P["vermillion"]}
markers = {8: "o", 16: "s", 32: "^", 64: "D"}

R = json.load(open(os.path.join(ROOT, "noise_results/mech_probe_results.json")))
MEAS_PROLONG = 3.2

fig, (ax, bx) = plt.subplots(1, 2, figsize=(9.6, 4.2))

# ---- (a) multiplicative sweep -------------------------------------------- #
mults = [m for m in R["mults"] if m > 0]
cms = [R["c_m_values"][str(m)] for m in mults]
for N in NS:
    y = [R["matrix"][f"mult{m}_N{N}"]["prolongation"] for m in mults]
    ax.plot(cms, y, color=colors[N], marker=markers[N], label=f"$N$ = {N}")
ax.axhline(1.0, color=P["black"], ls="--", lw=0.9)
ax.axhline(MEAS_PROLONG, color=P["purple"], ls=":", lw=1.2)
cphys = R["c_m_phys"]
per_n_cm = [R["physical_estimate"]["per_N"][str(N)]["c_m_median"] for N in NS]
ax.axvspan(min(per_n_cm), max(per_n_cm), color=P["grey"], alpha=0.18, lw=0)
ax.axvline(cphys, color=P["grey"], ls="-", lw=0.9)
ax.set_xscale("log")
ax.annotate(rf"physical $c_m\approx{cphys:.3f}$" + "\n(band: per-$N$ spread)",
            xy=(cphys, 0.62), xytext=(cphys * 1.25, 0.52), fontsize=7,
            color=P["grey"])
ax.annotate(rf"measured target $\approx{MEAS_PROLONG:g}$",
            xy=(0.16, MEAS_PROLONG), xytext=(0.16, MEAS_PROLONG - 0.22),
            fontsize=7, color=P["purple"])
ax.set_xlabel(r"multiplicative amplitude $c_m$   ($\sigma_{\mathrm{eff}}"
              r"= c_m\,[1-\min(\rho_1,\rho_2)]$)")
ax.set_ylabel("prolongation factor (median capture / deterministic)")
ax.set_title("(a) breath-locked multiplicative noise ($N$-free amplitude)")
ax.set_ylim(0.2, 3.45)
ax.legend(loc="center left", bbox_to_anchor=(0.02, 0.72), fontsize=7)

# ---- (b) WS-constants decomposition -------------------------------------- #
red_med = {}
for N in NS:
    p = os.path.join(ROOT, f"manifold_results/ws_decomposition_N{N}.jsonl")
    rows = [json.loads(l) for l in open(p) if l.strip()]
    act = {r["seed"]: r for r in rows if r["kind"] == "actual"}
    par = {r["seed"]: r for r in rows if r["kind"] == "partner"}
    ta = np.array([act[s]["t_abs"] for s in sorted(act)])
    tp = np.array([par[s]["t_abs"] for s in sorted(act)])
    tr = np.nanmedian([act[s]["t_capture_reduced"] for s in sorted(act)])
    red_med[N] = float(tr)
    bx.plot(ta, tp, markers[N], color=colors[N], ms=2.6, alpha=0.35, mew=0,
            label=f"$N$ = {N}")
    bx.plot(np.median(ta), np.median(tp), markers[N], color=colors[N], ms=9,
            mec=P["black"], mew=1.0, zorder=5)
lims = [8, 2300]
bx.plot(lims, lims, color=P["black"], ls="--", lw=0.9, zorder=1)
rm = np.median(list(red_med.values()))
bx.axhline(rm, color=P["purple"], ls=":", lw=1.2)
bx.annotate("reduced-flow median\n(partner would sit here if constants\n"
            "carried the prolongation)", xy=(lims[0] * 1.2, rm), fontsize=7,
            color=P["purple"], va="bottom",
            xytext=(lims[0] * 1.2, rm * 1.08))
bx.annotate(r"$t_{\mathrm{partner}}=t_{\mathrm{actual}}$",
            xy=(700, 700), xytext=(550, 1000), fontsize=7, color=P["black"],
            rotation=0)
bx.set_xscale("log")
bx.set_yscale("log")
bx.set_xlim(*lims)
bx.set_ylim(*lims)
bx.set_xlabel(r"actual-constants absorption time $t_{\mathrm{actual}}$ (s)")
bx.set_ylabel(r"uniformized-constants $t_{\mathrm{partner}}$ (s)")
bx.set_title("(b) WS-constants decomposition (large = per-$N$ median)")
bx.legend(loc="lower right", fontsize=7)

fig.tight_layout()
paths = style.savefig(fig, "fig_mech_probe", tight=True)
print("wrote", *paths)
