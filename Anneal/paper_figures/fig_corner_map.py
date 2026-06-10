"""
Corner-generality map figure (revision item CP1, reviewer Major Concern).

Renders paper_figures/fig_corner_map.{pdf,png} from the committed data file
paper/revision-data-gated/results_corner.json (produced by
tools/reduced-ode/corner_map_data.py — see its docstring for the full method).
Pure redraw: no dynamics computed here, no randomness, no manual
post-processing. Style via tools/paper-figures/style.py like every other
paper figure.

Content:
  - the four regimes of the (β, A) plane (no chimera / stable chimera /
    breathing / post-homoclinic) bounded by the saddle-node (Eq. 17), Hopf
    (Eq. 18) and freshly traced homoclinic curves;
  - the unstable-spiral escape rate σ(β, A) = Re λ of the chimera fixed point,
    shaded over the post-homoclinic region (grid of 500+ samples, linearly
    interpolated for display, clipped at the homoclinic curve) with labelled
    contours;
  - the Takens–Bogdanov point and the two shipped operating corners;
  - the 8 aging test points annotated with their censored-Weibull shape
    k (primary fitter; both fitters agree to <0.01 — see the JSON).

Regeneration chain:
  python3 tools/reduced-ode/corner_map_data.py     # heavy compute -> JSON
  python3 paper_figures/fig_corner_map.py          # this script
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools", "paper-figures"))

import style  # noqa: E402  (sets the Agg backend)
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

DATA = os.path.join(ROOT, "paper", "revision-data-gated", "results_corner.json")


def main():
    style.apply_style()
    with open(DATA) as f:
        d = json.load(f)

    c = d["curves"]
    beta = np.array(c["beta_dense"])
    A_SN = np.array(c["A_SN_series"])
    A_H = np.array(c["A_H_series"])
    hcb = np.array(c["beta_hc"])
    hcA = np.array(c["A_hc"])
    A_hc_dense = np.interp(beta, hcb, hcA)
    tb = c["tb_point"]
    pal = style.PALETTE

    fig, ax = plt.subplots(figsize=(7.0, 5.2))

    # ---- region fills (sub-homoclinic regimes, fig7 colour language) ----
    ax.fill_between(beta, 0.0, A_SN, color=pal["grey"], alpha=0.16, lw=0, zorder=0)
    ax.fill_between(beta, A_SN, A_H, color=pal["green"], alpha=0.16, lw=0, zorder=0)
    ax.fill_between(beta, A_H, A_hc_dense, color=pal["skyblue"], alpha=0.20,
                    lw=0, zorder=0)

    # ---- sigma field over the post-homoclinic region ----
    pts = [p for p in d["sigma_field"]["points"] if p["ok"]]
    pb = np.array([p["beta"] for p in pts])
    pA = np.array([p["A"] for p in pts])
    ps = np.array([p["sigma"] for p in pts])
    # display grid: linear interpolation of the regular sigma samples, clipped
    # at the homoclinic curve (and at the sampled domain edges)
    from scipy.interpolate import griddata
    gb, gA = np.meshgrid(np.linspace(hcb.min(), hcb.max(), 300),
                         np.linspace(pA.min() - 0.004, 0.55, 300))
    gs = griddata((pb, pA), ps, (gb, gA), method="linear")
    gs_near = griddata((pb, pA), ps, (gb, gA), method="nearest")
    gs = np.where(np.isfinite(gs), gs, gs_near)  # fill hull edge, then clip
    mask = gA < np.interp(gb, beta, A_hc_dense)
    gs = np.ma.masked_where(mask | ~np.isfinite(gs), gs)
    levels = MaxNLocator(nbins=9).tick_values(ps.min(), ps.max())
    cf = ax.contourf(gb, gA, gs, levels=levels, cmap="YlOrRd", alpha=0.75,
                     zorder=1)
    cl = ax.contour(gb, gA, gs, levels=levels[1:-1], colors="k",
                    linewidths=0.45, alpha=0.45, zorder=2)
    lab_levels = [l for l in levels[1:-1] if 0.010 <= l <= 0.042][::2]
    ax.clabel(cl, levels=lab_levels, fontsize=6, fmt="%.3f", inline=True,
              inline_spacing=2)
    cbar = fig.colorbar(cf, ax=ax, pad=0.015, fraction=0.046)
    cbar.set_label(r"unstable-spiral rate  $\sigma$  (s$^{-1}$)")
    cbar.ax.tick_params(labelsize=7)

    # ---- curves ----
    ax.plot(beta, A_SN, color=style.ROLES["sn"], lw=1.8, zorder=4)
    snn = c["sn_numeric"]
    ax.plot([p["beta"] for p in snn], [p["A"] for p in snn], ls="none",
            marker="o", mfc="none", mec=style.ROLES["sn"], mew=1.1, ms=4.5,
            zorder=5)
    ax.plot(beta, A_H, color=style.ROLES["hopf"], lw=1.8, ls="--", zorder=4)
    hnn = c["hopf_numeric"]
    ax.plot([p["beta"] for p in hnn], [p["A"] for p in hnn], ls="none",
            marker="s", mfc="none", mec=style.ROLES["hopf"], mew=1.1, ms=4.5,
            zorder=5)
    ax.plot(hcb, hcA, color=style.ROLES["homoclinic"], lw=1.8, ls="-.", zorder=4)
    ax.plot(hcb[::2], hcA[::2], ls="none", marker="D", mfc="none",
            mec=style.ROLES["homoclinic"], mew=1.0, ms=3.4, zorder=5)

    # ---- Takens-Bogdanov ----
    ax.plot([tb["beta"]], [tb["A"]], marker="o", ms=8, mfc=style.ROLES["tb"],
            mec="white", mew=1.0, zorder=6)
    ax.annotate(f"Takens–Bogdanov\n({tb['beta']}, {tb['A']})",
                xy=(tb["beta"], tb["A"]),
                xytext=(tb["beta"] - 0.010, tb["A"] - 0.105), ha="center",
                fontsize=7.5,
                arrowprops=dict(arrowstyle="->", color=style.ROLES["tb"], lw=0.9))

    # ---- region labels ----
    ax.text(0.115, 0.05, "no chimera", ha="center", va="center", fontsize=8,
            color="#555555")
    ax.text(0.046, 0.175, "stable\nchimera", ha="center", va="center",
            fontsize=8, color=pal["green"])
    ax.text(0.125, 0.305, "breathing\nchimera", ha="center", va="center",
            fontsize=8, color=pal["blue"])
    ax.text(0.124, 0.538, "post-homoclinic (sync only)", ha="center",
            va="center", fontsize=8.5, color=pal["vermillion"],
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none",
                      alpha=0.7))

    # ---- shipped corners ----
    corner_styles = {
        "stable": dict(marker="*", ms=14, color=style.ROLES["corner_stable"]),
        "transient": dict(marker="X", ms=10, color=style.ROLES["corner_transient"]),
    }
    for corner in d["corners"]:
        is_tr = abs(corner["A"] - 0.5) < 1e-9
        st = corner_styles["transient" if is_tr else "stable"]
        ax.plot([corner["beta"]], [corner["A"]], ls="none", mec="white",
                mew=1.0, zorder=9, **st)
    ax.annotate("operating corner", xy=(0.0525, 0.4975),
                xytext=(0.0685, 0.4615), fontsize=7.5,
                color=corner_styles["transient"]["color"], ha="left",
                arrowprops=dict(arrowstyle="->",
                                color=corner_styles["transient"]["color"],
                                lw=0.9), zorder=8)
    ax.annotate("stable corner (0.05, 0.2)", xy=(0.05, 0.2),
                xytext=(0.066, 0.155), fontsize=7.5,
                color=corner_styles["stable"]["color"], ha="left",
                arrowprops=dict(arrowstyle="->",
                                color=corner_styles["stable"]["color"],
                                lw=0.9), zorder=8)

    # ---- aging test points with fitted k ----
    # per-point label offsets (data coords), tuned only for legibility
    off = {(0.03, 0.461): (0.0045, -0.016), (0.03, 0.50): (0.0035, 0.0165),
           (0.05, 0.430): (0.0045, -0.016), (0.05, 0.50): (0.0065, -0.0185),
           (0.10, 0.374): (0.0045, -0.016), (0.10, 0.50): (0.0045, -0.0155),
           (0.18, 0.339): (-0.0045, -0.018), (0.16, 0.50): (0.0045, -0.0155)}
    for p in d["aging_points"]:
        b, A = p["beta"], p["A"]
        ax.plot([b], [A], ls="none", marker="o", ms=6.5, mfc="white",
                mec="k", mew=1.2, zorder=7)
        dx, dy = off.get((b, A), (0.004, -0.018))
        ax.text(b + dx, A + dy, f"$k={p['k_primary']:.2f}$", fontsize=7,
                color="k", ha="left" if dx > 0 else "right", va="center",
                zorder=8,
                bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none",
                          alpha=0.75))

    ax.set_xlim(0.02, 0.24)
    ax.set_ylim(0.0, 0.56)
    ax.set_xlabel(r"phase lag  $\beta$")
    ax.set_ylabel(r"coupling asymmetry  $A$")

    handles = [
        Line2D([], [], color=style.ROLES["sn"], lw=1.8,
               label="saddle-node (Eq. 17)"),
        Line2D([], [], color=style.ROLES["hopf"], lw=1.8, ls="--",
               label="Hopf (Eq. 18)"),
        Line2D([], [], color=style.ROLES["homoclinic"], lw=1.8, ls="-.",
               label="homoclinic (bisection)"),
        Line2D([], [], ls="none", marker="o", mfc="white", mec="k", mew=1.2,
               ms=6.5, label="aging test point (Weibull $k$)"),
        Line2D([], [], ls="none", **corner_styles["stable"], mec="white",
               label="shipped corner $A=0.2$"),
        Line2D([], [], ls="none", **corner_styles["transient"], mec="white",
               label="shipped corner $A=0.5$"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=7, framealpha=0.94,
              frameon=True, edgecolor="none")

    paths = style.savefig(fig, "fig_corner_map")
    plt.close(fig)
    print("wrote:", *[os.path.relpath(p, ROOT) for p in paths])


if __name__ == "__main__":
    main()
