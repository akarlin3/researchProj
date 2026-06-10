"""
Depth-dependence of the reduced-flow aging shape (revision item CP1 figure).

Renders paper_figures/fig_depth_aging.{pdf,png} from committed analysis output
only -- no integration, no fitting, no randomness:

  paper/revision-data-gated/results_depth.json   assembled 8-point table
      (both fitters' k and CIs, depth beyond homoclinic, band-relative depth),
      changepoint brackets, all-1400 sensitivity anchor
      (produced by tools/reduced-ode/depth_aging_fit.py from the raw capture
      times committed in paper/revision-data-gated/results_corner.json)

Content (Appendix C, "Generality of the post-homoclinic regime"):
  (a) censored-Weibull shape k versus raw depth A - A_hc(beta): profile-CI
      and bootstrap-CI error bars per point, the memoryless k = 1 line, the
      beta <= 0.10 onset bracket [0.0205, 0.0595] (a coverage gap -- the onset
      is bracketed, not resolved), the operating corner (star) and its
      all-1400 sensitivity refit (open diamond). The shallow beta = 0.18
      point ages at the same raw depth where beta <= 0.10 points do not.
  (b) the same k versus band-relative depth (A - A_hc)/(A_hc - A_H): a single
      threshold in [0.27, 0.35] band-widths separates all eight points,
      the quantitative form of the Takens-Bogdanov band-narrowing reading.

Run: python3 paper_figures/fig_depth_aging.py
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools", "paper-figures"))

import style  # noqa: E402  (sets the Agg backend)

style.apply_style()
import matplotlib.pyplot as plt  # noqa: E402

P = style.PALETTE
IN_JSON = os.path.join(ROOT, "paper", "revision-data-gated", "results_depth.json")

BETA_COLOR = {0.03: P["blue"], 0.05: P["vermillion"], 0.10: P["green"],
              0.16: P["purple"], 0.18: P["orange"]}


def draw_points(ax, rows, xkey, dx_boot):
    """One point per row: profile-CI bar (heavy) + bootstrap-CI bar (light,
    offset by dx_boot, multiplicative if the axis is log)."""
    logx = ax.get_xscale() == "log"
    for r in rows:
        x = r[xkey]
        col = BETA_COLOR[r["beta"]]
        is_corner = "OPERATING CORNER" in r["tag"]
        lo, hi = r["k_ci_primary"]
        ax.errorbar([x], [r["k_primary"]],
                    yerr=[[r["k_primary"] - lo], [hi - r["k_primary"]]],
                    fmt="*" if is_corner else "o",
                    ms=11 if is_corner else 4.5, color=col, lw=1.1, capsize=2,
                    zorder=4 if is_corner else 3)
        blo, bhi = r["k_ci_boot"]
        xb = x * dx_boot if logx else x + dx_boot
        ax.errorbar([xb], [r["k_boot"]],
                    yerr=[[r["k_boot"] - blo], [bhi - r["k_boot"]]],
                    fmt="none", ecolor=col, elinewidth=0.6, capsize=1.2,
                    alpha=0.55, zorder=2)


def main():
    with open(IN_JSON) as f:
        d = json.load(f)
    rows = d["table"]
    anchor = d["anchor_all1400"]
    br_depth = d["changepoint"]["raw_depth_beta_le_0p10"]["bracket"]
    br_rel = d["changepoint"]["rel_depth_all8"]["bracket"]

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(7.0, 2.9))

    # ----------------------------------------------------------- (a) raw depth
    ax_a.axhline(1.0, color=P["grey"], lw=0.8, ls="--")
    ax_a.text(0.98, 1.005, "memoryless $k=1$", transform=ax_a.get_yaxis_transform(),
              ha="right", va="bottom", fontsize=6.5, color=P["grey"])
    ax_a.axvspan(br_depth[0], br_depth[1], color=P["skyblue"], alpha=0.25, lw=0)
    ax_a.text(0.5 * (br_depth[0] + br_depth[1]), 1.605,
              "onset bracket\n($\\beta \\leq 0.10$;\nno coverage)",
              ha="center", va="top", fontsize=6.0, color=P["black"])
    draw_points(ax_a, rows, "depth_beyond_hc", dx_boot=0.0035)

    # all-1400 sensitivity refit of the operating corner (degenerate N=4 rows in)
    blo, bhi = anchor["k_ci_boot"]
    ax_a.errorbar([anchor["depth_beyond_hc"] - 0.0035], [anchor["k_boot"]],
                  yerr=[[anchor["k_boot"] - blo], [bhi - anchor["k_boot"]]],
                  fmt="D", mfc="none", ms=4, color=P["grey"], lw=0.8, capsize=1.5,
                  zorder=2)
    ax_a.annotate("corner, all-1400\nsensitivity", fontsize=6.0, color=P["grey"],
                  xy=(anchor["depth_beyond_hc"] - 0.0035, anchor["k_boot"]),
                  xytext=(0.102, 1.12), ha="left", va="top",
                  arrowprops=dict(arrowstyle="-", color=P["grey"], lw=0.6))

    b18 = next(r for r in rows if r["beta"] == 0.18)
    ax_a.annotate("$\\beta=0.18$: shallow\nbut aging ($W=0.015$)",
                  fontsize=6.0, color=BETA_COLOR[0.18],
                  xy=(b18["depth_beyond_hc"], b18["k_primary"]),
                  xytext=(0.002, 1.135), ha="left", va="center",
                  arrowprops=dict(arrowstyle="-", color=BETA_COLOR[0.18], lw=0.6))

    ax_a.set_xlabel(r"depth beyond homoclinic $A - A_{hc}(\beta)$")
    ax_a.set_ylabel(r"Weibull shape $\hat k$")
    ax_a.set_xlim(0.0, 0.195)
    ax_a.set_ylim(0.88, 1.62)
    ax_a.set_title("(a) aging vs. depth: onset only bracketed", fontsize=9)


    # -------------------------------------------------- (b) band-relative depth
    ax_b.axhline(1.0, color=P["grey"], lw=0.8, ls="--")
    ax_b.set_xscale("log")
    ax_b.axvspan(br_rel[0], br_rel[1], color=P["skyblue"], alpha=0.25, lw=0)
    ax_b.text(0.5 * (br_rel[0] + br_rel[1]), 0.905,
              "single threshold\nseparates all 8",
              ha="center", va="bottom", fontsize=6.0, color=P["black"])
    draw_points(ax_b, rows, "rel_depth", dx_boot=1.07)
    ax_b.annotate("$\\beta=0.18$", fontsize=6.0, color=BETA_COLOR[0.18],
                  xy=(b18["rel_depth"], b18["k_primary"]),
                  xytext=(b18["rel_depth"], 1.24), ha="center", va="top")
    ax_b.set_xlabel(r"band-relative depth $(A - A_{hc})\,/\,(A_{hc} - A_{H})$")
    ax_b.set_ylabel(r"Weibull shape $\hat k$")
    ax_b.set_ylim(0.88, 1.62)
    ax_b.set_title("(b) depth in breathing-band widths", fontsize=9)

    handles = [plt.Line2D([], [], marker="o", ls="none", ms=4.5,
                          color=BETA_COLOR[b], label=rf"$\beta={b:.2f}$")
               for b in (0.03, 0.05, 0.10, 0.16, 0.18)]
    handles.append(plt.Line2D([], [], marker="*", ls="none", ms=9,
                              color=BETA_COLOR[0.05], label="operating corner"))
    ax_b.legend(handles=handles, fontsize=6.0, loc="lower right", ncol=2,
                borderpad=0.4, handletextpad=0.3, columnspacing=0.8)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        out = os.path.join(HERE, f"fig_depth_aging.{ext}")
        fig.savefig(out)
        print("wrote", out)

    # console check against the JSON the caption will cite
    print("onset bracket (depth, beta<=0.10):", [round(v, 4) for v in br_depth])
    print("threshold bracket (rel depth, all 8):", [round(v, 3) for v in br_rel])


if __name__ == "__main__":
    main()
