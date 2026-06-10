"""
Cross-system label experiment figure (revision item CP3, \\label{fig:xlabel}).

Renders paper_figures/fig_xlabel.{pdf,png} from committed results only — no
integration, no relabeling, no randomness:

  tools/xlabel/results/xlabel_runs.jsonl          per-run t_fp / t_abs / t_camp
  paper/revision-data-gated/results_xlabel.json   per-cell aggregates + CIs

Content (the honest null: the graze/absorption gap does not transfer to the ring):
  (a) per-run first-passage label t_fp versus absorption-grade label t_abs for all
      900 runs in the three cells (log axes); virtually every point sits on the
      diagonal — the 3-10% that sit above it are the runs whose first sub-threshold
      excursion self-heals;
  (b) over-trigger fraction (first excursion recovers) per cell with bootstrap 95%
      CIs against the mean-field range (69-97% of detector firings; up to 98% of
      first crossings, Sec. 3);
  (c) exponential-MLE lifetime ratio tau_abs/tau_fp per cell with bootstrap 95%
      CIs against the mean-field factor ~2.2.

Run: python3 paper_figures/fig_xlabel.py
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

style.apply_style()
P = style.PALETTE

RUNS = os.path.join(ROOT, "tools", "xlabel", "results", "xlabel_runs.jsonl")
AGG = os.path.join(ROOT, "paper", "revision-data-gated", "results_xlabel.json")

CELL_ORDER = ["b0.110_N64", "b0.130_N64", "b0.130_N256"]
CELL_LABEL = {
    "b0.110_N64": r"$\beta=0.110$, $N=64$",
    "b0.130_N64": r"$\beta=0.130$, $N=64$",
    "b0.130_N256": r"$\beta=0.130$, $N=256$",
}
CELL_COLOR = {
    "b0.110_N64": P["blue"],
    "b0.130_N64": P["vermillion"],
    "b0.130_N256": P["green"],
}
CELL_MARKER = {"b0.110_N64": "o", "b0.130_N64": "s", "b0.130_N256": "^"}


def main():
    runs = {c: [] for c in CELL_ORDER}
    with open(RUNS) as f:
        for line in f:
            d = json.loads(line)
            runs[d["cell"]].append((d["t_fp"], d["t_abs"]))
    agg = json.load(open(AGG))["cells"]

    fig, (ax_a, ax_b, ax_c) = plt.subplots(
        1, 3, figsize=(7.0, 2.5), gridspec_kw={"width_ratios": [1.35, 1, 1]})

    # (a) per-run t_fp vs t_abs scatter, log-log, diagonal
    lo, hi = 20, 9000
    ax_a.plot([lo, hi], [lo, hi], color=P["grey"], lw=0.8, ls=":", zorder=1)
    for cell in CELL_ORDER:
        fp, ab = np.array(runs[cell]).T
        on = np.abs(ab - fp) < 1e-9
        ax_a.scatter(fp[on], ab[on], s=7, marker=CELL_MARKER[cell],
                     facecolors="none", edgecolors=CELL_COLOR[cell], lw=0.5,
                     alpha=0.45, zorder=2, label=CELL_LABEL[cell])
        ax_a.scatter(fp[~on], ab[~on], s=11, marker=CELL_MARKER[cell],
                     color=CELL_COLOR[cell], zorder=3)
    ax_a.set_xscale("log")
    ax_a.set_yscale("log")
    ax_a.set_xlim(lo, hi)
    ax_a.set_ylim(lo, hi)
    ax_a.set_xlabel(r"first-passage label $t_{\mathrm{fp}}$ (t.u.)")
    ax_a.set_ylabel(r"absorption label $t_{\mathrm{abs}}$ (t.u.)")
    ax_a.set_title("(a) the two labels coincide", fontsize=9)
    ax_a.legend(fontsize=5.8, loc="upper left", handletextpad=0.1,
                borderaxespad=0.2)

    x = np.arange(len(CELL_ORDER))

    # (b) over-trigger fraction per cell vs the mean-field range
    ax_b.axhspan(0.69, 0.97, color=P["grey"], alpha=0.25, lw=0)
    ax_b.text(0.97, 0.83, "mean field\n(69-97% of firings)",
              transform=ax_b.get_yaxis_transform(), fontsize=6,
              color=P["black"], ha="right", va="center")
    for i, cell in enumerate(CELL_ORDER):
        a = agg[cell]
        ci = a["over_trigger_ci95"]
        f = a["over_trigger_fraction"]
        ax_b.errorbar([i], [f], yerr=[[f - ci[0]], [ci[1] - f]],
                      fmt=CELL_MARKER[cell], color=CELL_COLOR[cell],
                      ms=5, capsize=3, lw=1)
    ax_b.set_xticks(x)
    ax_b.set_xticklabels([CELL_LABEL[c] for c in CELL_ORDER],
                         fontsize=5.8, rotation=20)
    ax_b.set_xlim(-0.6, 2.6)
    ax_b.set_ylim(0, 1.02)
    ax_b.set_ylabel("over-trigger fraction")
    ax_b.set_title("(b) first crossings absorb", fontsize=9)

    # (c) MLE lifetime ratio per cell vs the mean-field factor
    ax_c.axhline(1.0, color=P["grey"], lw=0.8, ls=":")
    ax_c.axhline(2.2, color=P["grey"], lw=0.9, ls="--")
    ax_c.text(0.97, 2.13, r"mean field $\approx 2.2\times$",
              transform=ax_c.get_yaxis_transform(), fontsize=6,
              color=P["black"], ha="right", va="top")
    for i, cell in enumerate(CELL_ORDER):
        a = agg[cell]
        ci = a["ratio_abs_fp_mle_ci95"]
        rt = a["ratio_abs_fp_mle"]
        ax_c.errorbar([i], [rt], yerr=[[rt - ci[0]], [ci[1] - rt]],
                      fmt=CELL_MARKER[cell], color=CELL_COLOR[cell],
                      ms=5, capsize=3, lw=1)
    ax_c.set_xticks(x)
    ax_c.set_xticklabels([CELL_LABEL[c] for c in CELL_ORDER],
                         fontsize=5.8, rotation=20)
    ax_c.set_xlim(-0.6, 2.6)
    ax_c.set_ylim(0.95, 2.4)
    ax_c.set_yticks([1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2])
    ax_c.set_ylabel(r"$\hat\tau_{\mathrm{abs}}\,/\,\hat\tau_{\mathrm{fp}}$")
    ax_c.set_title("(c) lifetimes move 2-4%", fontsize=9)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        out = os.path.join(HERE, f"fig_xlabel.{ext}")
        fig.savefig(out)
        print("wrote", out)

    # console check against the snippet's quoted numbers
    for cell in CELL_ORDER:
        a = agg[cell]
        print(f"{cell}: over-trigger {a['over_trigger_fraction']:.3f} "
              f"{a['over_trigger_ci95']}  ratio {a['ratio_abs_fp_mle']:.4f} "
              f"{a['ratio_abs_fp_mle_ci95']}")


if __name__ == "__main__":
    main()
