"""CP3 -- regenerate the benchmark figures (vector PDF, deterministic from seed).

Reuses the cached predictions and the arm constructors from benchmark.py /
conditional.py so every figure traces to the same numbers as the CP1/CP2
printouts. Writes into gauge/figures/.

Run:  python scripts/make_figures.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from gauge.baselines import build_predictions, PARAM_NAMES
from gauge.benchmark import evaluate, DISPLAY_SCALE
from gauge.conformal import interval_width, empirical_coverage
from gauge.conditional import _arms_for_param

FIG_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "gauge", "figures")
REP_ALPHA = 0.10
PARAM_LABEL = ("D", "D*", "f")


def fig_coverage_vs_nominal(R, arms, M):
    alphas = sorted(R["meta"]["alphas"])
    nominal = [1 - a for a in alphas]
    show = [("raw:MDN-DeepEnsemble", "raw MDN", "o", "-"),
            ("raw:Bayesian-MCMC", "raw Bayesian", "^", "-"),
            ("conformal:CQR-HGB", "conformal (CQR)", "s", "--"),
            ("conformalized:MDN-DeepEnsemble", "conformalized MDN", "D", "-.")]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for j, ax in enumerate(axes):
        ax.plot([0.6, 1.0], [0.6, 1.0], color="k", lw=1, ls=":",
                label="ideal (y=x)")
        for arm, lab, mk, ls in show:
            cov = [M[(arm, j, a)]["coverage"] for a in alphas]
            ax.plot(nominal, cov, marker=mk, ls=ls, label=lab)
        ax.set_title(f"{PARAM_LABEL[j]}: realized vs nominal coverage")
        ax.set_xlabel("nominal coverage (1 - alpha)")
        ax.set_ylabel("realized coverage")
        ax.set_xlim(0.65, 0.97)
        ax.set_ylim(0.4, 1.0)
        ax.grid(alpha=0.3)
        if j == 0:
            ax.legend(fontsize=8, loc="lower right")
    fig.suptitle("Model-based UQ under-covers (below diagonal); "
                 "conformal sits on the diagonal", y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "coverage_vs_nominal.pdf"),
                bbox_inches="tight")
    plt.close(fig)


def fig_sharpness_vs_snr(R, arms):
    a = REP_ALPHA
    snr_levels = sorted(set(int(s) for s in R["meta"]["snr_grid"]))
    test_snr = R["test_snr"]
    show = [("raw:MDN-DeepEnsemble", "raw MDN", "o"),
            ("conformal:CQR-HGB", "conformal (CQR)", "s"),
            ("conformalized:MDN-DeepEnsemble", "conformalized MDN", "D")]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for j, ax in enumerate(axes):
        for arm, lab, mk in show:
            lo, hi = arms[arm](j, a)
            w = interval_width(lo, hi) * DISPLAY_SCALE[j]
            med = [np.median(w[test_snr == s]) for s in snr_levels]
            ax.plot(snr_levels, med, marker=mk, label=lab)
        ax.set_title(f"{PARAM_LABEL[j]} interval width (alpha={a})")
        ax.set_xlabel("SNR (b=0)")
        units = "1e-3 mm^2/s" if j < 2 else "fraction"
        ax.set_ylabel(f"median width [{units}]")
        ax.set_xscale("log")
        ax.grid(alpha=0.3)
        if j == 0:
            ax.legend(fontsize=8)
    fig.suptitle("Conformalized MDN is the sharpest method with guaranteed "
                 "coverage", y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "sharpness_vs_snr.pdf"),
                bbox_inches="tight")
    plt.close(fig)


def fig_conditional_heatmap(R):
    a = REP_ALPHA
    snr_levels = sorted(set(int(s) for s in R["meta"]["snr_grid"]))
    test_snr = R["test_snr"]
    arms = _arms_for_param(R, 1, a)            # D* parameter
    dstar = R["test_true"][:, 1]
    edges = np.quantile(dstar, [1 / 3, 2 / 3])
    regime = np.digitize(dstar, edges)
    reg_names = ["lo D*", "mid D*", "hi D*"]
    panels = ["raw-MDN", "CQR (plain)", "CQR (Mondrian/SNR)", "conformalized-MDN"]

    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    for ax, name in zip(axes.ravel(), panels):
        lo, hi = arms[name]
        grid = np.full((3, len(snr_levels)), np.nan)
        for rg in range(3):
            for ci, s in enumerate(snr_levels):
                m = (regime == rg) & (test_snr == s)
                if m.sum():
                    grid[rg, ci] = empirical_coverage(lo[m], hi[m], dstar[m])
        im = ax.imshow(grid, vmin=0.70, vmax=1.0, cmap="RdYlGn", aspect="auto")
        # mark the nominal contour by annotating values
        for rg in range(3):
            for ci in range(len(snr_levels)):
                ax.text(ci, rg, f"{grid[rg, ci]:.2f}", ha="center", va="center",
                        fontsize=8)
        ax.set_xticks(range(len(snr_levels)))
        ax.set_xticklabels([f"SNR{s}" for s in snr_levels], fontsize=8)
        ax.set_yticks(range(3))
        ax.set_yticklabels(reg_names, fontsize=8)
        ax.set_title(f"{name}  (D*, nominal {1-a:.2f})", fontsize=10)
    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.8)
    cbar.set_label("conditional coverage (green=ok, red=under-covers)")
    fig.suptitle("Conditional coverage of D*: the hi-D* row under-covers across "
                 "all SNR (and resists SNR-Mondrian)", y=1.0)
    fig.savefig(os.path.join(FIG_DIR, "conditional_coverage_heatmap.pdf"),
                bbox_inches="tight")
    plt.close(fig)


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    R = build_predictions(force=os.environ.get("GAUGE_FORCE") == "1")
    arms, M = evaluate(R)
    fig_coverage_vs_nominal(R, arms, M)
    fig_sharpness_vs_snr(R, arms)
    fig_conditional_heatmap(R)
    print("figures written to", FIG_DIR)
    for f in sorted(os.listdir(FIG_DIR)):
        print("  ", f)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
