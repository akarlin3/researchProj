#!/usr/bin/env python
"""Generate the Matrix manuscript figures from seeded results (run-then-plot).

Two figures, committed as PDFs (regenerable):

  fig_synthetic.pdf  (a) the synthetic twin + the trust gate firing exactly on the
                          measured low-SNR zone (why AUROC = 1.00 is by construction);
                      (b) the closed-loop trajectory (mean true f, n_treat, mean dose).
  fig_ferry.pdf      (a) REAL anatomy (RTSTRUCT labels) + REAL dose geometry (RTDOSE);
                      (b) the F1 honest negative: held untrusted f-drop is 0 until real
                          dose geometry enters, then is positive (CI excludes 0).

Everything is synthetic/seeded except the Ferry anatomy+dose panels, which use the cached
public TCIA grids (or are skipped with a clear note if the substrate is unavailable). No
clinical claim; no patient blobs are committed. Run: <proteus python> paper/figures/make_figures.py
"""
from __future__ import annotations

import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))          # Matrix/
sys.path.insert(0, ROOT)

from matrix import MatrixConfig, Twin, Interfaces, run_loop, LoopState, TUMOR, TREAT
from matrix.fit import fit_scan
from matrix.loop import stage_scan, stage_posterior, stage_gates
from matrix.evaluate import convergence_series

LABEL_CMAP = ListedColormap(["#eaeaea", "#d1495b", "#3a6ea5"])  # normal / tumour / OAR
RESULTS = os.path.join(ROOT, "results")


def _single_state(cfg, ifaces):
    twin = Twin.build(cfg)
    rng = np.random.default_rng(cfg.seed + 12345)
    st = LoopState(iteration=0, n_voxels=cfg.n_voxels)
    stage_scan(twin, cfg, rng, st)
    stage_posterior(cfg, st, ifaces.ruler)
    stage_gates(cfg, st, ifaces.trust_gate, ifaces.action_gate)
    return twin, st


def fig_synthetic():
    cfg = MatrixConfig()
    ifaces = Interfaces.placeholders()
    twin, st = _single_state(cfg, ifaces)
    twin_run, states = run_loop(cfg, Interfaces.placeholders())
    series = convergence_series(states)

    labels = twin.labels.reshape(cfg.ny, cfg.nx)
    lowsnr = twin.lowsnr.reshape(cfg.ny, cfg.nx)
    untrust = (~st.trustworthy).reshape(cfg.ny, cfg.nx)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.2, 3.7))

    axA.imshow(labels, cmap=LABEL_CMAP, vmin=0, vmax=2, interpolation="nearest")
    # outline the measured low-SNR (untrustworthy) zone the gate keys on
    axA.contour(lowsnr.astype(float), levels=[0.5], colors="k", linewidths=1.2, linestyles="--")
    ys, xs = np.where(untrust)
    axA.scatter(xs, ys, marker="x", s=42, c="k", linewidths=1.4,
                label="trust gate: untrustworthy")
    axA.set_title("(a) twin + trust gate\n(x = flagged; dashed = low-SNR zone)", fontsize=9)
    axA.set_xticks([]); axA.set_yticks([])
    handles = [plt.Line2D([], [], marker="s", ls="", mfc="#d1495b", mec="none", label="tumour"),
               plt.Line2D([], [], marker="s", ls="", mfc="#3a6ea5", mec="none", label="OAR"),
               plt.Line2D([], [], marker="x", ls="", c="k", label="flagged untrustworthy")]
    axA.legend(handles=handles, loc="lower right", fontsize=6, framealpha=0.9)

    it = np.arange(len(series["mean_f_truth"]))
    axB.plot(it, series["mean_f_truth"], "-o", color="#d1495b", label="mean true $f$")
    axB.set_xlabel("iteration"); axB.set_ylabel("mean true $f$", color="#d1495b")
    axB.tick_params(axis="y", labelcolor="#d1495b")
    axB.set_title("(b) closed-loop trajectory", fontsize=9)
    ax2 = axB.twinx()
    ax2.plot(it, series["n_treat"], "-s", color="#3a6ea5", label="$n_{\\mathrm{treat}}$")
    ax2.plot(it, series["mean_dose"], "-^", color="#666", label="mean dose (Gy)")
    ax2.set_ylabel("$n_{\\mathrm{treat}}$ / mean dose (Gy)")
    lines = axB.get_lines() + ax2.get_lines()
    axB.legend(lines, [l.get_label() for l in lines], loc="center right", fontsize=7)

    fig.tight_layout()
    out = os.path.join(HERE, "fig_synthetic.pdf")
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"  wrote {os.path.relpath(out, ROOT)}")


def fig_ferry():
    fer = json.load(open(os.path.join(RESULTS, "RESULTS_FERRY_CP2.json")))
    sy = fer["synthetic"]; gf = fer["grounded_full"]; f1 = fer["F1_held_drop_real_dose"]

    # try the cached public grids for the anatomy/dose panels (no network, no blobs committed)
    labels = dose = None
    try:
        from matrix.ferry import load_substrate, FerryDataUnavailable
        try:
            sub = load_substrate(G=fer["config"]["nx"], allow_download=False)
            labels = np.asarray(sub.labels, float)
            dose = np.asarray(sub.dose_gy, float)
        except FerryDataUnavailable:
            print("  [note] Ferry cache unavailable; drawing F1 panel only "
                  "(anatomy/dose panels need the cached public grids).")
    except Exception as e:
        print(f"  [note] could not load Ferry substrate ({e!r}); F1 panel only.")

    ncols = 3 if labels is not None else 1
    fig, axes = plt.subplots(1, ncols, figsize=(3.4 * ncols, 3.5))
    axes = np.atleast_1d(axes)

    if labels is not None:
        a0 = axes[0]
        a0.imshow(labels, cmap=LABEL_CMAP, vmin=0, vmax=2, interpolation="nearest")
        a0.set_title("(a) REAL anatomy\n(RTSTRUCT labels)", fontsize=9)
        a0.set_xticks([]); a0.set_yticks([])
        a1 = axes[1]
        im = a1.imshow(dose, cmap="magma", interpolation="nearest")
        a1.set_title("(b) REAL dose geometry\n(RTDOSE, rescaled)", fontsize=9)
        a1.set_xticks([]); a1.set_yticks([])
        fig.colorbar(im, ax=a1, fraction=0.046, pad=0.04)
        axF = axes[2]
        tag = "(c)"
    else:
        axF = axes[0]; tag = ""

    # F1 honest-negative bar: held untrusted f-drop, synthetic vs grounded:full
    syn_held = sy["held_untrusted"][0]
    full_held = f1["value_ci"][0]
    full_lo, full_hi = f1["value_ci"][1], f1["value_ci"][2]
    trusted = gf["drop_trusted"][0]
    xs = [0, 1, 2]
    vals = [syn_held, full_held, trusted]
    errs = [[0, full_held - full_lo, 0], [0, full_hi - full_held, 0]]
    colors = ["#bbbbbb", "#d1495b", "#3a6ea5"]
    axF.bar(xs, vals, yerr=errs, capsize=5, color=colors)
    axF.set_xticks(xs)
    axF.set_xticklabels(["held\n(synthetic\ndose)", "held\n(REAL\ndose)", "treated\n(REAL\ndose)"],
                        fontsize=7)
    axF.set_ylabel("per-voxel $f$ drop")
    axF.set_title(f"{tag} F1: action-suppression\n$\\neq$ outcome-protection", fontsize=9)
    axF.axhline(0, color="k", lw=0.6)
    axF.set_ylim(0, max(vals) * 1.35)
    axF.annotate("held, but still\ndrops on real dose", xy=(1, full_held),
                 xytext=(1.62, full_held * 0.55), fontsize=6.5, ha="center",
                 arrowprops=dict(arrowstyle="->", lw=0.8))

    fig.tight_layout()
    out = os.path.join(HERE, "fig_ferry.pdf")
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"  wrote {os.path.relpath(out, ROOT)}")


def main():
    print("Matrix figures (run-then-plot from seeded results)")
    fig_synthetic()
    fig_ferry()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
