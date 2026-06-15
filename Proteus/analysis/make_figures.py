"""Regenerate the revision figures from committed analysis JSON.

  Fig 4 (left), enlarged baseline  -> figures/fig4_left_enlarged.{pdf,png}
        Exposed-cleft rate | triad for the PET-hydrolase branch, the AChE branch,
        and the ENLARGED random floor (Wilson 95% CI error bar). Annotated with the
        non-superiority result from analysis/equivalence_tost.py.

  Supplementary pLDDT figure          -> figures/figS_plddt_bands.{pdf,png}
        ECDF + box of global pLDDT by band (divergent <25%, divergent <20%,
        random baseline triad-bearers, near-homolog above-line), from
        analysis/plddt_confound.py. Shows the divergent tail is NOT lower-pLDDT
        than the baseline, so the exposed-cleft depletion is not a prediction-error
        artifact.

Reproduce, from the repo root:
    PYTHONPATH=src python analysis/make_figures.py
"""
from __future__ import annotations

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
WT = os.path.dirname(HERE)
PROC = os.path.join(WT, "data", "processed")
FIGS = os.path.join(WT, "figures")
os.makedirs(FIGS, exist_ok=True)

TEAL = "#2a9d8f"
RED = "#e76f51"
ORANGE = "#e9a23b"
GREY = "#6c757d"


def fig4_left():
    eq = json.load(open(os.path.join(PROC, "equivalence.json")))
    pet = eq["pet_branch"]["rate"] * 100
    ache = eq["ache_branch"]["rate"] * 100
    base = eq["baseline"]["rate"] * 100
    bw = [w * 100 for w in eq["baseline"]["wilson95"]]
    n2 = eq["baseline"]["n"]
    ub = eq["one_sided_95_upper_bound_on_enrichment"] * 100

    fig, ax = plt.subplots(figsize=(4.2, 4.0))
    xs = [0, 1, 2]
    bars = ax.bar(xs, [pet, ache, base], width=0.62,
                  color=[TEAL, RED, ORANGE], zorder=3)
    # Wilson CI on the enlarged floor
    ax.errorbar(2, base, yerr=[[base - bw[0]], [bw[1] - base]], fmt="none",
                ecolor="black", capsize=5, lw=1.4, zorder=4)
    # floor reference line
    ax.axhline(base, ls="--", color=ORANGE, lw=1.2, zorder=2)
    for x, v in zip(xs, [pet, ache, base]):
        ax.text(x, v + 1.5, f"{v:.1f}%", ha="center", va="bottom", fontsize=11)
    ax.text(0.5, base + 2.0,
            f"not enriched vs floor\n(1-sided 95% upper bound\non enrichment = {ub:+.1f} pp)",
            ha="center", va="bottom", fontsize=8.5, style="italic", color=GREY)
    ax.set_xticks(xs)
    ax.set_xticklabels(["PETase\nbranch", "AChE\nbranch",
                        f"random floor\n(enlarged, n={n2})"])
    ax.set_ylabel("exposed-cleft rate | triad (%)")
    ax.set_ylim(0, 70)
    ax.set_title("Top-300 per anchor vs enlarged baseline", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"fig4_left_enlarged.{ext}"), dpi=200)
    plt.close(fig)
    print(f"[fig] fig4_left_enlarged: PET {pet:.1f} AChE {ache:.1f} "
          f"floor {base:.1f} [{bw[0]:.1f},{bw[1]:.1f}]")


def _ecdf(vals):
    v = np.sort(np.asarray(vals, float))
    y = np.arange(1, len(v) + 1) / len(v)
    return v, y


def figS_plddt():
    pc = json.load(open(os.path.join(PROC, "plddt_confound.json")))
    a = pc["_arrays"]
    groups = [
        ("divergent <25%", np.array(a["divergent_global"]) * 100, TEAL),
        ("divergent <20%", np.array(a["divergent_sub20_global"]) * 100, "#1d6f64"),
        ("random baseline", np.array(a["baseline_global"]) * 100, ORANGE),
        ("near-homolog\nabove-line", np.array(a["near_homolog_global"]) * 100, RED),
    ]
    groups = [(n, v, c) for n, v, c in groups if len(v)]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(8.6, 4.0))
    # left: ECDF
    for name, v, c in groups:
        x, y = _ecdf(v)
        axL.step(x, y, where="post", color=c, lw=1.8,
                 label=f"{name.replace(chr(10),' ')} (n={len(v)}, md={np.median(v):.1f})")
    axL.set_xlabel("global pLDDT (%)")
    axL.set_ylabel("ECDF")
    axL.set_title("pLDDT by band (ECDF)", fontsize=11)
    axL.legend(fontsize=7.5, loc="upper left")
    axL.spines[["top", "right"]].set_visible(False)

    # right: box/violin
    data = [v for _, v, _ in groups]
    labels = [n for n, _, _ in groups]
    cols = [c for _, _, c in groups]
    parts = axR.violinplot(data, showmedians=True, showextrema=False)
    for pc_, c in zip(parts["bodies"], cols):
        pc_.set_facecolor(c)
        pc_.set_alpha(0.55)
    axR.set_xticks(range(1, len(labels) + 1))
    axR.set_xticklabels(labels, fontsize=8)
    axR.set_ylabel("global pLDDT (%)")
    axR.set_title("pLDDT by band (distribution)", fontsize=11)
    axR.spines[["top", "right"]].set_visible(False)

    dvb = next(t for t in pc["tests_global"]
               if t["comparison"] == "divergent vs baseline")
    msg = (f"divergent vs baseline: MWU p={dvb['mannwhitney_p']:.2g}, "
           f"Cliff's d={dvb['cliffs_delta_a_vs_b']:+.2f}\n"
           f"divergent NOT lower than baseline -> depletion is not a pLDDT artifact")
    fig.suptitle(msg, fontsize=9, y=0.02, va="bottom")
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"figS_plddt_bands.{ext}"), dpi=200)
    plt.close(fig)
    print(f"[fig] figS_plddt_bands: groups={[(n.replace(chr(10),' '),len(v)) for n,v,_ in groups]}")


def fig3():
    """Fig 3: (left) the 53-fold triad enrichment (UNCHANGED literals), (right) the
    exposed-cleft rate among random triad-bearers, updated to the enlarged floor."""
    eq = json.load(open(os.path.join(PROC, "equivalence.json")))
    base = eq["baseline"]["rate"] * 100
    bw = [w * 100 for w in eq["baseline"]["wilson95"]]
    n2 = eq["baseline"]["n"]
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(7.2, 4.0))
    axL.bar([0, 1], [1.9, 99.3], width=0.6, color=[GREY, TEAL], zorder=3)
    for x, v in zip([0, 1], [1.9, 99.3]):
        axL.text(x, v + 1.5, f"{v}%", ha="center", fontsize=11)
    axL.set_xticks([0, 1]); axL.set_xticklabels(["random", "fold-class\nhits"])
    axL.set_ylabel("catalytic-triad rate (%)"); axL.set_ylim(0, 108)
    axL.set_title("S4 — fold-class search enriches (53×)", fontsize=10.5)
    axL.spines[["top", "right"]].set_visible(False)
    axR.bar([0], [base], width=0.5, color=ORANGE, zorder=3)
    axR.errorbar(0, base, yerr=[[base - bw[0]], [bw[1] - base]], fmt="none",
                 ecolor="black", capsize=5, lw=1.4)
    axR.text(0, base + 2, f"{base:.1f}%", ha="center", fontsize=11)
    axR.text(0, base / 2, "about half of random\nserine hydrolases pass",
             ha="center", va="center", fontsize=8.5, style="italic", color="white")
    axR.set_xticks([0]); axR.set_xticklabels([f"random\ntriad-bearers\n(n={n2})"])
    axR.set_ylabel("clear exposed-cleft line (%)"); axR.set_ylim(0, 100)
    axR.set_title("S5 — cleft line is not selective", fontsize=10.5)
    axR.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"fig3_enlarged.{ext}"), dpi=200)
    plt.close(fig)
    print(f"[fig] fig3_enlarged: floor bar {base:.1f}% [{bw[0]:.1f},{bw[1]:.1f}], n={n2}")


def fig4_right():
    """Fig 4 (right): bits-stratified tail gradient vs the ENLARGED floor line."""
    bg = json.load(open(os.path.join(PROC, "bits_gradient_enlarged.json")))
    fl = bg["enlarged_floor"]["rate"] * 100
    tiers = bg["tiers"]
    labels = [t["label"] for t in tiers]
    rates = [t["rate"] * 100 for t in tiers]
    fig, ax = plt.subplots(figsize=(4.4, 4.0))
    xs = list(range(len(tiers)))
    ax.plot(xs, rates, "-o", color=TEAL, lw=2, ms=7, zorder=3)
    ax.axhline(fl, ls="--", color=ORANGE, lw=1.3, zorder=2)
    ax.text(len(xs) - 1, fl + 1.5, f"enlarged floor {fl:.1f}%", ha="right",
            color=ORANGE, fontsize=8.5)
    for x, t in zip(xs, tiers):
        v = t["rate"] * 100
        tag = f"{v:.0f}%"
        if t["above_floor"]:
            tag += "*"
        ax.text(x, v + 2.5, tag, ha="center", fontsize=9.5)
    top = tiers[0]
    ax.text(0.02, 0.06,
            f"top-25 above floor: {100*top['rate']:.0f}% vs {fl:.1f}%\n"
            f"Fisher p={top['fisher_p_vs_enlarged_floor']:.3g}, "
            f"{top['rr_vs_enlarged_floor']:.2f}×",
            transform=ax.transAxes, fontsize=8, style="italic", color=GREY,
            va="bottom")
    ax.set_xticks(xs)
    ax.set_xticklabels([l.replace("rank ", "rank\n") for l in labels], fontsize=8)
    ax.set_ylabel("exposed-cleft rate | triad (%)")
    ax.set_ylim(0, 100)
    ax.set_title("Gradient crosses floor between top-50/top-100", fontsize=9.5)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"fig4_right_enlarged.{ext}"), dpi=200)
    plt.close(fig)
    print(f"[fig] fig4_right_enlarged: floor {fl:.1f}%, top-25 "
          f"p={top['fisher_p_vs_enlarged_floor']:.3g}")


def main() -> int:
    fig4_left()
    fig4_right()
    fig3()
    figS_plddt()
    print(f"[fig] wrote figures to {FIGS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
