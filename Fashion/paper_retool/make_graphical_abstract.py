#!/usr/bin/env python
"""Graphical-abstract figure for the boundary-railing-first manuscript (NMR in Biomedicine).

Single legible panel at the NMR-in-Biomedicine graphical-abstract size
(50 mm wide x 60 mm high), drawn DIRECTLY from the same seeded run artifact that
feeds the numbers-gate -- nothing is hand-entered:

  * Sextant/results/railing_results.json   (the four cohorts of Table 1: rates + bootstrap CIs)

The four bars are exactly the four rows of the manuscript railing table, so the
graphical abstract can never drift from the table or from numbers.tex. Run:

    KMP_DUPLICATE_LIB_OK=TRUE conda run -n proteus python make_graphical_abstract.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent                       # repo root (../../)
RAIL = ROOT / "Sextant" / "results" / "railing_results.json"
OUT = HERE / "figures"
OUT.mkdir(parents=True, exist_ok=True)

MM = 1 / 25.4
C_PRIMARY = "#1b4965"      # OSIPI abdomen homogeneous ROI (the headline cohort)
C_GENERAL = "#5fa8d3"      # generalisation cohorts

# Cohort name in the JSON -> short label for the panel. Order = Table 1 order.
COHORTS = [
    ("abdomen_homogeneous", "Abdomen\n(homog. ROI)", C_PRIMARY),
    ("abdomen_full", "Abdomen\n(full ROI)", C_GENERAL),
    ("lihc_liver_4b_0_50_500_800", "Liver 4-$b$", C_GENERAL),
    ("lihc_liver_3b_50_400_800", "Liver 3-$b$", C_GENERAL),
]


def main():
    rail = json.loads(RAIL.read_text())
    by_name = {c["name"]: c for c in rail["cohorts"]}

    rows = []
    for name, label, color in COHORTS:
        c = by_name[name]
        pt = c["tight"]["frac_railed"] * 100.0
        lo = c["bootstrap_ci"]["lo"] * 100.0
        hi = c["bootstrap_ci"]["hi"] * 100.0
        rows.append((label, pt, lo, hi, color))

    plt.rcParams.update({"font.size": 5.0, "axes.linewidth": 0.5})
    fig, ax = plt.subplots(figsize=(50 * MM, 60 * MM), layout="constrained")

    y = list(range(len(rows)))[::-1]
    for yi, (label, pt, lo, hi, color) in zip(y, rows):
        ax.barh(yi, pt, color=color, height=0.62)
        ax.errorbar(pt, yi, xerr=[[pt - lo], [hi - pt]], fmt="none",
                    ecolor="black", capsize=1.3, lw=0.7)
        ax.text(1.5, yi, label, va="center", ha="left", fontsize=4.6,
                color="white", linespacing=0.9)
        # value label sits to the right of the CI whisker so it never overlaps the cap;
        # one-decimal precision matches Table 1 / numbers.tex exactly (54.7, 47.8, ...)
        ax.text(hi + 1.2, yi, f"{pt:.1f}%", va="center", ha="left", fontsize=5.0)

    ax.set_yticks([])
    ax.set_xlim(0, 95)
    ax.set_xticks([0, 25, 50, 75])
    ax.tick_params(axis="x", labelsize=4.6, length=2, pad=1)
    ax.set_xlabel(r"$D^{*}$ railing rate (\%)".replace(r"\%", "%"), fontsize=5.4)
    ax.set_title(r"$D^{*}$ rails across IVIM cohorts", fontsize=6.0, pad=3)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    fig.savefig(OUT / "graphical_abstract.pdf")          # exact 50x60 mm (no tight bbox)
    fig.savefig(OUT / "graphical_abstract.png", dpi=600)
    plt.close(fig)
    print(f"wrote graphical_abstract.{{pdf,png}} (50x60 mm) -> {OUT}")
    print("  bars (Table 1, Sextant railing_results.json):",
          ", ".join(f"{l.splitlines()[0]} {p:.1f}%" for l, p, *_ in rows))


if __name__ == "__main__":
    main()
