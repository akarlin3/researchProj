#!/usr/bin/env python
"""Build the HC2/CS2 robustness-battery figure (3 panels) from the result JSONs.

(A) in-silico misspecification isolation: railed fraction vs SNR under the
    well-specified bi-exponential and three misspecifications -- they overlap, and
    sit below the real-data 54.7% which the well-specified hard-corner recovers.
(B) phantom / known-truth f-sweep: railed fraction and relative D* recovery error
    rise together as perfusion fraction f falls.
(C) model criticism on real abdomen: criticised fraction (null baseline, all
    voxels, railed, non-railed) -- railed voxels are LESS criticised.

Run:  python Sextant/scripts/make_robustness_figure.py
"""
from __future__ import annotations

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_R = os.path.join(_ROOT, "results")
_OUT = os.path.join(_ROOT, "paper", "figures")


def main():
    misspec = json.load(open(os.path.join(_R, "misspecification_isolation.json")))
    phantom = json.load(open(os.path.join(_R, "phantom_recovery.json")))
    crit = json.load(open(os.path.join(_R, "model_criticism.json")))

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(13.5, 4.0))

    # ---- Panel A: misspecification isolation ----
    forwards = misspec["meta"]["forwards"]
    snrs = misspec["meta"]["snrs"]
    cells = {(c["forward"], c["snr"]): c for c in misspec["cells"]}
    labels = {"biexp_WS": "well-specified (correct model)", "triexp": "tri-exp tail",
              "kurtosis": "kurtosis", "noise_chi": "non-Rician noise"}
    markers = {"biexp_WS": "o", "triexp": "s", "kurtosis": "^", "noise_chi": "D"}
    for fw in forwards:
        y = [cells[(fw, s)]["frac_railed"] * 100 for s in snrs]
        axA.plot(snrs, y, marker=markers[fw], lw=(2.4 if fw == "biexp_WS" else 1.3),
                 label=labels[fw], color=("C3" if fw == "biexp_WS" else None),
                 zorder=(3 if fw == "biexp_WS" else 2))
    corner10 = cells[("biexp_WS", 10.0)]["frac_railed_hard_corner"] * 100
    axA.scatter([10], [corner10], s=120, marker="*", color="C3", zorder=5,
                label="well-specified, hard corner")
    axA.axhline(54.7, ls="--", color="0.4", lw=1.2)
    axA.text(60, 55.6, "real abdomen 54.7%", fontsize=8, color="0.3")
    axA.set_xscale("log")
    axA.set_xticks(snrs); axA.set_xticklabels([str(int(s)) for s in snrs])
    axA.set_xlabel("SNR"); axA.set_ylabel("NLLS $D^*$ railed (%)")
    axA.set_title("(A) Railing is estimator pathology,\nnot simulator–reality mismatch", fontsize=10)
    axA.legend(fontsize=7, loc="upper right"); axA.set_ylim(0, 62)

    # ---- Panel B: phantom f-sweep ----
    rows = phantom["B2_f_sweep"]["rows"]
    fvals = [r["f"] for r in rows]
    rail = [r["frac_railed"] * 100 for r in rows]
    err = [r["median_rel_dstar_error"] for r in rows]
    axB.plot(fvals, rail, "o-", color="C0", label="railed (%)")
    axB.set_xlabel("true perfusion fraction $f$")
    axB.set_ylabel("NLLS $D^*$ railed (%)", color="C0")
    axB.tick_params(axis="y", labelcolor="C0")
    axB.invert_xaxis()
    axB2 = axB.twinx()
    axB2.plot(fvals, err, "s--", color="C1", label="rel. $D^*$ error")
    axB2.set_ylabel("median rel. $D^*$ recovery error", color="C1")
    axB2.tick_params(axis="y", labelcolor="C1")
    axB.set_title("(B) Known-truth: railing tracks\n$D^*$ unrecoverability ($f$ decreasing →)", fontsize=10)

    # ---- Panel C: model criticism ----
    real = crit["real"]
    null = crit["validation"]["biexp_WS"]["flagged_frac"] * 100
    bars = {"well-spec.\nnull": null,
            "real\n(all)": real["criticised_frac"] * 100,
            "real |\nrailed": real["criticised_frac_among_railed"] * 100,
            "real |\nnon-railed": real["criticised_frac_among_nonrailed"] * 100}
    colors = ["0.6", "C0", "C2", "C3"]
    axC.bar(range(len(bars)), list(bars.values()), color=colors)
    axC.axhline(null, ls=":", color="0.5", lw=1)
    axC.set_xticks(range(len(bars))); axC.set_xticklabels(list(bars.keys()), fontsize=8)
    axC.set_ylabel("criticised / misspecified (%)")
    axC.set_title("(C) Real data: railed voxels are\nLESS model-misspecified", fontsize=10)
    for i, v in enumerate(bars.values()):
        axC.text(i, v + 0.3, f"{v:.1f}", ha="center", fontsize=8)

    fig.tight_layout()
    os.makedirs(_OUT, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(_OUT, f"robustness_battery.{ext}"), dpi=150,
                    bbox_inches="tight")
    print(f"[fig] wrote {_OUT}/robustness_battery.pdf/.png")


if __name__ == "__main__":
    main()
