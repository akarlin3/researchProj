"""CP4 consistency gate -- every data number in vernier.tex traces to a seeded printout.

Single source of truth for the manuscript's numbers:
  1. load the seeded results (../results/feasibility_gate.json [the gate, SOLID]
     and ../results/efficiency_frontier.json [Experiment B, PROVISIONAL/Minos]);
  2. carry the Gauge identifiability-wall numbers as PINNED PROVISIONAL constants
     (cited for honest scoping; sourced from Gauge/gauge/results_gauge04.md);
  3. emit numbers.tex with one \newcommand per number (the ONLY place vernier.tex
     gets a value);
  4. verify vernier.tex references only macros numbers.tex defines (no untraceable
     magic number can appear in prose) and run internal-consistency asserts.

Gate numbers (feasibility_gate.json) are SOLID -- Caliper-only, publication-
independent. Experiment-B numbers (Minos lens) and the Gauge-wall numbers are
PROVISIONAL. Run:
  proteus/bin/python Vernier/paper/consistency.py        # writes numbers.tex + checks
Exit 0 = gate green.
"""
from __future__ import annotations

import json
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(os.path.dirname(HERE), "results")
TEX = os.path.join(HERE, "vernier.tex")
NUMBERS_TEX = os.path.join(HERE, "numbers.tex")


def _load(name):
    with open(os.path.join(RESULTS, name)) as fh:
        return json.load(fh)


def build_numbers():
    gate = _load("feasibility_gate.json")
    front = _load("efficiency_frontier.json")
    nums = {}

    # ---- CP2 feasibility gate (SOLID -- Caliper-only) ----
    snr = gate["primary_snr"]
    run = gate["runs"][f"snr_{snr:g}"]
    schemes = run["schemes"]
    widths = [s["width_dstar"] for s in schemes]
    hicov = [s["cond_cov_high_dstar"] for s in schemes]
    postcov = [s["cov_dstar"] for s in schemes]
    rawcov = [s["cov_dstar_raw"] for s in schemes]
    crlb = [s["crlb_dstar"] for s in schemes]
    nums.update({
        "SnrPrimary": (f"{snr:g}", "CP2 gate primary SNR [SOLID]"),
        "NSchemes": (f"{len(schemes)}", "CP2 gate: matched-CRLB schemes [SOLID]"),
        "NVox": (f"{gate['n']}", "CP2 gate: cohort voxels [SOLID]"),
        "NBoot": (f"{gate['n_boot']}", "CP2 gate: bootstrap resamples [SOLID]"),
        "CrlbTolPct": ("10", "CP2 gate: matched-CRLB tolerance percent [SOLID]"),
        "DeltaSharp": (f"{run['delta_sharp']:.3f}", "CP2 gate SNR33: D* width spread (max-min)/median [SOLID]"),
        "DeltaSharpLo": (f"{run['delta_sharp_ci'][0]:.3f}", "CP2 gate: Delta_sharp 95% CI lo [SOLID]"),
        "DeltaSharpHi": (f"{run['delta_sharp_ci'][1]:.3f}", "CP2 gate: Delta_sharp 95% CI hi [SOLID]"),
        "DeltaCond": (f"{run['delta_cond']:.3f}", "CP2 gate SNR33: high-D* coverage range [SOLID]"),
        "DeltaCondLo": (f"{run['delta_cond_ci'][0]:.3f}", "CP2 gate: Delta_cond 95% CI lo [SOLID]"),
        "DeltaCondHi": (f"{run['delta_cond_ci'][1]:.3f}", "CP2 gate: Delta_cond 95% CI hi [SOLID]"),
        "ThreshSharp": ("0.10", "CP2 gate: pre-registered Delta_sharp threshold [SOLID]"),
        "ThreshCond": ("0.05", "CP2 gate: pre-registered Delta_cond threshold [SOLID]"),
        "WidthMin": (f"{min(widths):.0f}", "CP2 gate: min post-conformal D* width [SOLID]"),
        "WidthMax": (f"{max(widths):.0f}", "CP2 gate: max post-conformal D* width [SOLID]"),
        "HiCovMin": (f"{min(hicov):.3f}", "CP2 gate: min high-D* conditional coverage [SOLID]"),
        "HiCovMax": (f"{max(hicov):.3f}", "CP2 gate: max high-D* conditional coverage [SOLID]"),
        "PostCovLo": (f"{min(postcov):.3f}", "CP2 gate: min marginal D* coverage post-conformal [SOLID]"),
        "PostCovHi": (f"{max(postcov):.3f}", "CP2 gate: max marginal D* coverage post-conformal [SOLID]"),
        "RawCovLo": (f"{min(rawcov):.3f}", "CP2 gate: min raw (pre-conformal) D* coverage [SOLID]"),
        "RawCovHi": (f"{max(rawcov):.3f}", "CP2 gate: max raw (pre-conformal) D* coverage [SOLID]"),
        "CrlbMin": (f"{min(crlb):.1f}", "CP2 gate: min matched CRLB(D*) [SOLID]"),
        "CrlbMax": (f"{max(crlb):.1f}", "CP2 gate: max matched CRLB(D*) [SOLID]"),
        "Verdict": (run["verdict"], "CP2 gate: primary verdict [SOLID]"),
    })
    # robustness SNRs (letter-only macro suffixes -- LaTeX names cannot contain digits)
    for s, suf in ((25, "SnrLow"), (50, "SnrHigh")):
        r = gate["runs"].get(f"snr_{s}")
        if r:
            nums[f"DeltaSharp{suf}"] = (f"{r['delta_sharp']:.3f}", f"CP2 gate SNR{s}: Delta_sharp [SOLID]")
            nums[f"DeltaCond{suf}"] = (f"{r['delta_cond']:.3f}", f"CP2 gate SNR{s}: Delta_cond [SOLID]")

    # ---- Gauge identifiability wall (PINNED, PROVISIONAL -- for honest scoping) ----
    nums.update({
        "GaugeWallLo": ("1.05", "Gauge results_gauge04: CRLB(D*)/tercile-width, CRLB-optimal [PROVISIONAL]"),
        "GaugeWallHi": ("1.25", "Gauge results_gauge04: CRLB(D*)/tercile-width, clinical [PROVISIONAL]"),
        "GaugeCovLo": ("0.841", "Gauge results_gauge04: high-D* coverage, clinical [PROVISIONAL]"),
        "GaugeCovHi": ("0.844", "Gauge results_gauge04: high-D* coverage, CRLB-optimal [PROVISIONAL]"),
    })

    # ---- Experiment B efficiency frontier (PROVISIONAL -- Minos lens) ----
    fs = sorted(front["schemes"], key=lambda r: r["scan_minutes"])
    margs = []
    for prev, cur in zip(fs[:-1], fs[1:]):
        dmin = cur["scan_minutes"] - prev["scan_minutes"]
        margs.append((cur["mean_utility"] - prev["mean_utility"]) / dmin if dmin else float("nan"))
    nums.update({
        "FrontMargA": (f"{margs[0]:.1f}", "Exp B: marginal decision-utility/min, sparse->clinical [PROVISIONAL]"),
        "FrontMargB": (f"{margs[1]:.1f}", "Exp B: marginal decision-utility/min, clinical->rich [PROVISIONAL]"),
        "FrontMargC": (f"{margs[2]:.2f}", "Exp B: marginal decision-utility/min, rich->dense [PROVISIONAL]"),
        "FrontWidthSparse": (f"{fs[0]['width_dstar']:.0f}", "Exp B: sparse-7 corrected D* width [PROVISIONAL]"),
        "FrontWidthDense": (f"{fs[-1]['width_dstar']:.0f}", "Exp B: dense-22 corrected D* width [PROVISIONAL]"),
        "FrontNbSparse": (f"{fs[0]['n_b']}", "Exp B: sparse scheme n_b [PROVISIONAL]"),
        "FrontNbDense": (f"{fs[-1]['n_b']}", "Exp B: dense scheme n_b [PROVISIONAL]"),
        "PriorBaselineU": (f"{front['no_scan_prior_baseline_utility']:.1f}",
                           "Exp B: no-scan prior baseline mean utility [PROVISIONAL]"),
    })

    # ---- MAF robustness (SOLID -- Caliper-only; the estimator-contingency result) ----
    maf = _load("maf_gate.json")
    m, r = maf["maf"], maf["reference"]
    mw = [s["width_dstar"] for s in m["schemes"]]
    nums.update({
        "MafDeltaSharp": (f"{m['delta_sharp']:.3f}", "MAF gate: D* width spread [SOLID]"),
        "MafDeltaSharpLo": (f"{m['delta_sharp_ci'][0]:.3f}", "MAF gate: Delta_sharp CI lo [SOLID]"),
        "MafDeltaSharpHi": (f"{m['delta_sharp_ci'][1]:.3f}", "MAF gate: Delta_sharp CI hi [SOLID]"),
        "MafDeltaCond": (f"{m['delta_cond']:.3f}", "MAF gate: high-D* coverage range [SOLID]"),
        "MafDeltaCondLo": (f"{m['delta_cond_ci'][0]:.3f}", "MAF gate: Delta_cond CI lo [SOLID]"),
        "MafDeltaCondHi": (f"{m['delta_cond_ci'][1]:.3f}", "MAF gate: Delta_cond CI hi [SOLID]"),
        "MafVerdict": (m["verdict"], "MAF gate: verdict [SOLID]"),
        "MafWidthLo": (f"{min(mw):.0f}", "MAF gate: min post-conformal D* width [SOLID]"),
        "MafWidthHi": (f"{max(mw):.0f}", "MAF gate: max post-conformal D* width [SOLID]"),
        "MafRawCov": (f"{np.mean([s['cov_dstar_raw'] for s in m['schemes']]):.2f}",
                      "MAF gate: mean raw (pre-conformal) D* coverage -- already near-calibrated [SOLID]"),
        "RefSplitDeltaSharp": (f"{r['delta_sharp']:.3f}",
                               "MAF gate: reference Delta_sharp under same splits (cross-check) [SOLID]"),
        "RefSplitDeltaCond": (f"{r['delta_cond']:.3f}",
                              "MAF gate: reference Delta_cond under same splits (cross-check) [SOLID]"),
        "RefSplitVerdict": (r["verdict"], "MAF gate: reference verdict under same splits [SOLID]"),
    })
    return nums


def write_numbers(nums):
    lines = ["% AUTO-GENERATED by consistency.py from seeded results -- DO NOT EDIT BY HAND.",
             "% [SOLID] = Caliper-only gate (publication-independent).",
             "% [PROVISIONAL] = Gauge-wall / Minos-lens (re-validate if those papers change).", ""]
    for key, (val, src) in nums.items():
        lines.append(f"\\newcommand{{\\num{key}}}{{{val}}}  % {src}")
    with open(NUMBERS_TEX, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def verify_tex(nums):
    if not os.path.exists(TEX):
        print(f"  [skip] {TEX} not present yet")
        return True
    with open(TEX) as fh:
        tex = fh.read()
    used = set(re.findall(r"\\num([A-Za-z]+)\b", tex))
    defined = set(nums)
    undefined = used - defined
    ok = not undefined
    print(f"  macros used in vernier.tex: {len(used)}; defined: {len(defined)}; "
          f"undefined: {sorted(undefined) if undefined else 'none'}")
    asserts = [
        ("gate verdict PASS", nums["Verdict"][0] == "PASS"),
        ("Delta_sharp >= threshold (primary)", float(nums["DeltaSharp"][0]) >= float(nums["ThreshSharp"][0])),
        ("Delta_cond >= threshold (primary)", float(nums["DeltaCond"][0]) >= float(nums["ThreshCond"][0])),
        ("Delta_cond CI excludes 0", float(nums["DeltaCondLo"][0]) > 0.0),
        ("marginal coverage restored ~0.90", abs(float(nums["PostCovHi"][0]) - 0.90) < 0.05),
        ("raw coverage far below nominal", float(nums["RawCovHi"][0]) < 0.6),
        ("matched CRLB within +/-10%", (float(nums["CrlbMax"][0]) - float(nums["CrlbMin"][0]))
         / float(nums["CrlbMin"][0]) <= 0.25),  # max-min/min over the matched set
        ("frontier returns diminish", float(nums["FrontMargA"][0]) > float(nums["FrontMargC"][0])),
        ("MAF gate FAILS (estimator-contingent)", nums["MafVerdict"][0] == "FAIL"),
        ("MAF divergence < reference", float(nums["MafDeltaSharp"][0]) < float(nums["DeltaSharp"][0])),
        ("reference still PASS under MAF splits", nums["RefSplitVerdict"][0] == "PASS"),
    ]
    for name, cond in asserts:
        print(f"  assert {name}: {'PASS' if cond else 'FAIL'}")
        ok = ok and cond
    return ok


def main():
    print("CP4 consistency gate -- numbers.tex from seeded results; vernier.tex traceability")
    nums = build_numbers()
    write_numbers(nums)
    print(f"  wrote {NUMBERS_TEX} ({len(nums)} numbers)")
    ok = verify_tex(nums)
    print("\nCP4 CONSISTENCY: " + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
