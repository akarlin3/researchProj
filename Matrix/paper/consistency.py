#!/usr/bin/env python
"""Traceability gate for the Matrix manuscript (CP-style number consistency).

Regenerates ``numbers.tex`` from the seeded gate printouts —

  * ``results/RESULTS_CP2.json``       ruler + trust/action gates self-test (synthetic twin)
  * ``results/RESULTS_CP4.json``       closed-loop run (synthetic twin)
  * ``results/RESULTS_FERRY_CP2.json`` closed-loop run on REAL anatomy + dose geometry

— then verifies that every ``\\num*`` macro used in ``matrix.tex`` is defined and that the
manuscript's load-bearing claims hold internally (the honest negative, the AUROC scoping, the
byte-identity anchor, every headline CI). Run: ``python paper/consistency.py`` (exit 0 = PASS).

Every number in the manuscript traces to one of these seeded files; none is typed by hand.
The synthetic anchors regenerate offline via ``reproduce.sh``; the Ferry anchor regenerates
from the public TCIA dataset (network) and is committed so this gate runs offline.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "results"
TEX = HERE / "matrix.tex"
NUMBERS = HERE / "numbers.tex"

# The byte-identity anchor: matrix/loop.py is unchanged between synthetic and Ferry runs.
LOOP_PY_BLOB_SHA1 = "4a34806ac4fa55c0ce5453b9864d37c67abfda92"


def _load(name):
    return json.loads((RESULTS / name).read_text())


def build_numbers():
    cp2 = _load("RESULTS_CP2.json")
    cp4 = _load("RESULTS_CP4.json")
    fer = _load("RESULTS_FERRY_CP2.json")

    c = cp4["config"]
    rul = cp2["ruler"]; tg = cp2["trust_gate"]; ag = cp2["action_gate"]
    sup = cp4["suppression"]; con = cp4["convergence"]; tr = cp4["trajectory"]
    sub = fer["substrate"]; gf = fer["grounded_full"]; sy = fer["synthetic"]
    f1 = fer["F1_held_drop_real_dose"]

    def arrow(xs, f="{:d}"):
        # one math span (text-mode-safe): "$a \rightarrow b \rightarrow ...$"
        return "$" + r" \rightarrow ".join(f.format(x) for x in xs) + "$"

    nums = {
        # ---- run configuration (synthetic twin) -------------------------------
        "Nx": (f"{c['nx']}", "synthetic grid nx"),
        "Ny": (f"{c['ny']}", "synthetic grid ny"),
        "Nvox": (f"{c['n_voxels']}", "synthetic voxel count"),
        "Seed": (f"{c['seed']}", "master seed"),
        "Niter": (f"{c['n_iter']}", "loop iterations"),
        "SNR": (f"{c['snr']:.0f}", "nominal SNR"),
        "SNRlow": (f"{c['snr_low']:.0f}", "low-SNR (untrustworthy) zone SNR"),
        "DoseMin": (f"{c['dose_min']:.0f}", "dose floor (Gy)"),
        "DoseMax": (f"{c['dose_max']:.0f}", "dose ceiling (Gy)"),

        # ---- ruler (Fashion-shaped placeholder) -------------------------------
        "CovF": (f"{rul['coverage_f_95']:.3f}", "95% empirical coverage of f"),
        "CovFlo": (f"{rul['coverage_f_95_ci'][1]:.3f}", "95% coverage of f, CI lo"),
        "CovFhi": (f"{rul['coverage_f_95_ci'][2]:.3f}", "95% coverage of f, CI hi"),
        "ECEf": (f"{rul['ece_f']:.3f}", "expected calibration error of f"),

        # ---- trust gate (Minos-shaped placeholder) ----------------------------
        "AUROC": (f"{tg['auroc_sigmaf_vs_lowsnr']:.2f}",
                  "AUROC(calibrated sigma_f vs low-SNR label) -- separability by construction"),
        "FireLow": (f"{tg['fire_rate_lowsnr'][0]:.2f}", "fire rate inside low-SNR zone"),
        "FireGood": (f"{tg['fire_rate_goodsnr'][0]:.2f}", "fire rate outside low-SNR zone"),
        "FireGoodHi": (f"{tg['fire_rate_goodsnr'][2]:.3f}", "fire rate outside, CI hi"),

        # ---- action suppression (whole synthetic run) -------------------------
        "SupUngated": (f"{sup['act_rate_untrust_ungated'][0]:.3f}", "untrust action rate, ungated"),
        "SupUngatedLo": (f"{sup['act_rate_untrust_ungated'][1]:.3f}", "untrust ungated CI lo"),
        "SupUngatedHi": (f"{sup['act_rate_untrust_ungated'][2]:.3f}", "untrust ungated CI hi"),
        "SupGated": (f"{sup['act_rate_untrust_gated'][0]:.3f}", "untrust action rate, gated"),
        "SupTrust": (f"{sup['act_rate_trust_gated'][0]:.3f}", "trust action rate, gated"),
        "SupTrustLo": (f"{sup['act_rate_trust_gated'][1]:.3f}", "trust gated CI lo"),
        "SupTrustHi": (f"{sup['act_rate_trust_gated'][2]:.3f}", "trust gated CI hi"),

        # ---- convergence (synthetic twin) -------------------------------------
        "DropTrust": (f"{con['drop_trusted_tumour'][0]:.3f}", "trusted tumour f drop"),
        "DropTrustLo": (f"{con['drop_trusted_tumour'][1]:.3f}", "trusted f drop CI lo"),
        "DropTrustHi": (f"{con['drop_trusted_tumour'][2]:.3f}", "trusted f drop CI hi"),
        "DropHeld": (f"{con['drop_untrusted_tumour_held'][0]:.3f}", "untrusted (held) tumour f drop"),
        "TreatTraj": (arrow(tr["n_treat"]), "TREAT-action trajectory"),
        "FtruthFirst": (f"{tr['mean_f_truth'][0]:.3f}", "mean true f, iteration 0"),
        "FtruthLast": (f"{tr['mean_f_truth'][-1]:.3f}", "mean true f, final iteration"),

        # ---- Ferry: real anatomy + dose geometry ------------------------------
        "FerryDataset": (sub["collection"], "Ferry dataset"),
        "FerryDOI": (sub["doi"], "Ferry dataset DOI"),
        "FerryLicense": (sub["license"], "Ferry dataset license"),
        "FerryPatient": (sub["patient"].replace("_", r"\_"), "Ferry patient id (TeX-escaped)"),
        "FerryGrid": (f"{fer['config']['nx']}", "Ferry grid edge"),
        "FerryDoseLo": (f"{sub['dose_gy_range'][0]:.2f}", "real delivered dose min (Gy)"),
        "FerryDoseHi": (f"{sub['dose_gy_range'][1]:.1f}", "real delivered dose max (Gy)"),
        "FerrySliceZ": (f"{sub['slice_z_mm']}", "axial slice z (mm)"),
        "GTumorVox": (f"{gf['n_tumor']}", "grounded tumour voxel count"),
        "STumorVox": (f"{sy['n_tumor']}", "synthetic tumour voxel count (matched grid)"),
        "SmDropTrust": (f"{sy['drop_trusted'][0]:.3f}", "synthetic matched-grid trusted f drop"),
        "SmSupUngated": (f"{sy['sup_ungated'][0]:.3f}", "synthetic matched-grid untrust ungated rate"),
        "GTrust": (f"{gf['n_trusted']}", "grounded trusted tumour voxels"),
        "GUntrust": (f"{gf['n_untrusted']}", "grounded untrusted tumour voxels"),
        "GDropTrust": (f"{gf['drop_trusted'][0]:.3f}", "grounded trusted tumour f drop"),
        "GDropTrustLo": (f"{gf['drop_trusted'][1]:.3f}", "grounded trusted f drop CI lo"),
        "GDropTrustHi": (f"{gf['drop_trusted'][2]:.3f}", "grounded trusted f drop CI hi"),
        "GSupUngated": (f"{gf['sup_ungated'][0]:.3f}", "grounded untrust action rate, ungated"),
        "GSupUngatedLo": (f"{gf['sup_ungated'][1]:.3f}", "grounded ungated CI lo"),
        "GSupUngatedHi": (f"{gf['sup_ungated'][2]:.3f}", "grounded ungated CI hi"),
        "GSupGated": (f"{gf['sup_gated'][0]:.3f}", "grounded untrust action rate, gated"),
        "GTreatTraj": (arrow(gf["n_treat"]), "grounded TREAT-action trajectory"),

        # ---- F1 -- the honest negative (centerpiece) --------------------------
        "FOneHeld": (f"{f1['value_ci'][0]:.3f}", "F1 held untrusted f drop under REAL dose"),
        "FOneHeldLo": (f"{f1['value_ci'][1]:.3f}", "F1 held drop CI lo"),
        "FOneHeldHi": (f"{f1['value_ci'][2]:.3f}", "F1 held drop CI hi"),
        "FOneSyn": (f"{f1['synthetic_baseline_ci'][0]:.3f}", "F1 synthetic-baseline held drop"),

        # ---- byte-identity anchor ---------------------------------------------
        "LoopBlob": (LOOP_PY_BLOB_SHA1[:12], "matrix/loop.py git blob (truncated)"),
    }
    return nums, cp2, cp4, fer


def write_numbers(nums):
    lines = [
        "% AUTO-GENERATED by paper/consistency.py from seeded results -- DO NOT EDIT BY HAND.",
        "% Synthetic-twin anchors (CP2/CP4) regenerate offline via reproduce.sh; the Ferry",
        "% anchor regenerates from the public TCIA dataset (network) and is committed.",
        "% Calibration/trust/convergence numbers are PROVISIONAL on Fashion + Minos (\\PROV);",
        "% 'the loop closes' and the Ferry geometry result are SOLID for what they claim.",
        "",
    ]
    for k, (v, why) in nums.items():
        lines.append(f"\\newcommand{{\\num{k}}}{{{v}}}  % {why}")
    NUMBERS.write_text("\n".join(lines) + "\n")


def verify(nums, cp2, cp4, fer):
    ok = True
    # 1. every \num* used in the tex is defined
    if TEX.exists():
        used = set(re.findall(r"\\num([A-Za-z]+)\b", TEX.read_text()))
        undefined = used - set(nums)
        print(f"  macros used={len(used)} defined={len(nums)} undefined={sorted(undefined)}")
        if undefined:
            ok = False
    else:
        print("  (matrix.tex absent -- numbers.tex generated; skipping macro cross-check)")

    # 2. load-bearing internal-consistency asserts
    tg = cp2["trust_gate"]; rul = cp2["ruler"]
    con = cp4["convergence"]; sup = cp4["suppression"]; tr = cp4["trajectory"]
    gf = fer["grounded_full"]; sy = fer["synthetic"]; f1 = fer["F1_held_drop_real_dose"]
    asserts = [
        # ruler + trust gate (synthetic)
        ("ruler f-coverage near nominal (>=0.85)", rul["coverage_f_95"] >= 0.85),
        ("ruler ECE_f small (<0.08)", rul["ece_f"] < 0.08),
        ("AUROC == 1.00 (separability by construction)", tg["auroc_sigmaf_vs_lowsnr"] == 1.0),
        ("trust gate fires inside low-SNR zone (>0.6)", tg["fire_rate_lowsnr"][0] > 0.6),
        ("trust gate quiet outside (point <0.2)", tg["fire_rate_goodsnr"][0] < 0.2),
        # action suppression (synthetic)
        ("gated untrust action rate == 0", sup["act_rate_untrust_gated"][0] == 0.0),
        ("ungated untrust action rate > 0.1 (gate has effect)",
         sup["act_rate_untrust_ungated"][0] > 0.1),
        ("trustworthy voxels stay free to act (>0.1)", sup["act_rate_trust_gated"][0] > 0.1),
        # convergence (synthetic)
        ("trusted tumour f-drop CI excludes 0", con["drop_trusted_tumour"][1] > 0),
        ("untrusted (held) tumour f-drop == 0 on synthetic twin",
         con["drop_untrusted_tumour_held"][0] == 0.0),
        ("treatment converges (n_treat final == 0)", con["n_treat_final"] == 0),
        ("TREAT actions wind down (first > last)", tr["n_treat"][0] > tr["n_treat"][-1]),
        # Ferry grounded
        ("grounded loop closes: trusted f-drop CI excludes 0", gf["drop_trusted"][1] > 0),
        ("grounded gate suppresses (gated untrust rate == 0)", gf["sup_gated"][0] == 0.0),
        ("grounded treatment converges (n_treat final == 0)", gf["n_treat"][-1] == 0),
        ("grounded anatomy larger than synthetic", gf["n_tumor"] > sy["n_tumor"]),
        # F1 -- the honest negative (centerpiece)
        ("F1: held drop under REAL dose CI excludes 0 (>0)", f1["value_ci"][1] > 0),
        ("F1: synthetic-baseline held drop == 0", f1["synthetic_baseline_ci"][0] == 0.0),
        ("F1: real held-drop strictly exceeds synthetic baseline",
         f1["value_ci"][0] > f1["synthetic_baseline_ci"][0]),
        # byte-identity anchor
        ("loop.py blob == shipped (byte-unchanged across substrates)",
         fer.get("loop_py_blob_sha1") == LOOP_PY_BLOB_SHA1),
    ]
    for name, cond in asserts:
        print(f"  assert {name}: {'PASS' if cond else 'FAIL'}")
        ok = ok and bool(cond)
    return ok


def main():
    nums, cp2, cp4, fer = build_numbers()
    write_numbers(nums)
    print(f"wrote {NUMBERS} ({len(nums)} macros)")
    ok = verify(nums, cp2, cp4, fer)
    print("consistency gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
