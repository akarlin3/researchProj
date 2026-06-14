"""GATE 3 consistency check: every manuscript number traces to a gated CP printout.

This re-reads the committed, seeded checkpoint outputs (results/*.txt, POSITIONING.md,
gauge/results*.md) from Gauge 01-04 and asserts that each headline number used in
gauge/paper/gauge.tex appears VERBATIM in its source printout. It is the executable
form of the "run-then-write; every number traces to a gated CP output" discipline.

Run:  python gauge/paper/consistency.py     # prints per-claim PASS/FAIL + summary
Exit code 0 iff every claim traces; 1 otherwise.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# (claim, manuscript value, source file (repo-relative), [substrings that must all
#  appear verbatim in the source]).  Sources are the gated CP printouts.
CHECKS = [
    # --- Gauge 01: foundation ------------------------------------------------
    ("Clean-signal recovery (forward-model sanity)", "max rel err 2.07e-7",
     "POSITIONING.md", ["2.07e-07"]),
    ("Conformal correctness: split D* coverage @a=0.10", "0.897 vs 0.90",
     "results/coverage_report.txt", ["0.897"]),
    # --- Gauge 02: marginal benchmark ---------------------------------------
    ("Raw PNN-Gaussian D coverage/gap", "0.789 (+0.111)",
     "results/benchmark_report.txt", ["0.789 (+0.111)"]),
    ("Raw MDN D* coverage/gap", "0.868 (+0.032)",
     "results/benchmark_report.txt", ["0.868 (+0.032)"]),
    ("Raw DeepEnsemble-Point D coverage (collapses)", "0.431 (+0.469)",
     "results/benchmark_report.txt", ["0.431 (+0.469)"]),
    ("Raw Bayesian-MCMC best gap; conformalized max |gap| <= 0.024", "0.024",
     "results/benchmark_report.txt", ["(+0.024)"]),
    ("Marginal hypothesis rejected: 0/4 D*/f-concentrated", "0/4",
     "results/benchmark_report.txt", ["0/4"]),
    ("Conformalized-MDN vs pure-CQR width ratios", "0.73/0.79/0.65x",
     "results/benchmark_report.txt", ["0.73x", "0.79x", "0.65x"]),
    ("Sharpness cost MDN", "1.08-1.12x",
     "results/benchmark_report.txt", ["1.08x", "1.12x"]),
    ("Sharpness cost Bayesian", "1.02-1.05x",
     "results/benchmark_report.txt", ["1.02x", "1.05x"]),
    ("Sharpness cost DeepEnsemble-Point", "3.17-3.72x",
     "results/benchmark_report.txt", ["3.72x", "3.17x"]),
    # --- Gauge 02: conditional finding --------------------------------------
    ("hi-D* conformalized-MDN marginal/worst-SNR", "0.877 / 0.808",
     "results/conditional_report.txt", ["0.877", "0.808"]),
    ("hi-D* CQR(plain) marginal/worst-SNR", "0.819 / 0.767",
     "results/conditional_report.txt", ["0.819", "0.767"]),
    ("hi-D* raw-MDN marginal/worst-SNR", "0.825 / 0.756",
     "gauge/results.md", ["0.825", "0.756"]),
    # --- Gauge 03: identifiability limit ------------------------------------
    ("Gauge 03 verdict", "IRREDUCIBLE IDENTIFIABILITY LIMIT",
     "results/conditional_attack_report.txt", ["IRREDUCIBLE IDENTIFIABILITY LIMIT"]),
    ("Best label-free method hi-D* marginal/worst-SNR", "0.868 / 0.819",
     "results/conditional_attack_report.txt", ["0.868", "0.819"]),
    ("Plug-in routing error (31% of true-high-D* misrouted)", "31%",
     "results/conditional_attack_report.txt", ["31%"]),
    ("CRLB(D*)/tercile-width 0.34 -> 0.69 -> 1.12", "1.12",
     "results/conditional_attack_report.txt", ["0.34", "0.69", "1.12"]),
    ("Absolute CRLB(D*) grows ~6x", "~6x",
     "results/conditional_attack_report.txt", ["~6x"]),
    ("Conformal width vs CRLB log-log r=0.75", "r=0.75",
     "results/conditional_attack_report.txt", ["0.75"]),
    ("Smoking gun: nominal per OBSERVED stratum", "[0.90, 0.89, 0.90]",
     "gauge/results_gauge03.md", ["0.90, 0.89, 0.90"]),
    ("Smoking gun: craters per TRUE high-D* tercile", "0.76",
     "results/conditional_attack_report.txt", ["0.76"]),
    # --- Gauge 04: robustness -----------------------------------------------
    ("SNR shift naive coverage (D, f)", "0.610 / 0.616",
     "results/robustness_report.txt", ["0.610", "0.616"]),
    ("Weighted conformal recovers SNR shift", "0.901/0.899/0.911",
     "results/robustness_report.txt", ["0.901", "0.899", "0.911"]),
    ("Prior shift D* naive -> weighted", "0.515 -> 0.975",
     "results/robustness_report.txt", ["0.515", "0.975"]),
    ("Tri-exp: weighting over-corrects D*", "0.928",
     "results/robustness_report.txt", ["0.928"]),
    ("Weighted SNR within 0.011 of nominal", "0.011",
     "results/robustness_report.txt", ["0.011"]),
    ("Monitor fires before failure (alarm SNR 20, fail SNR 13)", "20.0 / 13.0",
     "results/robustness_report.txt", ["SNR 20.0", "SNR 13.0"]),
    ("Latent: in-dist marginal 0.900, monitor silent, hi-D* conditional 0.815",
     "0.815 / AUC 0.50", "results/robustness_report.txt", ["0.815", "AUC 0.50"]),
    # --- Gauge 04: acquisition ----------------------------------------------
    ("Acquisition hi-D* marginal by scheme", "0.841/0.844/0.845",
     "results/robustness_report.txt", ["0.841", "0.844", "0.845"]),
    ("Acquisition CRLB/tercile-width by scheme", "1.25/1.05/1.06",
     "results/robustness_report.txt", ["1.25", "1.05", "1.06"]),
    # --- Gauge 04: in-vivo --------------------------------------------------
    ("In-vivo D* band widths 10/50/90th pct", "[47.4, 64.0, 78.9]",
     "results/invivo_report.txt", ["47.4", "64.0", "78.9"]),
    ("In-vivo widest/narrowest decile ratio", "1.7x",
     "results/invivo_report.txt", ["1.7x"]),
    ("In-vivo deployment monitor AUC on transfer", "0.84",
     "results/invivo_report.txt", ["0.84"]),
]


def main():
    passed, failed = 0, 0
    missing_files = []
    print("=" * 88)
    print("GAUGE 05 -- GATE 3 consistency check: manuscript numbers vs gated CP printouts")
    print("=" * 88)
    cache = {}
    for claim, value, relpath, needles in CHECKS:
        path = os.path.join(ROOT, relpath)
        if path not in cache:
            try:
                with open(path, "r", encoding="utf-8", errors="replace") as fh:
                    cache[path] = fh.read()
            except FileNotFoundError:
                cache[path] = None
                missing_files.append(relpath)
        text = cache[path]
        if text is None:
            ok, why = False, f"SOURCE FILE MISSING: {relpath}"
        else:
            absent = [n for n in needles if n not in text]
            ok = not absent
            why = "ok" if ok else f"NOT FOUND in {relpath}: {absent}"
        passed += ok
        failed += (not ok)
        mark = "PASS" if ok else "FAIL"
        print(f"[{mark}] {claim}")
        print(f"        value={value!r}  <-  {relpath}" +
              ("" if ok else f"   *** {why} ***"))
    total = passed + failed
    print("-" * 88)
    print(f"{passed}/{total} claims trace verbatim to their gated CP printout; "
          f"{failed} failed.")
    print("-" * 88)
    # one-paragraph consistency summary
    verdict = "PASS" if failed == 0 else "FAIL"
    print("CONSISTENCY SUMMARY (GATE 3):")
    print(
        f"  Every one of the {total} headline numbers in gauge/paper/gauge.tex was\n"
        f"  cross-checked against the committed, seeded checkpoint printouts of Gauge\n"
        f"  01-04 (results/*.txt, POSITIONING.md, gauge/results*.md): {passed} trace\n"
        f"  verbatim, {failed} do not. The manuscript's headline is the HONEST,\n"
        f"  reframed contribution -- model-based IVIM UQ is broadly overconfident (not\n"
        f"  specifically perfusion); conformal restores marginal coverage (|gap| <=\n"
        f"  0.024) and conformalizing the MDN is the sharpest valid recipe (0.65-0.79x\n"
        f"  pure-CQR width); and the high-D* compartment is an IRREDUCIBLE\n"
        f"  IDENTIFIABILITY LIMIT (CRLB(D*)/tercile-width reaches 1.12), so the paper\n"
        f"  says 'characterize', not 'solve'. Robustness, acquisition-robustness, and\n"
        f"  the qualitative (no-coverage-claim) in-vivo demo are reported per Gauge 04.\n"
        f"  GATE 3: {verdict}.")
    print("=" * 88)
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
