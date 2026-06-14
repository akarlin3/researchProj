"""CP1 -- conformal coverage validation + sharpness characterization (GATE 1).

Builds the synthetic IVIM cohort, then for each parameter (D, D*, f) and a sweep
of alpha:

  * split conformal  wrapped around the NLLS point estimator
  * CQR              wrapped around gradient-boosted quantile regression

and reports realized vs nominal coverage on the held-out test split. GATE 1: if
realized marginal coverage does not track nominal within tolerance, the conformal
implementation is wrong -- STOP. Then characterizes sharpness (interval width)
per parameter and across SNR, flagging where the ill-posed D*/f compartment blows
up (the regime Gauge 02's benchmark will target).

Run:  python scripts/run_coverage.py
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from gauge.cohort import generate_cohort, DEFAULT_SNR_GRID, DEFAULT_SEED
from gauge.estimators import fit_nlls_batch, IVIMQuantileRegressor
from gauge.conformal import split_conformal, cqr, empirical_coverage, interval_width

# --- experiment configuration ------------------------------------------------
N_TRAIN, N_CAL, N_TEST = 4000, 2000, 3000
ALPHAS = (0.05, 0.10, 0.20, 0.30)
PARAMS = ("D", "D*", "f")
TOL_COV = 0.03          # GATE 1 tolerance on |realized - nominal| marginal coverage
DISPLAY_SCALE = (1e3, 1e3, 1.0)   # show D, D* widths in 1e-3 mm^2/s; f absolute
WIDTH_UNITS = ("1e-3 mm^2/s", "1e-3 mm^2/s", "(fraction)")
RESULTS_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")


class Tee:
    """Write to stdout and a file at once so the printout is also archived."""
    def __init__(self, path):
        self.f = open(path, "w")

    def __call__(self, *args):
        line = " ".join(str(a) for a in args)
        print(line)
        self.f.write(line + "\n")

    def close(self):
        self.f.close()


def build_predictions(cohort):
    """Return split-conformal point preds and CQR quantile preds on cal/test."""
    t0 = time.time()
    # Split-conformal base: NLLS point estimates on calibration + test.
    cal_point = fit_nlls_batch(cohort.signals["cal"], cohort.b)
    test_point = fit_nlls_batch(cohort.signals["test"], cohort.b)
    t_nlls = time.time() - t0

    # CQR base: gradient-boosted quantile regression trained on the train split.
    levels = sorted({a / 2 for a in ALPHAS} | {1 - a / 2 for a in ALPHAS})
    t1 = time.time()
    qreg = IVIMQuantileRegressor(levels, random_state=0).fit(
        cohort.signals["train"], cohort.params["train"])
    # Predict each needed quantile for each parameter on cal + test.
    cal_q = {(j, q): qreg.predict_quantile(cohort.signals["cal"], j, q)
             for j in range(3) for q in levels}
    test_q = {(j, q): qreg.predict_quantile(cohort.signals["test"], j, q)
              for j in range(3) for q in levels}
    t_cqr = time.time() - t1
    return cal_point, test_point, cal_q, test_q, t_nlls, t_cqr


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    out = Tee(os.path.join(RESULTS_DIR, "coverage_report.txt"))

    cohort = generate_cohort(N_TRAIN, N_CAL, N_TEST, seed=DEFAULT_SEED)
    cal_true = cohort.params["cal"]
    test_true = cohort.params["test"]
    test_snr = cohort.snr["test"]

    cal_point, test_point, cal_q, test_q, t_nlls, t_cqr = build_predictions(cohort)

    out("=" * 78)
    out("CP1 -- CONFORMAL COVERAGE VALIDATION (GATE 1)")
    out("=" * 78)
    out(f"cohort seed: {cohort.seed}   b-values: {cohort.b.size}   "
        f"SNR grid (b=0): {list(DEFAULT_SNR_GRID)}")
    out(f"split sizes: train={N_TRAIN}  calibration={N_CAL}  test={N_TEST}")
    out(f"alpha sweep: {list(ALPHAS)}  (nominal coverage = 1 - alpha)")
    out(f"base estimators: split-conformal<-NLLS, CQR<-grad-boosted quantile reg")
    out(f"timing: NLLS cal+test = {t_nlls:.1f}s,  CQR train+predict = {t_cqr:.1f}s")
    out("")

    # --- coverage tables, per method ----------------------------------------
    # store widths at a representative alpha for the sharpness section
    widths_by_method = {}
    gate_failures = []

    for method in ("split-conformal", "CQR"):
        out("-" * 78)
        out(f"[{method}]  realized vs nominal coverage on TEST  (n={N_TEST})")
        out("-" * 78)
        out(f"{'param':>5} | {'alpha':>5} | {'nominal':>7} | {'realized':>8} | "
            f"{'|diff|':>6} | {'med width':>10} | {'units':>11} | {'gate':>4}")
        out("-" * 78)
        widths_here = {}
        for j, pname in enumerate(PARAMS):
            for a in ALPHAS:
                if method == "split-conformal":
                    lo, hi, _ = split_conformal(
                        cal_point[:, j], cal_true[:, j], test_point[:, j], a)
                else:
                    qlo, qhi = a / 2, 1 - a / 2
                    lo, hi, _ = cqr(
                        cal_q[(j, qlo)], cal_q[(j, qhi)], cal_true[:, j],
                        test_q[(j, qlo)], test_q[(j, qhi)], a)
                cov = empirical_coverage(lo, hi, test_true[:, j])
                w = interval_width(lo, hi)
                medw = np.median(w) * DISPLAY_SCALE[j]
                diff = abs(cov - (1 - a))
                ok = diff <= TOL_COV
                if not ok:
                    gate_failures.append((method, pname, a, cov, 1 - a))
                widths_here[(j, a)] = w
                out(f"{pname:>5} | {a:>5.2f} | {1-a:>7.2f} | {cov:>8.3f} | "
                    f"{diff:>6.3f} | {medw:>10.3f} | {WIDTH_UNITS[j]:>11} | "
                    f"{'OK' if ok else 'FAIL':>4}")
            out("-" * 78)
        widths_by_method[method] = widths_here
        out("")

    # --- GATE 1 verdict ------------------------------------------------------
    out("=" * 78)
    if gate_failures:
        out(f"GATE 1: FAIL -- {len(gate_failures)} (method,param,alpha) cells "
            f"missed nominal by > {TOL_COV}:")
        for m, p, a, cov, nom in gate_failures:
            out(f"   {m} {p} alpha={a}: realized={cov:.3f} nominal={nom:.3f}")
        out("Do NOT proceed to sharpness interpretation; the implementation is "
            "suspect (exchangeability violation, calibration leak, or base bug).")
    else:
        out(f"GATE 1: PASS -- all {len(PARAMS)*len(ALPHAS)*2} cells within "
            f"{TOL_COV} of nominal across both methods, all parameters, all alpha.")
    out("=" * 78)
    out("")

    # --- sharpness vs SNR (only meaningful once GATE 1 passes) ---------------
    rep_alpha = 0.10
    out("-" * 78)
    out(f"SHARPNESS vs SNR  (median interval width, alpha={rep_alpha}, "
        f"nominal coverage {1-rep_alpha:.2f})")
    out("with conditional (per-SNR) realized coverage -- note: conformal "
        "guarantees MARGINAL,")
    out("not conditional, coverage; per-SNR deviations are expected and motivate "
        "Mondrian CP.")
    out("-" * 78)
    snr_levels = sorted(set(int(s) for s in DEFAULT_SNR_GRID))
    for method in ("split-conformal", "CQR"):
        out(f"[{method}]")
        header = f"{'param':>5} | " + " | ".join(
            f"SNR{ s:>3}".rjust(14) for s in snr_levels)
        out(header)
        out("  (each cell: median_width [cond.coverage])")
        for j, pname in enumerate(PARAMS):
            w = widths_by_method[method][(j, rep_alpha)]
            cells = []
            for s in snr_levels:
                m = np.isclose(test_snr, s)
                if m.sum() == 0:
                    cells.append("        n/a   ")
                    continue
                if method == "split-conformal":
                    lo, hi, _ = split_conformal(
                        cal_point[:, j], cal_true[:, j], test_point[:, j], rep_alpha)
                else:
                    qlo, qhi = rep_alpha / 2, 1 - rep_alpha / 2
                    lo, hi, _ = cqr(
                        cal_q[(j, qlo)], cal_q[(j, qhi)], cal_true[:, j],
                        test_q[(j, qlo)], test_q[(j, qhi)], rep_alpha)
                covs = empirical_coverage(lo[m], hi[m], test_true[m, j])
                medw = np.median(w[m]) * DISPLAY_SCALE[j]
                cells.append(f"{medw:>8.3f}[{covs:.2f}]")
            out(f"{pname:>5} | " + " | ".join(c.rjust(14) for c in cells))
        out(f"  width units: D,D* in 1e-3 mm^2/s ; f is fraction")
        out("")

    # --- compact sharpness summary: D* blow-up vs D --------------------------
    out("-" * 78)
    out("SHARPNESS SUMMARY (alpha=0.10): ratio of low-SNR to high-SNR median "
        "width,")
    out("and ratio of D*/f widths to D width -- where the ill-posed compartment "
        "blows up.")
    out("-" * 78)
    lo_snr, hi_snr = min(snr_levels), max(snr_levels)
    for method in ("split-conformal", "CQR"):
        out(f"[{method}]")
        # recompute medians per param at lo/hi SNR
        wref = {}
        for j, pname in enumerate(PARAMS):
            w = widths_by_method[method][(j, rep_alpha)]
            m_lo = np.isclose(test_snr, lo_snr)
            m_hi = np.isclose(test_snr, hi_snr)
            med_lo = np.median(w[m_lo]) * DISPLAY_SCALE[j]
            med_hi = np.median(w[m_hi]) * DISPLAY_SCALE[j]
            wref[pname] = (med_lo, med_hi)
            out(f"  {pname:>3}: width(SNR{lo_snr})={med_lo:8.3f}  "
                f"width(SNR{hi_snr})={med_hi:8.3f}  "
                f"low/high ratio={med_lo/med_hi:6.1f}x  ({WIDTH_UNITS[j]})")
        out("")

    out("=" * 78)
    out("NOTE: D and D*/f live on different scales; compare widths within a "
        "parameter across SNR,")
    out("not across parameters. The relative (width / median truth) blow-up of "
        "D* at low SNR is the")
    out("ill-posed-compartment signal Gauge 02's conformal-vs-model-based "
        "benchmark will target.")
    out("=" * 78)

    out.close()
    _make_figure(cohort, widths_by_method, test_snr, snr_levels, rep_alpha)
    return 1 if gate_failures else 0


def _make_figure(cohort, widths_by_method, test_snr, snr_levels, rep_alpha):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:  # pragma: no cover - plotting is a nicety
        print(f"(figure skipped: {e})")
        return
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for j, (pname, ax) in enumerate(zip(PARAMS, axes)):
        for method, marker in (("split-conformal", "o"), ("CQR", "s")):
            w = widths_by_method[method][(j, rep_alpha)]
            med = [np.median(w[np.isclose(test_snr, s)]) * DISPLAY_SCALE[j]
                   for s in snr_levels]
            ax.plot(snr_levels, med, marker=marker, label=method)
        ax.set_title(f"{pname} interval width (alpha={rep_alpha})")
        ax.set_xlabel("SNR (b=0)")
        ax.set_ylabel(f"median width [{WIDTH_UNITS[j]}]")
        ax.set_xscale("log")
        ax.grid(True, alpha=0.3)
        if j == 0:
            ax.legend()
    fig.tight_layout()
    path = os.path.join(RESULTS_DIR, "sharpness_vs_snr.png")
    fig.savefig(path, dpi=130)
    print(f"(figure written: {path})")


if __name__ == "__main__":
    raise SystemExit(main())
