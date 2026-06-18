"""Caliper benchmark report -- figures regenerated *solely* from the CSV.

    results/benchmark.csv  ->  examples/figures/*.png

Run (after `python -m caliper.benchmark` has written the CSV):

    python examples/benchmark_report.py
    python examples/benchmark_report.py --csv results/benchmark.csv --snr 80

This script reads **only** ``results/benchmark.csv`` -- it never re-runs an
estimator or touches the toolkit's compute path, so the figures are a faithful
rendering of the table the harness produced. Three figures per estimator:

1. calibration-vs-SNR -- marginal coverage vs SNR (where calibration degrades
   as SNR drops, and how each calibration method restores it);
2. coverage-by-D*-tercile bars -- raw vs CQR vs Mondrian at a fixed SNR (the
   high-D* conditional gap marginal CQR leaves behind);
3. Mondrian width cost -- the sharpness price of per-tercile validity at high D*.

INTERNAL TOOL DEMO, not a benchmark release. Every number plotted traces to a
row of the CSV, which itself traces to a fixed-seed run of caliper.benchmark.
"""
from __future__ import annotations

import argparse
import csv
import math
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")  # headless: write PNGs, never open a window

import matplotlib.pyplot as plt  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
DEFAULT_CSV = os.path.join(_ROOT, "results", "benchmark.csv")
DEFAULT_FIGDIR = os.path.join(_HERE, "figures")

# tercile order + display labels (must match caliper.benchmark.STRATUM_NAMES)
TERCILES = ("dstar_lo", "dstar_mid", "dstar_hi")
TERCILE_LABELS = ("low D*", "mid D*", "high D*")
_NUMERIC = {"snr", "n", "coverage", "coverage_gap", "mean_width", "ece",
            "mean_pinball", "mean_interval_score", "alpha", "nominal"}


# --------------------------------------------------------------------------- #
# CSV -> typed rows (the only input)
# --------------------------------------------------------------------------- #
def load_rows(path: str) -> list[dict]:
    if not os.path.exists(path):
        raise SystemExit(
            f"CSV not found: {path}\nRun `python -m caliper.benchmark` first.")
    rows = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            out = dict(r)
            for k in list(out):
                if k in _NUMERIC:
                    out[k] = float(out[k]) if out[k] not in ("", "nan") else math.nan
            out["seed"] = int(float(r["seed"]))
            rows.append(out)
    if not rows:
        raise SystemExit(f"CSV is empty: {path}")
    return rows


def _agg(rows, value, **eq):
    """Seed-averaged (mean, std) of ``value`` over rows matching ``eq``."""
    vals = [r[value] for r in rows
            if all(r[k] == v for k, v in eq.items())
            and not (isinstance(r[value], float) and math.isnan(r[value]))]
    if not vals:
        return math.nan, math.nan
    return float(np.mean(vals)), float(np.std(vals))


def _snrs(rows, estimator):
    return sorted({r["snr"] for r in rows if r["estimator"] == estimator})


def _calibs_present(rows, estimator):
    order = ["raw", "split", "CQR", "Mondrian"]
    have = {r["calibration"] for r in rows if r["estimator"] == estimator}
    return [c for c in order if c in have]


# --------------------------------------------------------------------------- #
# Figure 1: calibration vs SNR
# --------------------------------------------------------------------------- #
def fig_calibration_vs_snr(rows, estimator, nominal, outpath):
    snrs = _snrs(rows, estimator)
    calibs = _calibs_present(rows, estimator)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.4))

    # left: D* marginal coverage vs SNR, one line per calibration method
    for calib in calibs:
        ys, es = zip(*[_agg(rows, "coverage", estimator=estimator, calibration=calib,
                            param="Dstar", stratum="all", snr=s) for s in snrs])
        axL.errorbar(snrs, ys, yerr=es, marker="o", capsize=3, label=calib)
    axL.axhline(nominal, ls="--", color="k", lw=1, label=f"nominal {nominal:.2f}")
    axL.set_xscale("log")
    axL.set_xlabel("SNR (at b=0)")
    axL.set_ylabel("marginal coverage (D*, 90% interval)")
    axL.set_title("Calibration vs SNR: raw degrades, conformal restores")
    axL.set_ylim(0, 1.02)
    axL.legend(fontsize=8, loc="center right")
    axL.grid(alpha=0.3)

    # right: raw coverage vs SNR for all three parameters (where it degrades)
    for param in ("D", "f", "Dstar"):
        ys, es = zip(*[_agg(rows, "coverage", estimator=estimator, calibration="raw",
                            param=param, stratum="all", snr=s) for s in snrs])
        axR.errorbar(snrs, ys, yerr=es, marker="s", capsize=3, label=param)
    axR.axhline(nominal, ls="--", color="k", lw=1, label=f"nominal {nominal:.2f}")
    axR.set_xscale("log")
    axR.set_xlabel("SNR (at b=0)")
    axR.set_ylabel("raw marginal coverage")
    axR.set_title("Raw over-confidence by parameter")
    axR.set_ylim(0, 1.02)
    axR.legend(fontsize=8, loc="lower right")
    axR.grid(alpha=0.3)

    fig.suptitle(f"Caliper benchmark ({estimator}) -- INTERNAL TOOL DEMO, not a release",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outpath, dpi=120)
    plt.close(fig)
    return outpath


# --------------------------------------------------------------------------- #
# Figure 2: coverage by D* tercile (raw vs CQR vs Mondrian)
# --------------------------------------------------------------------------- #
def fig_coverage_by_tercile(rows, estimator, nominal, snr, outpath):
    methods = [m for m in ("raw", "CQR", "Mondrian")
               if m in _calibs_present(rows, estimator)]
    x = np.arange(len(TERCILES))
    w = 0.8 / len(methods)
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    for i, meth in enumerate(methods):
        ys, es = zip(*[_agg(rows, "coverage", estimator=estimator, calibration=meth,
                            param="Dstar", stratum=t, snr=snr) for t in TERCILES])
        ax.bar(x + i * w, ys, w, yerr=es, capsize=3, label=meth)
    ax.axhline(nominal, ls="--", color="k", lw=1, label=f"nominal {nominal:.2f}")
    ax.set_xticks(x + w * (len(methods) - 1) / 2)
    ax.set_xticklabels(TERCILE_LABELS)
    ax.set_ylabel("conditional coverage (D*)")
    ax.set_ylim(0, 1.05)
    ax.set_title(f"D* coverage by tercile @ SNR={snr:g} ({estimator})\n"
                 "marginal CQR leaves a high-D* gap; Mondrian equalizes")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(outpath, dpi=120)
    plt.close(fig)
    return outpath


# --------------------------------------------------------------------------- #
# Figure 3: Mondrian width cost at high D*
# --------------------------------------------------------------------------- #
def fig_mondrian_width_cost(rows, estimator, snr, outpath):
    snrs = _snrs(rows, estimator)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.4))

    # left: mean D* interval width by tercile, CQR vs Mondrian, at one SNR
    methods = [m for m in ("CQR", "Mondrian") if m in _calibs_present(rows, estimator)]
    x = np.arange(len(TERCILES))
    w = 0.8 / max(len(methods), 1)
    for i, meth in enumerate(methods):
        ys = [_agg(rows, "mean_width", estimator=estimator, calibration=meth,
                   param="Dstar", stratum=t, snr=snr)[0] for t in TERCILES]
        axL.bar(x + i * w, ys, w, label=meth)
    axL.set_xticks(x + w * (len(methods) - 1) / 2)
    axL.set_xticklabels(TERCILE_LABELS)
    axL.set_yscale("log")
    axL.set_ylabel("mean D* interval width (log)")
    axL.set_title(f"Width by tercile @ SNR={snr:g}")
    axL.legend(fontsize=8)
    axL.grid(alpha=0.3, axis="y", which="both")

    # right: high/low D* width ratio vs SNR (CQR ~1; Mondrian inflates)
    for meth in methods:
        ratios = []
        for s in snrs:
            hi = _agg(rows, "mean_width", estimator=estimator, calibration=meth,
                      param="Dstar", stratum="dstar_hi", snr=s)[0]
            lo = _agg(rows, "mean_width", estimator=estimator, calibration=meth,
                      param="Dstar", stratum="dstar_lo", snr=s)[0]
            ratios.append(hi / lo if lo and not math.isnan(lo) else math.nan)
        axR.plot(snrs, ratios, marker="o", label=meth)
    axR.axhline(1.0, ls="--", color="k", lw=1, label="equal width")
    axR.set_xscale("log")
    axR.set_xlabel("SNR (at b=0)")
    axR.set_ylabel("high-D* / low-D* width ratio")
    axR.set_title("Sharpness cost of conditional validity")
    axR.legend(fontsize=8)
    axR.grid(alpha=0.3)

    fig.suptitle(f"Mondrian width cost ({estimator}) -- the price of per-tercile "
                 "coverage at the identifiability wall", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(outpath, dpi=120)
    plt.close(fig)
    return outpath


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def main(argv=None):
    ap = argparse.ArgumentParser(description="Regenerate Caliper benchmark figures from CSV.")
    ap.add_argument("--csv", default=DEFAULT_CSV)
    ap.add_argument("--figdir", default=DEFAULT_FIGDIR)
    ap.add_argument("--snr", type=float, default=None,
                    help="SNR for the tercile/width figures (default: highest in CSV)")
    args = ap.parse_args(argv)

    rows = load_rows(args.csv)
    os.makedirs(args.figdir, exist_ok=True)
    estimators = sorted({r["estimator"] for r in rows})
    nominal = next(r["nominal"] for r in rows)

    written = []
    for est in estimators:
        snrs = _snrs(rows, est)
        snr = args.snr if args.snr is not None else snrs[-1]  # highest SNR: gap clearest
        suffix = f"_{est}" if len(estimators) > 1 else ""
        written.append(fig_calibration_vs_snr(
            rows, est, nominal, os.path.join(args.figdir, f"calibration_vs_snr{suffix}.png")))
        written.append(fig_coverage_by_tercile(
            rows, est, nominal, snr,
            os.path.join(args.figdir, f"coverage_by_dstar_tercile{suffix}.png")))
        written.append(fig_mondrian_width_cost(
            rows, est, snr, os.path.join(args.figdir, f"mondrian_width_cost{suffix}.png")))

    print("Caliper benchmark report -- figures regenerated from", args.csv)
    for p in written:
        print("  wrote", os.path.relpath(p, _ROOT))
    print("\n(reads only the CSV; every value traces to a fixed-seed benchmark run.)")


if __name__ == "__main__":
    main()
