"""caliper.benchmark -- a reproducible evaluation harness for the toolkit.

INTERNAL TOOL INFRASTRUCTURE. This is a reproducible eval/demo that exercises
the Caliper toolkit end-to-end; it is **not** a citable benchmark, leaderboard,
or dataset release. No external/clinical data -- synthetic only, fixed seeds.

What it does
------------
Sweeps a grid

    {estimator} x {calibration: raw, split, CQR, Mondrian} x {SNR} x {seed}

and, for each cell, scores the (corrected) predictive quantiles with the
calibration ruler (:mod:`caliper.metrics`) and adds D*-tercile *conditional*
coverage and mean interval width via :mod:`caliper.conformal`. The result is a
tidy long-form table written to ``results/benchmark.csv`` plus a printed
summary.

Design contract
---------------
Every calibration method is expressed as a map from the estimator's *raw*
quantile array ``(n, P, L)`` (plus a calibration split) to a *corrected*
quantile array of the same shape, so the ruler scores all four identically:

* ``raw``      -- identity (the estimator's own, deliberately over-confident,
                  quantiles).
* ``split``    -- split-conformal *residual* intervals (point +/- Q) assembled
                  across every symmetric level-pair, reusing
                  :class:`caliper.conformal.SplitConformalResidual`. The point is
                  the estimator's median predicted quantile, so this works for
                  any estimator under the ``predict_quantiles`` contract.
* ``CQR``      -- marginal conformalized quantile regression
                  (:class:`caliper.conformal.SplitConformalQuantile`).
* ``Mondrian`` -- group-conditional CQR
                  (:class:`caliper.conformal.MondrianConformalQuantile`), grouped
                  by the *true* D* tercile. This is an **oracle** stratification
                  (it peeks at the true D* a deployment would not know); it is
                  used here only to expose the conditional-coverage / width trade
                  at the identifiability wall.

Estimators are addressed only through ``predict_quantiles(signals, q_levels)``.
The numpy-only reference estimator is always available; the torch MAF is added
automatically **iff** torch is importable (detected without importing the torch
module, so a torch-free environment never trips an import error).

Reproducibility
---------------
All cohorts are drawn with seeds derived deterministically from ``(seed, SNR
index)``; the reference estimator is closed-form and deterministic. Re-running
``run_grid`` with the same arguments yields byte-identical rows
(:func:`check_reproducible`).

Run
---
    python -m caliper.benchmark              # default grid -> results/benchmark.csv
    python -m caliper.benchmark --quick      # tiny grid (smoke test)
    python -m caliper.benchmark --check      # assert fixed seeds reproduce the table
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import math
import os
from typing import Optional, Sequence

import numpy as np

from . import conformal as C
from . import metrics as M
from .estimator_reference import ReferenceIVIMEstimator
from .forward import DEFAULT_BVALUES, PARAM_NAMES, synthetic_cohort

__all__ = [
    "DEFAULT_LEVELS",
    "DEFAULT_SNRS",
    "DEFAULT_SEEDS",
    "DEFAULT_ALPHA",
    "CALIBRATIONS",
    "CSV_COLUMNS",
    "STRATUM_NAMES",
    "torch_available",
    "default_estimators",
    "run_grid",
    "write_csv",
    "summarize",
    "check_reproducible",
]

# --------------------------------------------------------------------------- #
# Grid + schema constants
# --------------------------------------------------------------------------- #
DEFAULT_LEVELS = np.array([0.05, 0.25, 0.5, 0.75, 0.95])
DEFAULT_SNRS = (10.0, 20.0, 40.0, 80.0)
DEFAULT_SEEDS = (0, 1, 2)
DEFAULT_ALPHA = 0.10
CALIBRATIONS = ("raw", "split", "CQR", "Mondrian")

DSTAR = PARAM_NAMES.index("Dstar")
# tercile labels of the *true* D* used for conditional coverage everywhere.
STRATUM_NAMES = {0: "dstar_lo", 1: "dstar_mid", 2: "dstar_hi"}
ALL_STRATUM = "all"

# tidy long-form schema: one row per (estimator, calibration, snr, seed, param,
# stratum). Marginal ruler metrics live on the ``all`` stratum; the dstar_*
# strata carry per-tercile coverage + width only.
CSV_COLUMNS = [
    "estimator",
    "calibration",
    "snr",
    "seed",
    "param",
    "stratum",
    "n",
    "coverage",
    "coverage_gap",
    "mean_width",
    "ece",
    "mean_pinball",
    "mean_interval_score",
    "alpha",
    "nominal",
]

# Default output location, resolved relative to the package so cwd does not
# matter: <Caliper>/results/benchmark.csv .
_PKG_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_CSV_PATH = os.path.join(_PKG_ROOT, "results", "benchmark.csv")


# --------------------------------------------------------------------------- #
# Estimator registry (estimators are used only via predict_quantiles)
# --------------------------------------------------------------------------- #
def torch_available() -> bool:
    """True iff torch can be imported, checked *without* importing it."""
    return importlib.util.find_spec("torch") is not None


def default_estimators() -> tuple[str, ...]:
    """Reference always; MAF appended iff torch is importable."""
    return ("reference",) + (("maf",) if torch_available() else ())


def _build_estimator(name: str, bvalues: np.ndarray, seed: int):
    """Return ``(estimator, needs_training)`` for a registry name."""
    if name == "reference":
        return ReferenceIVIMEstimator(bvalues=np.asarray(bvalues, dtype=float)), False
    if name == "maf":
        if not torch_available():
            raise RuntimeError("estimator 'maf' requires torch (pip install -e '.[estimator]')")
        from .estimator_maf import MAFPosterior  # lazy: torch-only

        return MAFPosterior(n_bvalues=int(np.asarray(bvalues).shape[0]), seed=seed), True
    raise ValueError(f"unknown estimator {name!r}; known: reference, maf")


def _cell_seeds(seed: int, snr_index: int) -> dict[str, int]:
    """Deterministic, collision-free cohort seeds for one grid cell."""
    base = 10_000 + seed * 100 + snr_index * 10
    return {"train": base, "cal": base + 1, "test": base + 2}


# --------------------------------------------------------------------------- #
# Calibration: raw quantiles -> corrected quantiles (n, P, L)
# --------------------------------------------------------------------------- #
def _split_residual_quantiles(levels, q_cal, y_cal, q_test) -> np.ndarray:
    """Build a full quantile array from split-conformal *residual* intervals.

    The point estimate is the estimator's median predicted quantile (the level
    nearest 0.5). For each symmetric level pair ``(j_lo, j_hi)`` with miss-rate
    ``a = 2 * levels[j_lo]`` we fit an absolute-residual offset on the
    calibration split (reusing :class:`SplitConformalResidual`) and place
    ``point -/+ Q`` in the two columns. The result has symmetric, input-constant
    widths -- the conformal baseline against which CQR's adaptive widths are read.
    """
    levels = np.asarray(levels, dtype=float)
    L = levels.shape[0]
    mid = int(np.argmin(np.abs(levels - 0.5)))
    point_cal = q_cal[:, :, mid]
    point_test = q_test[:, :, mid]
    out = q_test.copy()
    out[:, :, mid] = point_test  # keep an odd-length centre at the point
    for j in range(L // 2):
        j_lo, j_hi = j, L - 1 - j
        a = 2.0 * levels[j_lo]  # nominal miss-rate of this central pair
        lo, hi = C.SplitConformalResidual(alpha=a).calibrate_apply(point_cal, y_cal, point_test)
        out[:, :, j_lo] = lo
        out[:, :, j_hi] = hi
    return np.sort(out, axis=2)  # enforce non-crossing quantiles


def _apply_calibration(
    name, levels, q_cal, y_cal, q_test, groups_cal, groups_test
) -> np.ndarray:
    """Map raw quantiles to corrected quantiles for one calibration method."""
    if name == "raw":
        return q_test.copy()
    if name == "split":
        return _split_residual_quantiles(levels, q_cal, y_cal, q_test)
    if name == "CQR":
        return C.SplitConformalQuantile(levels).calibrate_apply(q_cal, y_cal, q_test)
    if name == "Mondrian":
        return C.MondrianConformalQuantile(levels).calibrate_apply(
            q_cal, y_cal, groups_cal, q_test, groups_test
        )
    raise ValueError(f"unknown calibration {name!r}; known: {CALIBRATIONS}")


# --------------------------------------------------------------------------- #
# Scoring one cell -> rows
# --------------------------------------------------------------------------- #
def _row(estimator, calibration, snr, seed, param, stratum, alpha, *, n,
         coverage, mean_width, coverage_gap=None, ece=None, mean_pinball=None,
         mean_interval_score=None) -> dict:
    nominal = 1.0 - alpha
    if coverage_gap is None and coverage is not None:
        coverage_gap = coverage - nominal
    return {
        "estimator": estimator,
        "calibration": calibration,
        "snr": float(snr),
        "seed": int(seed),
        "param": param,
        "stratum": stratum,
        "n": int(n),
        "coverage": _f(coverage),
        "coverage_gap": _f(coverage_gap),
        "mean_width": _f(mean_width),
        "ece": _f(ece),
        "mean_pinball": _f(mean_pinball),
        "mean_interval_score": _f(mean_interval_score),
        "alpha": float(alpha),
        "nominal": float(nominal),
    }


def _f(x) -> Optional[float]:
    return None if x is None else float(x)


def _score_cell(estimator, calibration, snr, seed, levels, alpha, q_corr,
                y_test, strata) -> list[dict]:
    """Score a corrected quantile array into long-form rows (all + D* terciles)."""
    rows: list[dict] = []
    n_test = y_test.shape[0]
    scores = M.score_quantiles(y_test, q_corr, levels, alpha=alpha,
                               param_names=PARAM_NAMES)
    for p, s in enumerate(scores):
        # marginal ("all") row carries the full ruler metric set.
        rows.append(_row(
            estimator, calibration, snr, seed, s.name, ALL_STRATUM, alpha,
            n=n_test, coverage=s.coverage, coverage_gap=s.coverage_gap,
            mean_width=s.sharpness, ece=s.ece, mean_pinball=s.mean_pinball,
            mean_interval_score=s.mean_interval_score,
        ))
        # D*-tercile conditional coverage + width, via the conformal helper.
        lo, hi = M.central_interval(q_corr[:, p, :], levels, alpha)
        by_stratum = C.conditional_coverage_by_strata(y_test[:, p], lo, hi, strata)
        for g, sc in sorted(by_stratum.items()):
            rows.append(_row(
                estimator, calibration, snr, seed, s.name,
                STRATUM_NAMES.get(g, f"dstar_g{g}"), alpha,
                n=sc.n, coverage=sc.coverage, mean_width=sc.mean_width,
            ))
    return rows


# --------------------------------------------------------------------------- #
# The grid
# --------------------------------------------------------------------------- #
def run_grid(
    estimators: Optional[Sequence[str]] = None,
    calibrations: Sequence[str] = CALIBRATIONS,
    snrs: Sequence[float] = DEFAULT_SNRS,
    seeds: Sequence[int] = DEFAULT_SEEDS,
    levels=DEFAULT_LEVELS,
    alpha: float = DEFAULT_ALPHA,
    noise: str = "rician",
    n_cal: int = 4000,
    n_test: int = 9000,
    n_train: int = 8000,
    bvalues=DEFAULT_BVALUES,
    verbose: bool = True,
) -> list[dict]:
    """Run the full sweep and return tidy long-form rows.

    For each ``(estimator, snr, seed)`` the raw quantiles are predicted *once* on
    a calibration and a test cohort drawn at that SNR; all calibration methods
    then correct the *same* raw quantiles, so they are strictly comparable. The
    conditional-coverage stratification is the true-D* tercile of the test
    cohort throughout.
    """
    estimators = tuple(estimators) if estimators is not None else default_estimators()
    levels = np.asarray(levels, dtype=float)
    bvalues = np.asarray(bvalues, dtype=float)
    rows: list[dict] = []

    for est_name in estimators:
        for si, snr in enumerate(snrs):
            for seed in seeds:
                cs = _cell_seeds(seed, si)
                est, needs_train = _build_estimator(est_name, bvalues, seed)
                if needs_train:
                    train = synthetic_cohort(n=n_train, bvalues=bvalues, snr=snr,
                                             noise=noise, seed=cs["train"])
                    est.fit(train.signals, train.params)
                cal = synthetic_cohort(n=n_cal, bvalues=bvalues, snr=snr,
                                       noise=noise, seed=cs["cal"])
                test = synthetic_cohort(n=n_test, bvalues=bvalues, snr=snr,
                                        noise=noise, seed=cs["test"])
                q_cal = est.predict_quantiles(cal.signals, levels)
                q_test = est.predict_quantiles(test.signals, levels)
                groups_cal = M.tercile_groups(cal.params[:, DSTAR])
                strata = M.tercile_groups(test.params[:, DSTAR])  # oracle D* terciles
                for calib in calibrations:
                    q_corr = _apply_calibration(calib, levels, q_cal, cal.params,
                                                q_test, groups_cal, strata)
                    rows.extend(_score_cell(est_name, calib, snr, seed, levels,
                                            alpha, q_corr, test.params, strata))
                if verbose:
                    d = next(r for r in rows[::-1] if r["param"] == "Dstar"
                             and r["stratum"] == ALL_STRATUM)
                    print(f"  [{est_name}] SNR={snr:>5g} seed={seed} done "
                          f"(D* raw->Mondrian last cov={d['coverage']:.3f})")
    return rows


# --------------------------------------------------------------------------- #
# CSV I/O (numpy + stdlib only)
# --------------------------------------------------------------------------- #
def _csv_value(v) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        if math.isnan(v):
            return ""
        return f"{v:.10g}"
    return str(v)


def write_csv(rows: Sequence[dict], path: str = DEFAULT_CSV_PATH) -> str:
    """Write rows to ``path`` (creating parent dirs). Returns the path."""
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(CSV_COLUMNS)
        for r in rows:
            w.writerow([_csv_value(r[c]) for c in CSV_COLUMNS])
    return path


# --------------------------------------------------------------------------- #
# Printed summary
# --------------------------------------------------------------------------- #
def _mean(vals) -> float:
    vals = [v for v in vals if v is not None and not (isinstance(v, float) and math.isnan(v))]
    return sum(vals) / len(vals) if vals else float("nan")


def _select(rows, **eq):
    return [r for r in rows if all(r[k] == v for k, v in eq.items())]


def summarize(rows: Sequence[dict], param: str = "Dstar") -> str:
    """Render the headline patterns for one parameter (default D*) as text."""
    estimators = sorted({r["estimator"] for r in rows})
    calibs = [c for c in CALIBRATIONS if any(r["calibration"] == c for r in rows)]
    snrs = sorted({r["snr"] for r in rows})
    seeds = sorted({r["seed"] for r in rows})
    alpha = rows[0]["alpha"] if rows else DEFAULT_ALPHA
    nominal = 1.0 - alpha

    out: list[str] = []
    out.append("=== Caliper benchmark summary (INTERNAL TOOL DEMO, not a release) ===")
    out.append(f"param={param}  nominal central coverage={nominal:.3f} (alpha={alpha:.3f})")
    out.append(f"estimators={estimators}  SNRs={[f'{s:g}' for s in snrs]}  "
               f"seeds={seeds}  (seed-averaged below)")

    for est in estimators:
        out.append(f"\n--- estimator: {est} ---")
        # marginal coverage: calibration (rows) x SNR (cols)
        head = f"{'marginal cov':>14} | " + " | ".join(f"SNR{ s:>6g}" for s in snrs)
        out.append(head)
        out.append("-" * len(head))
        for calib in calibs:
            cells = []
            for s in snrs:
                sel = _select(rows, estimator=est, calibration=calib, snr=s,
                              param=param, stratum=ALL_STRATUM)
                cells.append(f"{_mean([r['coverage'] for r in sel]):>9.3f}")
            out.append(f"{calib:>14} | " + " | ".join(cells))

        # D*-tercile conditional coverage at the most-stressed (lowest) SNR
        s_lo = snrs[0]
        out.append(f"\nconditional coverage by D* tercile @ SNR={s_lo:g} "
                   f"(low/mid/high):")
        for calib in calibs:
            tri = []
            for g in (0, 1, 2):
                sel = _select(rows, estimator=est, calibration=calib, snr=s_lo,
                              param=param, stratum=STRATUM_NAMES[g])
                tri.append(f"{_mean([r['coverage'] for r in sel]):.3f}")
            out.append(f"{calib:>14}  {'/'.join(tri)}")

        # width inflation: high/low D* width ratio for CQR vs Mondrian, by SNR
        out.append("\nhigh/low D* mean-width ratio (CQR vs Mondrian), by SNR:")
        rhead = f"{'':>14} | " + " | ".join(f"SNR{s:>6g}" for s in snrs)
        out.append(rhead)
        for calib in ("CQR", "Mondrian"):
            if calib not in calibs:
                continue
            cells = []
            for s in snrs:
                hi = _mean([r["mean_width"] for r in _select(
                    rows, estimator=est, calibration=calib, snr=s, param=param,
                    stratum=STRATUM_NAMES[2])])
                lo = _mean([r["mean_width"] for r in _select(
                    rows, estimator=est, calibration=calib, snr=s, param=param,
                    stratum=STRATUM_NAMES[0])])
                ratio = hi / lo if lo and not math.isnan(lo) and lo != 0 else float("nan")
                cells.append(f"{ratio:>8.2f}x")
            out.append(f"{calib:>14} | " + " | ".join(cells))

    return "\n".join(out)


# --------------------------------------------------------------------------- #
# Reproducibility
# --------------------------------------------------------------------------- #
def _rows_equal(a: Sequence[dict], b: Sequence[dict]) -> bool:
    if len(a) != len(b):
        return False
    for ra, rb in zip(a, b):
        if set(ra) != set(rb):
            return False
        for k in ra:
            va, vb = ra[k], rb[k]
            if isinstance(va, float) and isinstance(vb, float):
                if math.isnan(va) and math.isnan(vb):
                    continue
                if va != vb:
                    return False
            elif va != vb:
                return False
    return True


def check_reproducible(**grid_kwargs) -> bool:
    """Run the grid twice and assert the rows are byte-identical."""
    grid_kwargs.setdefault("verbose", False)
    a = run_grid(**grid_kwargs)
    b = run_grid(**grid_kwargs)
    ok = _rows_equal(a, b)
    if not ok:
        raise AssertionError("benchmark grid is not reproducible across runs")
    return ok


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description="Caliper evaluation harness (tool demo).")
    ap.add_argument("--out", default=DEFAULT_CSV_PATH, help="output CSV path")
    ap.add_argument("--estimators", nargs="+", default=None,
                    help="estimator names (default: reference [+maf if torch])")
    ap.add_argument("--snrs", nargs="+", type=float, default=list(DEFAULT_SNRS))
    ap.add_argument("--seeds", nargs="+", type=int, default=list(DEFAULT_SEEDS))
    ap.add_argument("--n-cal", type=int, default=4000)
    ap.add_argument("--n-test", type=int, default=9000)
    ap.add_argument("--quick", action="store_true",
                    help="tiny grid (1 SNR pair, 1 seed, small cohorts) for smoke tests")
    ap.add_argument("--check", action="store_true",
                    help="assert fixed seeds reproduce the table, then exit")
    args = ap.parse_args(argv)

    kw = dict(estimators=args.estimators, snrs=tuple(args.snrs),
              seeds=tuple(args.seeds), n_cal=args.n_cal, n_test=args.n_test)
    if args.quick:
        kw.update(snrs=(20.0, 80.0), seeds=(0,), n_cal=800, n_test=1500)

    if args.check:
        check_reproducible(**kw)
        print("[benchmark] reproducibility check passed (two runs identical).")
        return

    print("Caliper benchmark -- INTERNAL TOOL DEMO (synthetic only, not a release)\n")
    rows = run_grid(**kw)
    path = write_csv(rows, args.out)
    print(f"\nwrote {len(rows)} rows -> {path}\n")
    print(summarize(rows))


if __name__ == "__main__":
    main()
