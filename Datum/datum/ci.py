r"""Bootstrap confidence intervals for the load-bearing reference numbers.

No voxel-level bootstrap-CI helper exists in Caliper or Gauge, so Datum provides
one here. It is a standard nonparametric bootstrap over the test voxels: resample
the test set with replacement ``n_boot`` times (seeded), recompute the ruler
metrics on each resample, and report percentile intervals. CIs are attached to the
load-bearing numbers -- marginal coverage-gap and ECE per parameter, and the
per-D\*-tercile coverage (the identifiability-wall regime).

These CIs are themselves ruler-derived, hence PROVISIONAL until the ruler locks.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from datum import _paths

_paths.ensure_deps()

from caliper.metrics import central_interval, ece_quantile, empirical_coverage  # noqa: E402


@dataclass
class Interval:
    """A point estimate with a percentile bootstrap CI."""
    point: float
    lo: float
    hi: float
    level: float = 0.95

    def as_tuple(self) -> tuple[float, float, float]:
        return (self.point, self.lo, self.hi)


def _pct_ci(point: float, samples: list[float], level: float) -> Interval:
    if not samples:
        return Interval(point=point, lo=float("nan"), hi=float("nan"), level=level)
    a = (1.0 - level) / 2.0
    lo, hi = np.percentile(np.asarray(samples, dtype=float), [100 * a, 100 * (1 - a)])
    return Interval(point=float(point), lo=float(lo), hi=float(hi), level=level)


def bootstrap_reference(y_true, q_corr, q_levels, alpha, strata, *,
                        n_boot: int = 1000, seed: int = 0, level: float = 0.95):
    """Bootstrap CIs for one corrected-quantile cell.

    Parameters
    ----------
    y_true : (n, P) ground truth (Caliper convention).
    q_corr : (n, P, L) corrected predicted quantiles.
    q_levels : (L,) ascending quantile levels in (0, 1).
    alpha : central-interval miss rate (so nominal = 1 - alpha).
    strata : (n,) integer group labels (D* terciles) for conditional coverage.
    n_boot, seed, level : bootstrap resamples, RNG seed, CI level.

    Returns a dict:
      {"marginal": {p: {"coverage_gap": Interval, "ece": Interval, "coverage": Interval}},
       "by_stratum": {(p, g): Interval (coverage)}}
    where p is the parameter index and g a stratum label.
    """
    y_true = np.asarray(y_true, dtype=float)
    q_corr = np.asarray(q_corr, dtype=float)
    q_levels = np.asarray(q_levels, dtype=float)
    strata = np.asarray(strata)
    n, P, _ = q_corr.shape
    nominal = 1.0 - alpha
    groups = sorted(int(g) for g in np.unique(strata))
    rng = np.random.default_rng(seed)

    # Point estimates (full sample).
    point = {p: {} for p in range(P)}
    point_strat = {}
    lo_hi = {}
    for p in range(P):
        lo, hi = central_interval(q_corr[:, p, :], q_levels, alpha)
        lo_hi[p] = (lo, hi)
        cov = empirical_coverage(y_true[:, p], lo, hi)
        point[p]["coverage"] = cov
        point[p]["coverage_gap"] = cov - nominal
        point[p]["ece"] = ece_quantile(y_true[:, p], q_corr[:, p, :], q_levels)
        for g in groups:
            m = strata == g
            point_strat[(p, g)] = (empirical_coverage(y_true[m, p], lo[m], hi[m])
                                   if m.any() else float("nan"))

    # Bootstrap accumulation.
    samp = {p: {"coverage_gap": [], "ece": [], "coverage": []} for p in range(P)}
    samp_strat = {(p, g): [] for p in range(P) for g in groups}
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        ys, qs, gs = y_true[idx], q_corr[idx], strata[idx]
        for p in range(P):
            lo, hi = central_interval(qs[:, p, :], q_levels, alpha)
            cov = empirical_coverage(ys[:, p], lo, hi)
            samp[p]["coverage"].append(cov)
            samp[p]["coverage_gap"].append(cov - nominal)
            samp[p]["ece"].append(ece_quantile(ys[:, p], qs[:, p, :], q_levels))
            for g in groups:
                m = gs == g
                if m.any():
                    samp_strat[(p, g)].append(empirical_coverage(ys[m, p], lo[m], hi[m]))

    marginal = {
        p: {
            "coverage": _pct_ci(point[p]["coverage"], samp[p]["coverage"], level),
            "coverage_gap": _pct_ci(point[p]["coverage_gap"], samp[p]["coverage_gap"], level),
            "ece": _pct_ci(point[p]["ece"], samp[p]["ece"], level),
        }
        for p in range(P)
    }
    by_stratum = {
        (p, g): _pct_ci(point_strat[(p, g)], samp_strat[(p, g)], level)
        for p in range(P) for g in groups
    }
    return {"marginal": marginal, "by_stratum": by_stratum, "n_boot": n_boot}
