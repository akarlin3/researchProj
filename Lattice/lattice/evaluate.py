"""Standardised calibration-evaluation interface for the Lattice DRO.

Lattice is the *reference object* (ground-truth cohorts + generators); the
canonical *scorer* is **Caliper** (`caliper.metrics`, coverage / ECE /
sharpness). To keep the dependency strictly one-way (Caliper and the papers
consume Lattice; Lattice imports nothing back), this module deliberately does
**not** import Caliper. Instead it defines:

1. the estimator contract any UQ method must satisfy to be scored on Lattice, and
2. :func:`to_scorer_inputs`, an adapter that hands a cohort's ground truth and a
   method's predicted quantiles to *any* quantile scorer in exactly the shape
   ``caliper.metrics.score_quantiles`` expects.

A worked example wiring this to Caliper lives in ``examples/`` (Caliper is an
optional example-only dependency, never a core import). The two tiny helpers
below (:func:`interval_coverage`, :func:`mean_sharpness`) exist only for
dependency-free smoke checks; production scoring should use Caliper.
"""

from __future__ import annotations

from typing import Protocol, runtime_checkable, Sequence

import numpy as np

__all__ = [
    "IVIMQuantileEstimator",
    "DEFAULT_QUANTILE_LEVELS",
    "to_scorer_inputs",
    "interval_coverage",
    "mean_sharpness",
    "central_interval",
]

# Symmetric quantile grid for a nominal 0.90 central interval (and the median).
DEFAULT_QUANTILE_LEVELS = (0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)


@runtime_checkable
class IVIMQuantileEstimator(Protocol):
    """Contract for any UQ method scored on Lattice.

    Implementations predict, for each input signal, a quantile of each
    parameter at each requested level.
    """

    def predict_quantiles(
        self, signals: np.ndarray, q_levels: Sequence[float]
    ) -> np.ndarray:
        """Return an ``(n, n_params, n_levels)`` array of predicted quantiles.

        ``signals`` is ``(n, n_b)``; ``q_levels`` is ascending in (0, 1); the
        parameter axis follows the cohort's ``param_names`` order
        (``(D, Dstar, f)``).
        """
        ...


def to_scorer_inputs(cohort, q_pred: np.ndarray, q_levels=DEFAULT_QUANTILE_LEVELS) -> dict:
    """Package a cohort + predicted quantiles for a quantile calibration scorer.

    The returned dict matches the keyword arguments of
    ``caliper.metrics.score_quantiles(y_true, q_pred, q_levels, param_names=...)``,
    so a caller scores Lattice with Caliper in one line::

        from caliper import metrics as M
        M.score_quantiles(**lattice.evaluate.to_scorer_inputs(cohort, q_pred))
    """
    q_pred = np.asarray(q_pred, dtype=float)
    y_true = cohort.params
    n, n_params = y_true.shape
    expected = (n, n_params, len(q_levels))
    if q_pred.shape != expected:
        raise ValueError(
            f"q_pred shape {q_pred.shape} != expected {expected} "
            f"(n, n_params, n_levels)"
        )
    return {
        "y_true": y_true,
        "q_pred": q_pred,
        "q_levels": np.asarray(q_levels, dtype=float),
        "param_names": list(cohort.param_names),
    }


def central_interval(q_pred: np.ndarray, q_levels, alpha: float = 0.10):
    """Extract the central ``1 - alpha`` ``(lo, hi)`` interval from quantiles."""
    q_levels = np.asarray(q_levels, dtype=float)
    lo_level, hi_level = alpha / 2.0, 1.0 - alpha / 2.0
    i_lo = int(np.argmin(np.abs(q_levels - lo_level)))
    i_hi = int(np.argmin(np.abs(q_levels - hi_level)))
    return q_pred[..., i_lo], q_pred[..., i_hi]


def interval_coverage(y_true: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """Per-parameter empirical coverage of ``[lo, hi]`` (dependency-free smoke check)."""
    inside = (y_true >= lo) & (y_true <= hi)
    return inside.mean(axis=0)


def mean_sharpness(lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """Per-parameter mean interval width (dependency-free smoke check)."""
    return (hi - lo).mean(axis=0)
