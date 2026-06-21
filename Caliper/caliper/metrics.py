"""caliper.metrics -- the calibration ruler.

A small, dependency-light (numpy-only) toolkit for *scoring* the calibration of
predicted quantiles. It is deliberately model-agnostic: it knows nothing about
IVIM, normalizing flows, or conformal prediction. Feed it true values and
predicted quantiles and it returns coverage, calibration error, sharpness, and
conditional (group-wise) coverage.

The central data contract is the quantile array produced by an estimator's
``predict_quantiles``:

    q_pred : ndarray, shape (n, n_params, n_levels)
    q_levels : ndarray, shape (n_levels,)   ascending, in (0, 1)
    y_true : ndarray, shape (n, n_params)

Everything here is pure numpy and deterministic. Run ``python metrics.py`` for a
self-test plus a short demo.

Scope (retooled Fashion). This ruler is a *scoped secondary* instrument: it
scores intervals against **ground truth**, so it lives on synthetic / DRO data
and cannot be applied to a real scan that has no labels. Under the retooled
Fashion the load-bearing readout is *conditional* coverage by D* tercile (the
high-D* identifiability wall), not a marginal severity; the assumption-free
primary -- boundary railing on real data -- needs no ruler at all. Intervals fed
in should already follow the honest-CRLB SD convention (wide where D* is
unidentified); see ``caliper.baselines`` ``sd_convention``.

Author: Caliper project. Numpy-only by design -- do not add heavy deps here.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np

__all__ = [
    "pinball_loss",
    "interval_score",
    "empirical_coverage",
    "quantile_calibration",
    "ece_quantile",
    "sharpness",
    "central_interval",
    "tercile_groups",
    "conditional_coverage",
    "ParamScore",
    "score_quantiles",
    "format_scorecard",
]


# --------------------------------------------------------------------------- #
# Proper scoring primitives
# --------------------------------------------------------------------------- #
def pinball_loss(y_true, q_pred, q) -> np.ndarray:
    """Quantile (pinball) loss for a single quantile level ``q`` in (0, 1).

    Elementwise; broadcast-compatible. Lower is better.
    """
    y_true = np.asarray(y_true, dtype=float)
    q_pred = np.asarray(q_pred, dtype=float)
    diff = y_true - q_pred
    return np.maximum(q * diff, (q - 1.0) * diff)


def interval_score(lo, hi, y, alpha: float) -> np.ndarray:
    """Interval score for a central (1 - alpha) prediction interval [lo, hi].

    Proper scoring rule (Gneiting & Raftery 2007): width plus a miss penalty
    of (2/alpha) times the shortfall on whichever side is violated. Lower is
    better. Elementwise / broadcastable.
    """
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    y = np.asarray(y, dtype=float)
    width = hi - lo
    below = (2.0 / alpha) * np.maximum(lo - y, 0.0)
    above = (2.0 / alpha) * np.maximum(y - hi, 0.0)
    return width + below + above


# --------------------------------------------------------------------------- #
# Coverage & calibration
# --------------------------------------------------------------------------- #
def empirical_coverage(y_true, lo, hi) -> float:
    """Fraction of ``y_true`` falling within [lo, hi] (inclusive)."""
    y_true = np.asarray(y_true, dtype=float)
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    return float(np.mean((y_true >= lo) & (y_true <= hi)))


def quantile_calibration(y_true, q_pred_level, q_level: float) -> float:
    """Empirical CDF check for one quantile level.

    Returns the fraction of observations at or below the predicted ``q_level``
    quantile. A perfectly calibrated quantile gives a value equal to
    ``q_level``.
    """
    y_true = np.asarray(y_true, dtype=float)
    q_pred_level = np.asarray(q_pred_level, dtype=float)
    return float(np.mean(y_true <= q_pred_level))


def ece_quantile(y_true, q_pred, q_levels) -> float:
    """Expected calibration error over a set of quantile levels.

    For each level we compute |empirical_below - nominal| and average. This is
    the standard quantile-regression calibration error; 0 is perfect.

    Parameters
    ----------
    y_true : (n,)
    q_pred : (n, n_levels)
    q_levels : (n_levels,)
    """
    y_true = np.asarray(y_true, dtype=float)
    q_pred = np.asarray(q_pred, dtype=float)
    q_levels = np.asarray(q_levels, dtype=float)
    gaps = []
    for j, ql in enumerate(q_levels):
        emp = np.mean(y_true <= q_pred[:, j])
        gaps.append(abs(emp - ql))
    return float(np.mean(gaps))


def sharpness(lo, hi) -> float:
    """Mean width of the prediction intervals (smaller = sharper)."""
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    return float(np.mean(hi - lo))


# --------------------------------------------------------------------------- #
# Interval extraction & conditioning
# --------------------------------------------------------------------------- #
def central_interval(q_pred, q_levels, alpha: float):
    """Extract a central (1 - alpha) interval from a quantile array.

    Picks the predicted quantiles nearest to alpha/2 and 1 - alpha/2.

    Parameters
    ----------
    q_pred : (..., n_levels)
    q_levels : (n_levels,) ascending
    alpha : nominal miss rate (e.g. 0.1 for a 90% interval)

    Returns
    -------
    lo, hi : arrays with the trailing ``n_levels`` axis removed.
    """
    q_levels = np.asarray(q_levels, dtype=float)
    q_pred = np.asarray(q_pred, dtype=float)
    lo_target = alpha / 2.0
    hi_target = 1.0 - alpha / 2.0
    j_lo = int(np.argmin(np.abs(q_levels - lo_target)))
    j_hi = int(np.argmin(np.abs(q_levels - hi_target)))
    return q_pred[..., j_lo], q_pred[..., j_hi]


def tercile_groups(x) -> np.ndarray:
    """Assign each element of ``x`` to a tercile (0, 1, 2) by its value.

    Useful for probing conditional coverage across, e.g., low/mid/high D*.
    """
    x = np.asarray(x, dtype=float)
    q1, q2 = np.quantile(x, [1.0 / 3.0, 2.0 / 3.0])
    groups = np.zeros(x.shape, dtype=int)
    groups[x > q1] = 1
    groups[x > q2] = 2
    return groups


def conditional_coverage(y_true, lo, hi, groups) -> dict[int, float]:
    """Empirical coverage within each group label.

    Returns a dict {group_label: coverage}.
    """
    y_true = np.asarray(y_true, dtype=float)
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    groups = np.asarray(groups)
    out: dict[int, float] = {}
    for g in np.unique(groups):
        m = groups == g
        out[int(g)] = empirical_coverage(y_true[m], lo[m], hi[m])
    return out


# --------------------------------------------------------------------------- #
# Per-parameter scorecard
# --------------------------------------------------------------------------- #
@dataclass
class ParamScore:
    """Calibration scorecard for a single output parameter."""

    name: str
    alpha: float
    coverage: float
    ece: float
    sharpness: float
    mean_pinball: float
    mean_interval_score: float
    conditional: dict[int, float] = field(default_factory=dict)

    @property
    def nominal(self) -> float:
        """Nominal central coverage, ``1 - alpha``."""
        return 1.0 - self.alpha

    @property
    def coverage_gap(self) -> float:
        """Signed marginal coverage gap (empirical - nominal)."""
        return self.coverage - self.nominal


def score_quantiles(
    y_true,
    q_pred,
    q_levels,
    alpha: float = 0.1,
    param_names: Optional[Sequence[str]] = None,
    conditioning: Optional[np.ndarray] = None,
) -> list[ParamScore]:
    """Score predicted quantiles against truth, per parameter.

    Parameters
    ----------
    y_true : (n, n_params)
    q_pred : (n, n_params, n_levels)
    q_levels : (n_levels,) ascending in (0, 1)
    alpha : central-interval miss rate for coverage/sharpness/interval-score.
    param_names : optional list of length n_params.
    conditioning : optional (n, n_params) or (n,) array; terciles of it define
        the conditional-coverage groups. If 1-D, the same grouping is used for
        every parameter. If None, conditioning uses each parameter's own truth.

    Returns
    -------
    list of ParamScore, one per parameter.
    """
    y_true = np.asarray(y_true, dtype=float)
    q_pred = np.asarray(q_pred, dtype=float)
    q_levels = np.asarray(q_levels, dtype=float)
    if y_true.ndim != 2:
        raise ValueError("y_true must be (n, n_params)")
    if q_pred.ndim != 3:
        raise ValueError("q_pred must be (n, n_params, n_levels)")
    n, n_params = y_true.shape
    if q_pred.shape[:2] != (n, n_params):
        raise ValueError("q_pred and y_true disagree on (n, n_params)")
    if q_pred.shape[2] != q_levels.shape[0]:
        raise ValueError("q_pred n_levels != len(q_levels)")
    if not np.all(np.diff(q_levels) > 0):
        raise ValueError("q_levels must be strictly ascending")
    if param_names is None:
        param_names = [f"p{j}" for j in range(n_params)]

    scores: list[ParamScore] = []
    for j in range(n_params):
        yj = y_true[:, j]
        qj = q_pred[:, j, :]  # (n, n_levels)
        lo, hi = central_interval(qj, q_levels, alpha)
        cov = empirical_coverage(yj, lo, hi)
        ece = ece_quantile(yj, qj, q_levels)
        shp = sharpness(lo, hi)
        # mean pinball averaged over all levels
        pin = float(
            np.mean([np.mean(pinball_loss(yj, qj[:, k], q_levels[k]))
                     for k in range(len(q_levels))])
        )
        iscore = float(np.mean(interval_score(lo, hi, yj, alpha)))

        if conditioning is None:
            cond_src = yj
        else:
            cond = np.asarray(conditioning, dtype=float)
            cond_src = cond[:, j] if cond.ndim == 2 else cond
        groups = tercile_groups(cond_src)
        cond_cov = conditional_coverage(yj, lo, hi, groups)

        scores.append(
            ParamScore(
                name=param_names[j],
                alpha=alpha,
                coverage=cov,
                ece=ece,
                sharpness=shp,
                mean_pinball=pin,
                mean_interval_score=iscore,
                conditional=cond_cov,
            )
        )
    return scores


def format_scorecard(scores: Sequence[ParamScore], title: str = "scorecard") -> str:
    """Render a list of ParamScore as a fixed-width text table."""
    lines = []
    lines.append(f"=== {title} ===")
    nominal = 1.0 - scores[0].alpha if scores else 0.0
    lines.append(f"nominal central coverage = {nominal:.3f} (alpha={scores[0].alpha:.3f})"
                 if scores else "(empty)")
    header = (f"{'param':>8} {'cover':>7} {'gap':>8} {'ECE':>7} "
              f"{'sharp':>9} {'pinball':>9} {'intvl':>9}")
    lines.append(header)
    lines.append("-" * len(header))
    for s in scores:
        lines.append(
            f"{s.name:>8} {s.coverage:>7.3f} {s.coverage_gap:>+8.3f} "
            f"{s.ece:>7.3f} {s.sharpness:>9.4g} {s.mean_pinball:>9.4g} "
            f"{s.mean_interval_score:>9.4g}"
        )
    # conditional coverage block
    lines.append("")
    lines.append("conditional coverage by tercile (group 0=low .. 2=high):")
    for s in scores:
        cc = " ".join(f"g{g}={v:.3f}" for g, v in sorted(s.conditional.items()))
        lines.append(f"{s.name:>8}  {cc}")
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# Self-test + demo
# --------------------------------------------------------------------------- #
def _self_test() -> None:
    """Deterministic invariants. Raises AssertionError on failure."""
    rng = np.random.default_rng(0)

    # pinball: known values
    assert np.isclose(pinball_loss(1.0, 0.0, 0.5), 0.5)
    assert np.isclose(pinball_loss(0.0, 1.0, 0.5), 0.5)
    assert np.isclose(pinball_loss(2.0, 0.0, 0.9), 0.9 * 2.0)
    assert np.isclose(pinball_loss(0.0, 2.0, 0.9), 0.1 * 2.0)

    # interval score: inside -> width; outside -> width + penalty
    assert np.isclose(interval_score(0.0, 2.0, 1.0, alpha=0.1), 2.0)
    assert np.isclose(interval_score(0.0, 2.0, 3.0, alpha=0.1), 2.0 + 20.0)
    assert np.isclose(interval_score(0.0, 2.0, -1.0, alpha=0.1), 2.0 + 20.0)

    # empirical coverage trivial cases
    y = np.array([0.0, 1.0, 2.0])
    assert empirical_coverage(y, np.full(3, -1.0), np.full(3, 3.0)) == 1.0
    assert empirical_coverage(y, np.full(3, 5.0), np.full(3, 6.0)) == 0.0

    # central_interval picks correct columns
    qlev = np.array([0.05, 0.5, 0.95])
    qp = np.array([[1.0, 2.0, 3.0]])
    lo, hi = central_interval(qp, qlev, alpha=0.1)
    assert np.isclose(lo[0], 1.0) and np.isclose(hi[0], 3.0)

    # well-specified Gaussian: predicted true quantiles -> coverage ~ nominal,
    # ECE ~ 0. This is the core calibration invariant.
    n = 20000
    mu = rng.normal(size=(n, 1))
    y_true = mu + rng.normal(size=(n, 1))  # sigma = 1 around mu
    levels = np.array([0.05, 0.25, 0.5, 0.75, 0.95])
    from math import sqrt  # noqa
    from statistics import NormalDist
    z = np.array([NormalDist().inv_cdf(p) for p in levels])
    q_pred = mu[:, :, None] + z[None, None, :]  # (n,1,n_levels)
    # Condition on the covariate mu (independent of the noise), NOT on the
    # label y: a well-specified model has ~nominal coverage in every covariate
    # tercile. (Grouping by the realized outcome would induce a selection
    # artifact -- a real effect we exploit deliberately in the IVIM probe.)
    scores = score_quantiles(y_true, q_pred, levels, alpha=0.1,
                             param_names=["gauss"], conditioning=mu)
    s = scores[0]
    assert abs(s.coverage - 0.9) < 0.02, f"coverage {s.coverage}"
    assert s.ece < 0.02, f"ece {s.ece}"
    for g, v in s.conditional.items():
        assert abs(v - 0.9) < 0.05, f"cond g{g}={v}"

    # tercile grouping balanced
    g = tercile_groups(rng.normal(size=9000))
    counts = np.bincount(g)
    assert np.all(np.abs(counts - 3000) < 200)

    print("[metrics self-test] all invariants passed")


def _demo() -> None:
    """Print a scorecard for a deliberately over-confident estimator."""
    rng = np.random.default_rng(1)
    n = 5000
    mu = rng.normal(size=(n, 2))
    # true noise sigma = 1, but the "model" believes sigma = 0.5 -> overconfident
    y_true = mu + rng.normal(size=(n, 2))
    levels = np.array([0.05, 0.25, 0.5, 0.75, 0.95])
    from statistics import NormalDist
    z = np.array([NormalDist().inv_cdf(p) for p in levels])
    q_pred = mu[:, :, None] + 0.5 * z[None, None, :]
    scores = score_quantiles(y_true, q_pred, levels, alpha=0.1,
                             param_names=["alpha", "beta"])
    print(format_scorecard(scores, title="demo: over-confident estimator"))
    print("\n(note: coverage < nominal and gap < 0 -> intervals too tight)")


if __name__ == "__main__":
    _self_test()
    print()
    _demo()
