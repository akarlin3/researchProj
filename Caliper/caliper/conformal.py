"""caliper.conformal -- split-conformal / CQR coverage correction.

Numpy-only and estimator-agnostic. Works on the quantile arrays produced by any
estimator exposing ``predict_quantiles(signals, q_levels) -> (n, n_params,
n_levels)``.

The method is conformalized quantile regression (CQR; Romano, Patterson &
Candes 2019) applied independently to each symmetric pair of quantile levels.
For a central pair (q_lo, q_hi) with miss-rate ``a = 2 * q_lo``, the calibration
conformity score is

    E_i = max(q_lo(x_i) - y_i,  y_i - q_hi(x_i)),

and the correction is the finite-sample-adjusted (1 - a) empirical quantile of
the {E_i}. The corrected interval is [q_lo - Q, q_hi + Q]. This restores
*marginal* coverage to nominal under exchangeability; it does not by itself fix
*conditional* coverage.

This module also provides three companions to plain CQR:

* ``SplitConformalResidual`` -- the simplest split-conformal interval, built from
  an estimator's *point* prediction and the absolute-residual nonconformity
  ``s_i = |y_i - point_i|``; the correction is the ``ceil((n+1)(1-a))/n``
  empirical quantile of ``{s_i}`` and the interval is ``point +/- Q``.
* ``MondrianConformalQuantile`` -- group-conditional CQR: an independent CQR
  correction is fit *within each caller-supplied group*, buying a per-group
  (Mondrian) coverage guarantee that plain CQR's marginal guarantee does not.
* ``conditional_coverage_by_strata`` -- a diagnostic returning empirical coverage
  *and mean interval width* within each stratum, the lens for reading where
  marginal calibration hides conditional miscoverage and at what width cost.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .metrics import empirical_coverage

__all__ = [
    "conformity_scores",
    "conformal_offset",
    "SplitConformalQuantile",
    "SplitConformalResidual",
    "MondrianConformalQuantile",
    "StratumCoverage",
    "conditional_coverage_by_strata",
    "format_strata_table",
]


def conformity_scores(q_lo, q_hi, y) -> np.ndarray:
    """CQR conformity scores E_i = max(q_lo - y, y - q_hi)."""
    q_lo = np.asarray(q_lo, dtype=float)
    q_hi = np.asarray(q_hi, dtype=float)
    y = np.asarray(y, dtype=float)
    return np.maximum(q_lo - y, y - q_hi)


def conformal_offset(scores, alpha: float) -> float:
    """Finite-sample (1 - alpha) quantile of the conformity scores.

    Uses the standard CQR level ceil((n+1)(1-alpha))/n, clipped to [0, 1].
    """
    scores = np.asarray(scores, dtype=float)
    n = scores.shape[0]
    if n == 0:
        return 0.0
    level = np.ceil((n + 1) * (1.0 - alpha)) / n
    level = min(max(level, 0.0), 1.0)
    return float(np.quantile(scores, level, method="higher"))


def _symmetric_pairs(n_levels: int):
    """Yield (j_lo, j_hi) index pairs symmetric about the centre."""
    for j in range(n_levels // 2):
        yield j, n_levels - 1 - j


@dataclass
class SplitConformalQuantile:
    """Split-conformal / CQR wrapper over an estimator's quantile output.

    Usage
    -----
    >>> cq = SplitConformalQuantile(q_levels)
    >>> cq.calibrate(q_pred_cal, y_cal)      # q_pred_cal: (n_cal, P, L)
    >>> q_corr = cq.apply(q_pred_test)       # (n_test, P, L)

    The correction is per-parameter and per symmetric level-pair. The median
    level (if present) is left unchanged. Output quantiles are re-sorted along
    the level axis to guarantee monotonicity.
    """

    q_levels: np.ndarray

    def __post_init__(self) -> None:
        self.q_levels = np.asarray(self.q_levels, dtype=float)
        if not np.all(np.diff(self.q_levels) > 0):
            raise ValueError("q_levels must be strictly ascending")
        # offsets_[p][(j_lo, j_hi)] = Q ; filled by calibrate()
        self.offsets_: list[dict[tuple[int, int], float]] = []
        self._n_params: int | None = None

    def calibrate(self, q_pred_cal, y_cal) -> "SplitConformalQuantile":
        """Fit per-parameter, per-level-pair CQR offsets on a calibration split.

        ``q_pred_cal`` is ``(n_cal, P, L)`` raw quantiles; ``y_cal`` is
        ``(n_cal, P)`` truth. Returns ``self``.
        """
        q_pred_cal = np.asarray(q_pred_cal, dtype=float)
        y_cal = np.asarray(y_cal, dtype=float)
        n, P, L = q_pred_cal.shape
        if y_cal.shape != (n, P):
            raise ValueError("y_cal must be (n_cal, n_params)")
        if L != self.q_levels.shape[0]:
            raise ValueError("q_pred_cal level axis != len(q_levels)")
        self._n_params = P
        self.offsets_ = []
        for p in range(P):
            off: dict[tuple[int, int], float] = {}
            for j_lo, j_hi in _symmetric_pairs(L):
                a = 2.0 * self.q_levels[j_lo]  # nominal miss-rate of this pair
                scores = conformity_scores(
                    q_pred_cal[:, p, j_lo], q_pred_cal[:, p, j_hi], y_cal[:, p]
                )
                off[(j_lo, j_hi)] = conformal_offset(scores, a)
            self.offsets_.append(off)
        return self

    def apply(self, q_pred) -> np.ndarray:
        """Apply the fitted offsets to ``(n, P, L)`` quantiles; re-sorted monotone."""
        if not self.offsets_:
            raise RuntimeError("call calibrate() before apply()")
        q_pred = np.asarray(q_pred, dtype=float)
        n, P, L = q_pred.shape
        if P != self._n_params:
            raise ValueError("n_params mismatch with calibration")
        out = q_pred.copy()
        for p in range(P):
            for (j_lo, j_hi), Q in self.offsets_[p].items():
                out[:, p, j_lo] = q_pred[:, p, j_lo] - Q
                out[:, p, j_hi] = q_pred[:, p, j_hi] + Q
        # enforce monotonic, non-crossing quantiles along the level axis
        out = np.sort(out, axis=2)
        return out

    def calibrate_apply(self, q_pred_cal, y_cal, q_pred_test) -> np.ndarray:
        """Convenience: calibrate on one split, correct another."""
        return self.calibrate(q_pred_cal, y_cal).apply(q_pred_test)


@dataclass
class SplitConformalResidual:
    """Plain split-conformal intervals from a point predictor.

    Nonconformity is the absolute residual ``s_i = |y_i - point_i|`` (per
    parameter); the correction ``Q`` is the ``ceil((n+1)(1-alpha))/n`` empirical
    quantile of the calibration scores and the prediction interval is
    ``point +/- Q``. Symmetric and width-constant across the test set -- the
    baseline conformal method against which CQR's input-adaptive widths are
    judged.

    Usage
    -----
    >>> sc = SplitConformalResidual(alpha=0.1)
    >>> sc.calibrate(point_cal, y_cal)        # both (n_cal, P)
    >>> lo, hi = sc.apply(point_test)         # both (n_test, P)
    """

    alpha: float = 0.1

    def __post_init__(self) -> None:
        if not 0.0 < self.alpha < 1.0:
            raise ValueError("alpha must lie in (0, 1)")
        self.offsets_: np.ndarray | None = None  # (P,)

    def calibrate(self, point_cal, y_cal) -> "SplitConformalResidual":
        """Fit per-parameter absolute-residual offsets on a calibration split.

        ``point_cal`` and ``y_cal`` are both ``(n_cal, P)``. Returns ``self``.
        """
        point_cal = np.asarray(point_cal, dtype=float)
        y_cal = np.asarray(y_cal, dtype=float)
        if point_cal.shape != y_cal.shape or point_cal.ndim != 2:
            raise ValueError("point_cal and y_cal must share shape (n_cal, n_params)")
        scores = np.abs(y_cal - point_cal)  # (n, P)
        self.offsets_ = np.array(
            [conformal_offset(scores[:, p], self.alpha) for p in range(scores.shape[1])]
        )
        return self

    def apply(self, point):
        """Return ``(lo, hi)`` = ``point -/+ Q`` for ``(n, P)`` point predictions."""
        if self.offsets_ is None:
            raise RuntimeError("call calibrate() before apply()")
        point = np.asarray(point, dtype=float)
        if point.ndim != 2 or point.shape[1] != self.offsets_.shape[0]:
            raise ValueError("point must be (n_test, n_params) matching calibration")
        lo = point - self.offsets_[None, :]
        hi = point + self.offsets_[None, :]
        return lo, hi

    def calibrate_apply(self, point_cal, y_cal, point_test):
        """Convenience: calibrate on one split, return (lo, hi) for another."""
        return self.calibrate(point_cal, y_cal).apply(point_test)


@dataclass
class MondrianConformalQuantile:
    """Group-conditional (Mondrian) CQR.

    Fits an independent :class:`SplitConformalQuantile` correction *within each
    group label* supplied at calibration time, so the (1 - alpha) coverage
    guarantee holds conditionally on group membership rather than only
    marginally. A pooled (all-data) CQR is also fit as a fallback for any group
    unseen during calibration.

    Usage
    -----
    >>> mq = MondrianConformalQuantile(q_levels)
    >>> mq.calibrate(q_pred_cal, y_cal, groups_cal)   # groups_cal: (n_cal,)
    >>> q_corr = mq.apply(q_pred_test, groups_test)    # groups_test: (n_test,)

    The price of the per-group guarantee is per-group width: groups whose true
    error is large (e.g. the high-D* tercile of IVIM) are corrected with a
    larger offset, so Mondrian restores their coverage only by inflating their
    intervals. ``conditional_coverage_by_strata`` is the tool for reading that
    trade.
    """

    q_levels: np.ndarray

    def __post_init__(self) -> None:
        self.q_levels = np.asarray(self.q_levels, dtype=float)
        if not np.all(np.diff(self.q_levels) > 0):
            raise ValueError("q_levels must be strictly ascending")
        self.group_models_: dict[int, SplitConformalQuantile] = {}
        self.global_: SplitConformalQuantile | None = None

    def calibrate(self, q_pred_cal, y_cal, groups_cal) -> "MondrianConformalQuantile":
        """Fit an independent CQR correction within each calibration group label.

        ``groups_cal`` is ``(n_cal,)`` integer labels (e.g. true-D* terciles). A
        pooled CQR is also fit as a fallback for groups unseen at test time.
        Returns ``self``.
        """
        q_pred_cal = np.asarray(q_pred_cal, dtype=float)
        y_cal = np.asarray(y_cal, dtype=float)
        groups_cal = np.asarray(groups_cal)
        if groups_cal.shape[0] != q_pred_cal.shape[0]:
            raise ValueError("groups_cal must have one label per calibration row")
        self.group_models_ = {}
        for g in np.unique(groups_cal):
            m = groups_cal == g
            self.group_models_[int(g)] = SplitConformalQuantile(self.q_levels).calibrate(
                q_pred_cal[m], y_cal[m]
            )
        self.global_ = SplitConformalQuantile(self.q_levels).calibrate(q_pred_cal, y_cal)
        return self

    def apply(self, q_pred, groups) -> np.ndarray:
        """Correct ``(n, P, L)`` quantiles per test-row group (pooled fallback)."""
        if self.global_ is None:
            raise RuntimeError("call calibrate() before apply()")
        q_pred = np.asarray(q_pred, dtype=float)
        groups = np.asarray(groups)
        if groups.shape[0] != q_pred.shape[0]:
            raise ValueError("groups must have one label per test row")
        out = q_pred.copy()
        for g in np.unique(groups):
            m = groups == g
            model = self.group_models_.get(int(g), self.global_)
            out[m] = model.apply(q_pred[m])
        return out

    def calibrate_apply(self, q_pred_cal, y_cal, groups_cal, q_pred_test, groups_test):
        """Convenience: calibrate per group on one split, correct another."""
        return self.calibrate(q_pred_cal, y_cal, groups_cal).apply(q_pred_test, groups_test)


@dataclass
class StratumCoverage:
    """Coverage and mean interval width within a single stratum."""

    stratum: int
    n: int
    coverage: float
    mean_width: float


def conditional_coverage_by_strata(y_true, lower, upper, strata) -> dict[int, StratumCoverage]:
    """Empirical coverage *and* mean interval width within each stratum.

    Parameters
    ----------
    y_true, lower, upper : 1-D arrays of equal length (one output parameter).
    strata : integer group label per observation (e.g. ``metrics.tercile_groups``
        of a covariate such as true D*).

    Returns
    -------
    dict mapping stratum label -> :class:`StratumCoverage`. Coverage reuses
    ``metrics.empirical_coverage`` (the canonical ruler); width is the mean of
    ``upper - lower`` within the stratum.
    """
    y_true = np.asarray(y_true, dtype=float)
    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    strata = np.asarray(strata)
    if not (y_true.shape == lower.shape == upper.shape == strata.shape):
        raise ValueError("y_true, lower, upper, strata must share shape (n,)")
    out: dict[int, StratumCoverage] = {}
    for g in np.unique(strata):
        m = strata == g
        out[int(g)] = StratumCoverage(
            stratum=int(g),
            n=int(m.sum()),
            coverage=empirical_coverage(y_true[m], lower[m], upper[m]),
            mean_width=float(np.mean(upper[m] - lower[m])),
        )
    return out


def format_strata_table(
    per_method: dict[str, dict[int, StratumCoverage]],
    stratum_names: dict[int, str] | None = None,
    title: str = "conditional coverage by stratum",
    nominal: float | None = None,
) -> str:
    """Render a {method -> {stratum -> StratumCoverage}} mapping as a text table.

    Columns are coverage and mean width per stratum; rows are methods. Width
    ratios (high/low stratum) are appended per method when there are >=2 strata.
    """
    methods = list(per_method)
    strata = sorted({s for d in per_method.values() for s in d})
    names = stratum_names or {s: f"g{s}" for s in strata}

    lines = [f"=== {title} ==="]
    if nominal is not None:
        lines.append(f"nominal coverage = {nominal:.3f}")
    head = f"{'method':>22} | " + " | ".join(
        f"{names[s] + ' cov':>9} {'width':>9}" for s in strata
    )
    lines.append(head)
    lines.append("-" * len(head))
    for meth in methods:
        d = per_method[meth]
        cells = []
        for s in strata:
            sc = d.get(s)
            if sc is None:
                cells.append(f"{'--':>9} {'--':>9}")
            else:
                cells.append(f"{sc.coverage:>9.3f} {sc.mean_width:>9.4g}")
        lines.append(f"{meth:>22} | " + " | ".join(cells))
    # width ratio high/low per method
    if len(strata) >= 2:
        lo_s, hi_s = strata[0], strata[-1]
        lines.append("")
        lines.append(f"width ratio ({names[hi_s]} / {names[lo_s]}):")
        for meth in methods:
            d = per_method[meth]
            if lo_s in d and hi_s in d and d[lo_s].mean_width > 0:
                ratio = d[hi_s].mean_width / d[lo_s].mean_width
                lines.append(f"{meth:>22}  {ratio:>6.2f}x")
    return "\n".join(lines)


if __name__ == "__main__":
    # Self-contained sanity demo: an over-confident Gaussian estimator that
    # conformal restores to nominal marginal coverage.
    from statistics import NormalDist

    from caliper import metrics as M

    rng = np.random.default_rng(0)
    levels = np.array([0.05, 0.25, 0.5, 0.75, 0.95])
    z = np.array([NormalDist().inv_cdf(p) for p in levels])

    def make(n):
        mu = rng.normal(size=(n, 1))
        y = mu + rng.normal(size=(n, 1))           # true sigma = 1
        q = mu[:, :, None] + 0.5 * z[None, None, :]  # believes sigma = 0.5
        return y, q

    y_cal, q_cal = make(2000)
    y_te, q_te = make(4000)
    raw = M.score_quantiles(y_te, q_te, levels, alpha=0.1)[0]
    cq = SplitConformalQuantile(levels).calibrate(q_cal, y_cal)
    q_corr = cq.apply(q_te)
    cor = M.score_quantiles(y_te, q_corr, levels, alpha=0.1)[0]
    print(f"raw  coverage={raw.coverage:.3f} gap={raw.coverage_gap:+.3f}")
    print(f"conf coverage={cor.coverage:.3f} gap={cor.coverage_gap:+.3f}")
