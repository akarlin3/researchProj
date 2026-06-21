"""Calibration ruler: coverage / ECE / sharpness, re-derived independently.

Standard, published quantities, re-derived here from their definitions -- NOT
imported from ``caliper.metrics`` (the forbidden module). Estimators expose a single
contract::

    predict_quantiles(signals, q_levels) -> q_pred,  shape (n, n_params, n_levels)

with ``q_levels`` ascending in (0, 1). For a central interval at nominal level
``1 - alpha`` the bounds are the ``alpha/2`` and ``1 - alpha/2`` predicted quantiles.

Definitions:
* **coverage**  = mean( lo <= y <= hi )
* **ECE**       = mean over the symmetric central levels of |empirical_coverage - nominal|
                  (Fashion's definition: "mean |coverage - nominal| across levels").
* **sharpness** = mean interval width hi - lo (must be read alongside coverage).
* **pinball / interval score** = proper-scoring cross-checks (Gneiting & Raftery 2007).
* **conditional coverage** by D*-tercile exposes the high-D* under-coverage.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np


def _idx(q_levels, level, tol=1e-9):
    q = np.asarray(q_levels, dtype=float)
    hits = np.where(np.abs(q - level) < tol)[0]
    if len(hits) == 0:
        raise ValueError(f"level {level} not in q_levels {q}")
    return int(hits[0])


def central_levels(q_levels):
    """Symmetric (lo, hi, nominal) triples from a symmetric quantile grid.

    e.g. levels {0.05,.1,.25,.5,.75,.9,.95} -> (0.05,0.95,0.90), (0.1,0.9,0.8),
    (0.25,0.75,0.50).
    """
    q = sorted(set(np.round(np.asarray(q_levels, dtype=float), 12)))
    out = []
    for lo in q:
        hi = round(1.0 - lo, 12)
        if lo < 0.5 and hi in q:
            out.append((lo, hi, round(hi - lo, 12)))
    return sorted(out, key=lambda t: t[2])


def empirical_coverage(y, lo, hi):
    y = np.asarray(y, dtype=float)
    return float(np.mean((y >= np.asarray(lo)) & (y <= np.asarray(hi))))


def sharpness(lo, hi):
    return float(np.mean(np.asarray(hi) - np.asarray(lo)))


def pinball_loss(y, q_pred_level, tau):
    y = np.asarray(y, dtype=float); q = np.asarray(q_pred_level, dtype=float)
    d = y - q
    return float(np.mean(np.maximum(tau * d, (tau - 1.0) * d)))


def interval_score(y, lo, hi, alpha):
    """Gneiting & Raftery (2007) interval score (lower is better)."""
    y = np.asarray(y, dtype=float); lo = np.asarray(lo); hi = np.asarray(hi)
    width = hi - lo
    below = (2.0 / alpha) * (lo - y) * (y < lo)
    above = (2.0 / alpha) * (y - hi) * (y > hi)
    return float(np.mean(width + below + above))


def tercile_groups(x):
    """Assign each element of ``x`` to a tercile label {0,1,2} by its own quantiles."""
    x = np.asarray(x, dtype=float)
    e1, e2 = np.quantile(x, [1 / 3, 2 / 3])
    g = np.zeros(len(x), dtype=int)
    g[x > e1] = 1
    g[x > e2] = 2
    return g


def conditional_coverage(y, lo, hi, groups):
    y = np.asarray(y, dtype=float)
    out = {}
    for gv in sorted(set(np.asarray(groups).tolist())):
        m = np.asarray(groups) == gv
        if m.any():
            out[int(gv)] = empirical_coverage(y[m], np.asarray(lo)[m], np.asarray(hi)[m])
    return out


@dataclass
class ParamScore:
    name: str
    alpha: float                 # central-interval miss rate (nominal = 1 - alpha)
    coverage: float
    ece: float
    sharpness: float
    mean_pinball: float
    mean_interval_score: float
    conditional: dict = field(default_factory=dict)

    @property
    def nominal(self):
        return 1.0 - self.alpha

    @property
    def coverage_gap(self):
        return self.coverage - self.nominal


def score_quantiles(y_true, q_pred, q_levels, alpha=0.05, param_names=None,
                    conditioning=None):
    """Score predicted quantiles for each parameter.

    Parameters
    ----------
    y_true : (n, P) ground truth.
    q_pred : (n, P, L) predicted quantiles at ``q_levels``.
    q_levels : (L,) ascending in (0, 1); must contain ``alpha/2`` and ``1-alpha/2``.
    alpha : central-interval miss rate (default 0.05 -> nominal 0.95, Fashion's headline).
    conditioning : optional (n, P) or (n,) values to tercile-stratify coverage by
        (e.g. true D*); if (n,P), column p conditions parameter p.
    """
    y_true = np.asarray(y_true, dtype=float)
    q_pred = np.asarray(q_pred, dtype=float)
    q_levels = np.asarray(q_levels, dtype=float)
    n, P = y_true.shape
    names = list(param_names) if param_names else [f"p{p}" for p in range(P)]
    lo_i, hi_i = _idx(q_levels, alpha / 2), _idx(q_levels, 1 - alpha / 2)
    clev = central_levels(q_levels)

    scores = []
    for p in range(P):
        y = y_true[:, p]
        lo, hi = q_pred[:, p, lo_i], q_pred[:, p, hi_i]
        # ECE: mean |empirical central coverage - nominal| across symmetric levels.
        gaps = []
        for lvl_lo, lvl_hi, nom in clev:
            cl, ch = q_pred[:, p, _idx(q_levels, lvl_lo)], q_pred[:, p, _idx(q_levels, lvl_hi)]
            gaps.append(abs(empirical_coverage(y, cl, ch) - nom))
        ece = float(np.mean(gaps)) if gaps else float("nan")
        pin = float(np.mean([pinball_loss(y, q_pred[:, p, j], q_levels[j])
                             for j in range(len(q_levels))]))
        cond = {}
        if conditioning is not None:
            cvals = np.asarray(conditioning, dtype=float)
            cvals = cvals[:, p] if cvals.ndim == 2 else cvals
            cond = conditional_coverage(y, lo, hi, tercile_groups(cvals))
        scores.append(ParamScore(
            name=names[p], alpha=float(alpha),
            coverage=empirical_coverage(y, lo, hi), ece=ece,
            sharpness=sharpness(lo, hi), mean_pinball=pin,
            mean_interval_score=interval_score(y, lo, hi, alpha),
            conditional=cond))
    return scores


def format_scorecard(scores):
    rows = ["param   nominal  coverage   gap     ece    sharpness",
            "-----   -------  --------  ------  ------  ---------"]
    for s in scores:
        rows.append(f"{s.name:6s}  {s.nominal:6.2f}  {s.coverage:8.3f}  "
                    f"{s.coverage_gap:+6.3f}  {s.ece:6.3f}  {s.sharpness:9.4g}")
    return "\n".join(rows)
