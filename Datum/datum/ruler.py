"""Datum's ruler adapter -- a thin, read-only wrapper over Caliper's metrics.

Datum does not own a ruler; it scores methods on **Fashion's calibration ruler as
packaged by Caliper** (``caliper.metrics``). This module is the single chokepoint
through which all scoring flows, so that (a) the ruler version is the one pinned in
``datum.manifest``, and (b) every number that comes out is marked PROVISIONAL while
the ruler is in review.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from datum import _paths

_paths.ensure_deps()  # make caliper/gauge importable inside the monorepo

from caliper.metrics import score_quantiles  # noqa: E402  (after bootstrap)

from datum.manifest import RULER  # noqa: E402
from datum.provisional import Provisional, stamp  # noqa: E402

# IVIM parameter order is fixed by the substrate: (D, D*, f).
PARAM_NAMES = ("D", "Dstar", "f")

# The metrics Datum reports, each ruler-derived and therefore PROVISIONAL.
REPORTED_METRICS = ("coverage", "coverage_gap", "ece", "sharpness",
                    "mean_pinball", "mean_interval_score")


@dataclass
class Scorecard:
    """Per-parameter ruler scores, with every number stamped PROVISIONAL."""
    alpha: float
    nominal: float
    per_param: dict          # param_name -> {metric: Provisional}
    conditional: dict        # param_name -> {group_label: coverage}

    def headline(self) -> dict:
        """The load-bearing numbers (coverage gap per parameter)."""
        return {p: m["coverage_gap"] for p, m in self.per_param.items()}


def score(y_true, q_pred, q_levels, alpha: float = 0.10,
          conditioning=None, param_names=PARAM_NAMES) -> Scorecard:
    """Score predicted quantiles against ground truth via Caliper's ruler.

    Parameters mirror ``caliper.metrics.score_quantiles``:
      y_true       : (n, P) ground-truth parameters.
      q_pred       : (n, P, L) predicted quantiles.
      q_levels     : (L,) quantile levels in (0, 1), ascending.
      conditioning : optional (n,) or (n, P) values for per-group coverage
                     (Datum uses D* for the identifiability-wall terciles).

    Returns a ``Scorecard`` in which every metric is a ``Provisional`` -- the
    ruler is in review, so these are never final reference numbers.
    """
    y_true = np.asarray(y_true, dtype=float)
    q_pred = np.asarray(q_pred, dtype=float)
    scores = score_quantiles(y_true, q_pred, np.asarray(q_levels, dtype=float),
                             alpha=alpha, param_names=list(param_names),
                             conditioning=conditioning)
    per_param, conditional = {}, {}
    for s in scores:
        per_param[s.name] = {
            "coverage": stamp(s.coverage, "coverage"),
            "coverage_gap": stamp(s.coverage_gap, "coverage_gap"),
            "ece": stamp(s.ece, "ece"),
            "sharpness": stamp(s.sharpness, "sharpness"),
            "mean_pinball": stamp(s.mean_pinball, "mean_pinball"),
            "mean_interval_score": stamp(s.mean_interval_score, "mean_interval_score"),
        }
        conditional[s.name] = dict(s.conditional)
    nominal = 1.0 - alpha
    return Scorecard(alpha=alpha, nominal=nominal,
                     per_param=per_param, conditional=conditional)


def ruler_id() -> str:
    """Human-readable id of the pinned ruler, for provenance in outputs."""
    return (f"{RULER['name']} v{RULER['version']} @ {RULER['commit']} "
            f"({RULER['manuscript_status']})")


__all__ = ["score", "Scorecard", "ruler_id", "PARAM_NAMES", "REPORTED_METRICS",
           "Provisional"]
