"""Caliper ruler demo -- the model-agnostic calibration ruler, numpy only.

Run:  python examples/ruler_demo.py      (no torch, no IVIM required)

The ruler (:mod:`caliper.metrics`) scores *any* estimator's predicted quantiles
against truth -- it knows nothing about IVIM, flows, or conformal prediction.
Feed it ``(y_true, q_pred, q_levels)`` and it returns coverage, calibration
error (ECE), sharpness, and conditional (group-wise) coverage.

This script shows the ruler distinguishing a *well-specified* estimator (reports
its true spread) from an *over-confident* one (reports too-narrow quantiles) on
the same toy target. Every number printed is produced by this fixed-seed run.
"""
from __future__ import annotations

from statistics import NormalDist

import numpy as np

from caliper import metrics as M

LEVELS = np.array([0.05, 0.25, 0.5, 0.75, 0.95])
ALPHA = 0.10  # 90% central intervals
_Z = np.array([NormalDist().inv_cdf(p) for p in LEVELS])


def _gaussian_quantiles(mu, sigma):
    """Predicted quantiles for a Gaussian estimator: (n, 1, n_levels)."""
    return mu[:, :, None] + sigma * _Z[None, None, :]


def main() -> None:
    print("Caliper ruler demo -- model-agnostic calibration scoring (numpy only)\n")
    rng = np.random.default_rng(0)
    n = 20000

    # A covariate mu and an outcome y with true conditional spread sigma = 1.0.
    mu = rng.normal(size=(n, 1))
    y_true = mu + rng.normal(0.0, 1.0, size=(n, 1))

    # Two estimators of the SAME outcome, scored by the SAME ruler:
    #   * well-specified: reports the true sigma = 1.0
    #   * over-confident: reports sigma = 0.5 (quantiles too narrow)
    cases = {
        "well-specified (sigma=1.0)": _gaussian_quantiles(mu, 1.0),
        "over-confident (sigma=0.5)": _gaussian_quantiles(mu, 0.5),
    }
    for title, q_pred in cases.items():
        # Condition conditional-coverage on the covariate mu (independent of the
        # noise): a well-specified model is ~nominal in every covariate tercile.
        scores = M.score_quantiles(y_true, q_pred, LEVELS, alpha=ALPHA,
                                   param_names=["y"], conditioning=mu)
        print(M.format_scorecard(scores, title=title))
        print()

    print("Reading it: the well-specified estimator lands at coverage ~0.90 with "
          "ECE ~0;\nthe over-confident one under-covers (coverage < 0.90, gap < 0) "
          "because its\nreported quantiles are too narrow -- exactly the failure "
          "conformal prediction\nrepairs (see examples/conformal_demo.py). Same "
          "ruler, any estimator.")


if __name__ == "__main__":
    main()
