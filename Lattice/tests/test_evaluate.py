"""Calibration-evaluation interface gates (the CP3 contract)."""

import numpy as np
import pytest

from lattice import make_cohort, to_scorer_inputs, DEFAULT_QUANTILE_LEVELS
from lattice.evaluate import (
    IVIMQuantileEstimator,
    central_interval,
    interval_coverage,
    mean_sharpness,
)


class _OracleQuantiles:
    """A trivial estimator: returns truth-centred Gaussian-ish quantiles."""

    def __init__(self, truth, sd):
        self.truth = truth
        self.sd = sd

    def predict_quantiles(self, signals, q_levels):
        from math import sqrt
        # inverse-normal via erfinv-free approximation is overkill; use a linear
        # spread keyed to the level rank -- enough to exercise the contract.
        levels = np.asarray(q_levels)
        z = (levels - 0.5) * 4.0  # maps (0,1) -> roughly (-2, 2)
        n, p = self.truth.shape
        out = np.empty((n, p, len(levels)))
        for j in range(p):
            out[:, j, :] = self.truth[:, [j]] + self.sd[j] * z[None, :]
        return out


def test_protocol_runtime_checkable():
    est = _OracleQuantiles(np.zeros((3, 3)), [1, 1, 1])
    assert isinstance(est, IVIMQuantileEstimator)


def test_to_scorer_inputs_shapes():
    c = make_cohort("biexp", n=50, seed=3)
    est = _OracleQuantiles(c.params, sd=[1e-4, 5e-3, 0.03])
    q = est.predict_quantiles(c.signals, DEFAULT_QUANTILE_LEVELS)
    payload = to_scorer_inputs(c, q)
    assert payload["y_true"].shape == (50, 3)
    assert payload["q_pred"].shape == (50, 3, len(DEFAULT_QUANTILE_LEVELS))
    assert payload["param_names"] == ["D", "Dstar", "f"]
    assert len(payload["q_levels"]) == len(DEFAULT_QUANTILE_LEVELS)


def test_to_scorer_inputs_rejects_bad_shape():
    c = make_cohort("biexp", n=10)
    with pytest.raises(ValueError):
        to_scorer_inputs(c, np.zeros((10, 3, 2)))  # wrong n_levels


def test_central_interval_and_coverage():
    c = make_cohort("biexp", n=400, seed=5)
    # Wide oracle intervals -> coverage near 1; tiny -> coverage near 0.
    wide = _OracleQuantiles(c.params, sd=[1.0, 1.0, 1.0])
    q = wide.predict_quantiles(c.signals, DEFAULT_QUANTILE_LEVELS)
    lo, hi = central_interval(q, DEFAULT_QUANTILE_LEVELS, alpha=0.10)
    cov = interval_coverage(c.params, lo, hi)
    assert np.all(cov > 0.99)
    sharp = mean_sharpness(lo, hi)
    assert np.all(sharp > 0)
