"""Tests for the NLLS IVIM base estimator.

The headline test is the CP0 forward-model sanity gate: fitting clean
(noise-free) signals generated from known (D, D*, f) must return the truth.
If clean recovery fails, the generator/estimator is broken.
"""
import numpy as np
import pytest

from gauge.forward import ivim_signal, DEFAULT_B_VALUES
from gauge.estimators import fit_nlls, fit_nlls_batch


CLEAN_CASES = [
    (0.8e-3, 15e-3, 0.10),
    (1.5e-3, 50e-3, 0.20),
    (1.2e-3, 30e-3, 0.30),
    (2.5e-3, 90e-3, 0.35),
    (1.0e-3, 10e-3, 0.05),
]


@pytest.mark.parametrize("D,Dstar,f", CLEAN_CASES)
def test_clean_recovery_returns_truth(D, Dstar, f):
    s = ivim_signal(DEFAULT_B_VALUES, D, Dstar, f, S0=1.0)
    est = fit_nlls(s, DEFAULT_B_VALUES)
    assert est["D"] == pytest.approx(D, rel=1e-2)
    assert est["Dstar"] == pytest.approx(Dstar, rel=1e-2)
    assert est["f"] == pytest.approx(f, rel=1e-2)


def test_clean_recovery_grid_max_error_small():
    rng = np.random.default_rng(0)
    max_rel = 0.0
    for _ in range(60):
        D = rng.uniform(0.5e-3, 3.0e-3)
        Dstar = rng.uniform(10e-3, 100e-3)
        f = rng.uniform(0.05, 0.40)
        s = ivim_signal(DEFAULT_B_VALUES, D, Dstar, f, S0=1.0)
        est = fit_nlls(s, DEFAULT_B_VALUES)
        for truth, key in [(D, "D"), (Dstar, "Dstar"), (f, "f")]:
            max_rel = max(max_rel, abs(est[key] - truth) / truth)
    assert max_rel < 1e-2


def test_batch_fit_shape_and_order():
    rng = np.random.default_rng(1)
    truths = np.array([[1.0e-3, 20e-3, 0.1], [1.8e-3, 60e-3, 0.25]])
    signals = np.stack([ivim_signal(DEFAULT_B_VALUES, *t) for t in truths])
    out = fit_nlls_batch(signals, DEFAULT_B_VALUES)
    assert out.shape == (2, 3)  # columns: D, Dstar, f
    np.testing.assert_allclose(out, truths, rtol=1e-2)
