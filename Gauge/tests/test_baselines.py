"""Contract + determinism smoke tests for the model-based baselines.

Kept small/fast (few epochs, tiny cohort). We do NOT assert calibration here --
calibration quality is exactly what CP1 measures honestly. We assert the
*interface*: every baseline emits finite per-parameter predictive draws of the
right shape, deterministically, and produces ordered quantile bands.
"""
import numpy as np
import pytest

from gauge.cohort import generate_cohort
from gauge.baselines import (
    SingleGaussianPNN, MDNDeepEnsemble, DeepEnsemblePoint, BayesianIVIM_MCMC,
    quantiles_from_samples, _nlls_init_and_noise,
)

pytest.importorskip("torch")


@pytest.fixture(scope="module")
def tiny():
    c = generate_cohort(200, 80, 80, snr_grid=(20, 50), seed=1)
    return c


@pytest.mark.parametrize("ctor", [
    lambda: SingleGaussianPNN(seed=0, epochs=40),
    lambda: MDNDeepEnsemble(n_members=2, n_comp=2, seed=0, epochs=40),
    lambda: DeepEnsemblePoint(n_members=2, seed=0, epochs=40),
])
def test_nn_baseline_samples_shape_finite_deterministic(tiny, ctor):
    model = ctor().fit(tiny.signals["train"], tiny.params["train"])
    s1 = model.predict_samples(tiny.signals["test"], n=64,
                               rng=np.random.default_rng(0))
    s2 = model.predict_samples(tiny.signals["test"], n=64,
                               rng=np.random.default_rng(0))
    assert s1.shape == (80, 3, 64)
    assert np.all(np.isfinite(s1))
    np.testing.assert_array_equal(s1, s2)  # same rng seed -> identical draws


def test_quantiles_are_ordered(tiny):
    model = DeepEnsemblePoint(n_members=2, seed=0, epochs=40).fit(
        tiny.signals["train"], tiny.params["train"])
    samp = model.predict_samples(tiny.signals["test"], n=128,
                                 rng=np.random.default_rng(0))
    lo, hi = quantiles_from_samples(samp, alpha=0.1)
    assert lo.shape == hi.shape == (80, 3)
    assert np.all(hi >= lo)


def test_mdn_uncertainty_split_nonnegative(tiny):
    model = MDNDeepEnsemble(n_members=3, n_comp=2, seed=0, epochs=40).fit(
        tiny.signals["train"], tiny.params["train"])
    al, ep = model.uncertainty_split(tiny.signals["test"])
    assert al.shape == ep.shape == (80, 3)
    assert np.all(al >= 0) and np.all(ep >= 0)


def test_bayesian_mcmc_shape_and_acceptance(tiny):
    b = tiny.b
    sig = tiny.signals["test"]
    theta, s0, sigma = _nlls_init_and_noise(sig, b)
    bayes = BayesianIVIM_MCMC(seed=0, n_samples=120, burn=200, thin=2)
    samp = bayes.predict_samples_for(sig, b, theta, s0, sigma,
                                     np.random.default_rng(0))
    assert samp.shape == (80, 3, 120)
    assert np.all(np.isfinite(samp))
    assert 0.05 < bayes.accept_rate < 0.95
    # posterior draws must respect the physiological prior bounds
    from gauge.cohort import D_RANGE, DSTAR_RANGE, F_RANGE
    assert samp[:, 0].min() >= D_RANGE[0] - 1e-9
    assert samp[:, 1].max() <= DSTAR_RANGE[1] + 1e-9
    assert samp[:, 2].min() >= F_RANGE[0] - 1e-9
