"""Correctness tests for split-conformal and CQR.

These run on toy distributions where the coverage answer is known analytically,
so they verify the *guarantee* itself (finite-sample marginal coverage under
exchangeability) independent of any IVIM specifics. The decisive test is that a
deliberately biased base predictor STILL achieves nominal coverage -- that is
exactly what distribution-free conformal prediction promises.
"""
import numpy as np
import pytest

from gauge.conformal import (
    conformal_quantile,
    split_conformal,
    cqr,
    empirical_coverage,
    interval_width,
)


def test_conformal_quantile_is_kth_order_statistic():
    # n=9, alpha=0.1 -> k = ceil(10*0.9) = 9 -> 9th smallest = max.
    scores = np.array([3, 1, 4, 1, 5, 9, 2, 6, 5], float)
    assert conformal_quantile(scores, 0.1) == 9.0
    # alpha=0.5 -> k = ceil(10*0.5) = 5 -> 5th smallest.
    assert conformal_quantile(scores, 0.5) == np.sort(scores)[4]


def test_conformal_quantile_infinite_when_calibration_too_small():
    scores = np.arange(5, dtype=float)
    # k = ceil(6 * 0.99) = 6 > 5 -> no finite quantile -> inf.
    assert conformal_quantile(scores, 0.01) == np.inf


@pytest.mark.parametrize("alpha", [0.05, 0.10, 0.20, 0.30])
def test_split_conformal_tracks_nominal_on_gaussian(alpha):
    rng = np.random.default_rng(0)
    mu, sigma = 5.0, 2.0
    n_cal, n_test = 4000, 20000
    cal_true = rng.normal(mu, sigma, n_cal)
    test_true = rng.normal(mu, sigma, n_test)
    cal_pred = np.full(n_cal, mu)          # correct point predictor
    test_pred = np.full(n_test, mu)
    lo, hi, _ = split_conformal(cal_pred, cal_true, test_pred, alpha)
    cov = empirical_coverage(lo, hi, test_true)
    assert cov == pytest.approx(1 - alpha, abs=0.02)


@pytest.mark.parametrize("alpha", [0.05, 0.10, 0.20])
def test_split_conformal_covers_despite_biased_base(alpha):
    # The distribution-free guarantee: a biased point predictor still yields
    # valid marginal coverage (the calibrated radius absorbs the bias).
    rng = np.random.default_rng(1)
    mu, sigma, bias = 0.0, 1.0, 4.0
    n_cal, n_test = 4000, 20000
    cal_true = rng.normal(mu, sigma, n_cal)
    test_true = rng.normal(mu, sigma, n_test)
    cal_pred = np.full(n_cal, mu + bias)
    test_pred = np.full(n_test, mu + bias)
    lo, hi, _ = split_conformal(cal_pred, cal_true, test_pred, alpha)
    cov = empirical_coverage(lo, hi, test_true)
    assert cov == pytest.approx(1 - alpha, abs=0.02)


def test_split_conformal_marginal_guarantee_in_expectation():
    # Average realized coverage over many calibration draws must be >= 1-alpha
    # (the finite-sample lower bound holds in expectation over calibration sets).
    alpha = 0.10
    mu, sigma = 1.0, 0.5
    covs = []
    for s in range(40):
        rng = np.random.default_rng(100 + s)
        cal_true = rng.normal(mu, sigma, 1000)
        test_true = rng.normal(mu, sigma, 4000)
        lo, hi, _ = split_conformal(np.full(1000, mu), cal_true,
                                    np.full(4000, mu), alpha)
        covs.append(empirical_coverage(lo, hi, test_true))
    assert np.mean(covs) >= 1 - alpha - 0.005


@pytest.mark.parametrize("alpha", [0.05, 0.10, 0.20])
def test_cqr_covers_with_miscalibrated_heteroscedastic_band(alpha):
    # Heteroscedastic noise with an x-dependent but wrongly-scaled base band.
    # CQR's additive correction must restore marginal coverage.
    rng = np.random.default_rng(2)
    n_cal, n_test = 6000, 20000
    x_cal = rng.uniform(0, 1, n_cal)
    x_test = rng.uniform(0, 1, n_test)
    sig = lambda x: 0.5 + 2.0 * x
    cal_true = rng.normal(0.0, sig(x_cal))
    test_true = rng.normal(0.0, sig(x_test))
    # Base "quantile" band: shape-aware but only a 1-sigma guess (miscalibrated).
    band = lambda x: sig(x)
    lo, hi, _ = cqr(-band(x_cal), band(x_cal), cal_true,
                    -band(x_test), band(x_test), alpha)
    cov = empirical_coverage(lo, hi, test_true)
    assert cov == pytest.approx(1 - alpha, abs=0.02)


def test_cqr_intervals_widen_where_noise_grows():
    # Sanity: CQR keeps the heteroscedastic shape (wider band at large x).
    rng = np.random.default_rng(3)
    x_cal = rng.uniform(0, 1, 6000)
    x_test = np.array([0.05, 0.95])
    sig = lambda x: 0.5 + 2.0 * x
    cal_true = rng.normal(0.0, sig(x_cal))
    lo, hi, _ = cqr(-sig(x_cal), sig(x_cal), cal_true,
                    -sig(x_test), sig(x_test), 0.1)
    w = interval_width(lo, hi)
    assert w[1] > w[0]


def test_empirical_coverage_and_width_basics():
    lo = np.array([0.0, 0.0, 0.0])
    hi = np.array([1.0, 1.0, 1.0])
    true = np.array([0.5, 2.0, -1.0])
    assert empirical_coverage(lo, hi, true) == pytest.approx(1 / 3)
    np.testing.assert_array_equal(interval_width(lo, hi), np.ones(3))
