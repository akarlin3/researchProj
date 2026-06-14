"""Tests for the Gauge 04 robustness building blocks (CP0 / CP1).

These verify the *mechanisms* on toy/analytic cases -- the forward-model
extensions, the weighted (covariate-shift) conformal fix, the Minos-style
deployment monitor, and the CRLB acquisition design -- without running the full
(slow) NLLS pipeline. The end-to-end shift/coverage numbers are produced and
gate-checked by ``python -m gauge.robustness``.
"""
import numpy as np
import pytest

from gauge.forward import (ivim_signal, ivim_signal_triexp, add_gaussian_noise,
                           add_rician_noise, design_crlb_dstar,
                           DEFAULT_B_VALUES)
from gauge.conformal import (weighted_split_conformal, conformal_quantile,
                             empirical_coverage)
from gauge.monitor import DeploymentMonitor
from gauge.robustness import (_importance_weights, crlb_optimal_scheme,
                              CLINICAL_B)


# --------------------------------------------------------------------------- #
# forward-model extensions
# --------------------------------------------------------------------------- #
def test_triexp_reduces_to_biexp_at_g_zero():
    b = DEFAULT_B_VALUES
    bi = ivim_signal(b, 1.5e-3, 40e-3, 0.2)
    tri = ivim_signal_triexp(b, 1.5e-3, 40e-3, 0.2, 200e-3, 0.0)
    assert np.allclose(bi, tri)


def test_triexp_perturbs_signal_when_g_positive():
    b = DEFAULT_B_VALUES
    bi = ivim_signal(b, 1.5e-3, 40e-3, 0.2)
    tri = ivim_signal_triexp(b, 1.5e-3, 40e-3, 0.2, 200e-3, 0.25)
    # the fast third compartment lifts the low-b signal but not b=0.
    assert tri[0] == pytest.approx(bi[0])
    assert np.any(np.abs(tri - bi) > 1e-3)


def test_gaussian_noise_is_unbiased_and_deterministic():
    rng = np.random.default_rng(0)
    clean = np.full(5000, 0.5)
    noisy = add_gaussian_noise(clean, snr=20.0, rng=rng, S0=1.0)
    assert abs(noisy.mean() - 0.5) < 0.01           # no Rician floor / bias
    assert np.std(noisy) == pytest.approx(0.05, rel=0.1)
    # determinism: same seed -> same draw
    a = add_gaussian_noise(clean, 20.0, np.random.default_rng(1))
    b = add_gaussian_noise(clean, 20.0, np.random.default_rng(1))
    assert np.array_equal(a, b)


def test_rician_floor_lifts_low_signal_above_gaussian():
    # At low signal the Rician mean exceeds the (unbiased) Gaussian mean.
    rng_r = np.random.default_rng(2)
    rng_g = np.random.default_rng(2)
    clean = np.full(20000, 0.05)
    r = add_rician_noise(clean, 10.0, rng_r).mean()
    g = add_gaussian_noise(clean, 10.0, rng_g).mean()
    assert r > g + 0.01


# --------------------------------------------------------------------------- #
# CRLB acquisition design
# --------------------------------------------------------------------------- #
def test_design_crlb_dstar_prefers_low_b_dense_scheme():
    rng = np.random.default_rng(0)
    n = 200
    D = rng.uniform(0.8e-3, 2.0e-3, n)
    Dstar = rng.uniform(10e-3, 100e-3, n)
    f = rng.uniform(0.1, 0.3, n)
    snr = np.full(n, 30.0)
    lowb = np.array([0, 5, 10, 20, 40, 60, 80, 100, 200, 800], float)
    highb = np.array([0, 100, 200, 300, 400, 500, 600, 700, 750, 800], float)
    _, hi_low, _ = design_crlb_dstar(lowb, D, Dstar, f, snr)
    _, hi_high, _ = design_crlb_dstar(highb, D, Dstar, f, snr)
    # D* (fast compartment) needs dense low-b sampling to be identifiable.
    assert hi_low < hi_high


def test_crlb_optimal_scheme_includes_b0_and_correct_count():
    b = crlb_optimal_scheme(n_b=11, seed=0)
    assert len(b) == 11
    assert b[0] == 0.0
    assert np.all(np.diff(b) > 0)                   # sorted, unique
    # a sensible D* design concentrates samples at low b.
    assert (b < 100).sum() >= 6


# --------------------------------------------------------------------------- #
# weighted (covariate-shift) conformal
# --------------------------------------------------------------------------- #
def test_weighted_split_reduces_to_unweighted_with_uniform_weights():
    rng = np.random.default_rng(0)
    scores = rng.exponential(size=500)
    w_cal = np.ones(500)
    w_test = np.ones(10)
    q_plain = conformal_quantile(scores, 0.1)
    q_w = weighted_split_conformal(scores, w_cal, w_test, 0.1)
    assert np.allclose(q_w, q_plain)


def test_weighted_conformal_restores_coverage_under_covariate_shift():
    # Heteroscedastic truth y|x ~ N(0, (1+3x)^2); calibration x~U(0,1), test x
    # shifted toward large x (bigger noise). Plain split under-covers; weighting
    # by the known density ratio restores ~nominal coverage.
    rng = np.random.default_rng(0)
    alpha = 0.1
    n = 6000
    xc = rng.uniform(0, 1, n)
    sc = np.abs(rng.normal(0, 1 + 3 * xc))           # |residual| scores
    xt = rng.beta(5, 2, n)                            # shifted toward 1
    yt = rng.normal(0, 1 + 3 * xt)
    # plain split radius
    q = conformal_quantile(sc, alpha)
    cov_plain = float(np.mean(np.abs(yt) <= q))
    # weighted: density ratio U(0,1)->Beta(5,2): w(x) = beta.pdf(x;5,2)
    from scipy.stats import beta as beta_dist
    w_cal = beta_dist.pdf(xc, 5, 2)
    w_test = beta_dist.pdf(xt, 5, 2)
    qj = weighted_split_conformal(sc, w_cal, w_test, alpha)
    cov_w = float(np.mean(np.abs(yt) <= qj))
    assert cov_plain < 0.87                          # plain under-covers
    assert abs(cov_w - (1 - alpha)) < 0.03           # weighting restores it


# --------------------------------------------------------------------------- #
# Minos-style deployment monitor (observable vs hidden)
# --------------------------------------------------------------------------- #
def test_monitor_silent_on_indistribution_fires_on_shift():
    rng = np.random.default_rng(0)
    cal_feat = rng.normal(0, 1, (1500, 6))
    cal_resid = rng.exponential(1.0, 1500)
    mon = DeploymentMonitor(seed=0).fit(cal_feat, cal_resid)
    # in-distribution test: should not fire; AUC ~ 0.5 (hidden-blind regime)
    id_feat = rng.normal(0, 1, (1500, 6))
    id_resid = rng.exponential(1.0, 1500)
    r_id = mon.evaluate(id_feat, id_resid)
    assert not r_id["fires"]
    assert r_id["auc"] < 0.6
    # observable shift (mean + scale move): should fire with high AUC
    sh_feat = rng.normal(1.5, 1.4, (1500, 6))
    sh_resid = rng.exponential(2.5, 1500)
    r_sh = mon.evaluate(sh_feat, sh_resid)
    assert r_sh["fires"]
    assert r_sh["auc"] > 0.8


def test_importance_weights_separate_shifted_populations():
    rng = np.random.default_rng(0)
    cal = rng.normal(0, 1, (1200, 4))
    test = rng.normal(2.0, 1, (1200, 4))             # clearly shifted
    w_cal, w_test, _ = _importance_weights(cal, test, seed=0)
    # test-like calibration points (high) get up-weighted vs the bulk.
    assert w_cal.shape == (1200,)
    assert w_test.shape == (1200,)
    assert np.median(w_test) > np.median(w_cal)


def test_clinical_scheme_is_a_valid_bvalue_set():
    assert CLINICAL_B[0] == 0.0
    assert np.all(np.diff(CLINICAL_B) > 0)
    assert CLINICAL_B.max() == 800.0
