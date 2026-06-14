"""Correctness tests for the Gauge 03 conditional-attack primitives.

Two families of guarantee are checked on toy data where the answer is known:

* the feature-conditional conformal methods (localized / conditional conformal)
  must (a) keep ~nominal MARGINAL coverage and (b) actually improve CONDITIONAL
  coverage on feature-dependent heteroscedastic noise, where a single global
  radius is known to fail; and
* the IVIM Fisher/CRLB machinery must match a finite-difference Jacobian and
  reproduce the qualitative identifiability facts (D* gets harder as D* grows,
  easier as SNR grows).
"""
import numpy as np
import pytest

from gauge.conformal import (
    conformal_quantile,
    split_conformal,
    empirical_coverage,
    weighted_conformal_quantile,
    localized_conformal,
    conditional_conformal,
)
from gauge.forward import (
    ivim_jacobian,
    ivim_signal,
    crlb,
    crlb_dstar_batch,
    DEFAULT_B_VALUES,
)


# --------------------------------------------------------------------------- #
# weighted conformal quantile
# --------------------------------------------------------------------------- #
def test_weighted_quantile_reduces_to_unweighted():
    # Uniform weights + matching self-weight == the standard conformal quantile.
    rng = np.random.default_rng(0)
    scores = rng.gamma(2.0, 1.0, size=200)
    for alpha in (0.05, 0.1, 0.2, 0.5):
        w = np.ones_like(scores)
        got = weighted_conformal_quantile(scores, w, alpha, self_weight=1.0)
        assert got == pytest.approx(conformal_quantile(scores, alpha))


def test_weighted_quantile_infinite_when_unreachable():
    scores = np.arange(5, dtype=float)
    # tiny calibration mass relative to a huge self-weight -> 1-alpha unreachable
    w = np.ones_like(scores)
    assert weighted_conformal_quantile(scores, w, 0.01, self_weight=100.0) == np.inf


# --------------------------------------------------------------------------- #
# localized conformal (Guan 2023)
# --------------------------------------------------------------------------- #
def test_localized_conformal_marginal_coverage():
    rng = np.random.default_rng(1)
    n_cal, n_test = 4000, 8000
    x_cal = rng.uniform(0, 1, (n_cal, 1))
    x_test = rng.uniform(0, 1, (n_test, 1))
    sig = lambda x: 0.3 + 2.0 * x[:, 0]
    y_cal = rng.normal(0.0, sig(x_cal))
    y_test = rng.normal(0.0, sig(x_test))
    scores = np.abs(y_cal)                              # predictor = 0
    t = localized_conformal(scores, x_cal, x_test, 0.1, bandwidth=0.25)
    cov = empirical_coverage(-t, t, y_test)
    assert cov == pytest.approx(0.9, abs=0.03)


def test_localized_conformal_beats_global_conditionally():
    # On feature-dependent noise, a global radius under-covers at high-sigma x
    # and over-covers at low-sigma x; LCP should flatten the per-bin coverage.
    rng = np.random.default_rng(2)
    n_cal, n_test = 5000, 10000
    x_cal = rng.uniform(0, 1, (n_cal, 1))
    x_test = rng.uniform(0, 1, (n_test, 1))
    sig = lambda x: 0.2 + 3.0 * x[:, 0]
    y_cal = rng.normal(0.0, sig(x_cal))
    y_test = rng.normal(0.0, sig(x_test))
    scores = np.abs(y_cal)

    lo_g, hi_g, _ = split_conformal(np.zeros(n_cal), y_cal,
                                    np.zeros(n_test), 0.1)
    t_l = localized_conformal(scores, x_cal, x_test, 0.1, bandwidth=0.2)

    edges = np.quantile(x_test[:, 0], [1 / 3, 2 / 3])
    reg = np.digitize(x_test[:, 0], edges)

    def worst_gap(lo, hi):
        return max(abs(0.9 - empirical_coverage(lo[reg == r], hi[reg == r],
                                                y_test[reg == r]))
                   for r in range(3))

    assert worst_gap(-t_l, t_l) < worst_gap(lo_g, hi_g)
    assert worst_gap(-t_l, t_l) < 0.05


# --------------------------------------------------------------------------- #
# conditional conformal (Gibbs, Cherian & Candes 2023)
# --------------------------------------------------------------------------- #
def test_conditional_conformal_marginal_and_conditional():
    rng = np.random.default_rng(3)
    n_cal, n_test = 6000, 12000
    x_cal = rng.uniform(0, 1, (n_cal, 1))
    x_test = rng.uniform(0, 1, (n_test, 1))
    sig = lambda x: 0.2 + 3.0 * x[:, 0]
    y_cal = rng.normal(0.0, sig(x_cal))
    y_test = rng.normal(0.0, sig(x_test))
    scores = np.abs(y_cal)
    # quantile-regress the score on the feature; threshold tracks local scale
    t = np.clip(conditional_conformal(scores, x_cal, x_test, 0.1), 0.0, None)
    cov = empirical_coverage(-t, t, y_test)
    assert cov == pytest.approx(0.9, abs=0.03)

    edges = np.quantile(x_test[:, 0], [1 / 3, 2 / 3])
    reg = np.digitize(x_test[:, 0], edges)
    worst = max(abs(0.9 - empirical_coverage(-t[reg == r], t[reg == r],
                                             y_test[reg == r]))
                for r in range(3))
    assert worst < 0.05


# --------------------------------------------------------------------------- #
# IVIM Fisher information / CRLB
# --------------------------------------------------------------------------- #
def test_jacobian_matches_finite_difference():
    b = DEFAULT_B_VALUES
    D, Dstar, f, S0 = 1.5e-3, 30e-3, 0.2, 1.0
    J = ivim_jacobian(b, D, Dstar, f, S0)
    base = ivim_signal(b, D, Dstar, f, S0=S0)
    steps = {0: 1e-7, 1: 1e-7, 2: 1e-5, 3: 1e-5}
    args = [D, Dstar, f, S0]
    for k, h in steps.items():
        up = list(args)
        up[k] += h
        fd = (ivim_signal(b, up[0], up[1], up[2], S0=up[3]) - base) / h
        np.testing.assert_allclose(J[:, k], fd, rtol=1e-3, atol=1e-6)


def test_crlb_dstar_grows_with_dstar_and_shrinks_with_snr():
    b = DEFAULT_B_VALUES
    # under-identification: larger true D* -> larger CRLB(D*)
    lo = crlb(b, 1.5e-3, 20e-3, 0.2, snr=30)["Dstar"]
    hi = crlb(b, 1.5e-3, 95e-3, 0.2, snr=30)["Dstar"]
    assert hi > lo
    # more information: higher SNR -> smaller CRLB(D*)
    noisy = crlb(b, 1.5e-3, 90e-3, 0.2, snr=10)["Dstar"]
    clean = crlb(b, 1.5e-3, 90e-3, 0.2, snr=100)["Dstar"]
    assert clean < noisy


def test_crlb_batch_matches_scalar():
    b = DEFAULT_B_VALUES
    D = np.array([1.0e-3, 1.5e-3, 2.0e-3])
    Dstar = np.array([20e-3, 50e-3, 90e-3])
    f = np.array([0.1, 0.2, 0.3])
    snr = np.array([10.0, 30.0, 100.0])
    batch = crlb_dstar_batch(b, D, Dstar, f, snr)
    for i in range(3):
        scalar = crlb(b, D[i], Dstar[i], f[i], snr=snr[i])["Dstar"]
        assert batch[i] == pytest.approx(scalar, rel=1e-6)
