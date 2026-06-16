"""Tests for the alternative-forward-model cohort (circularity + envelope).

Fast tests only: the forward-model physics, the continuity-to-bi-exponential
limit, the effective-D* surrogate, the dispersion-Jacobian CRLB, and determinism
of the cohort draw. The full Arm-1/Arm-2 pipeline is exercised by the seeded
report + multiseed sweep, not in unit tests.
"""
import numpy as np
import pytest

from gauge.forward import (ivim_signal, ivim_signal_dispersion,
                           ivim_signal_lognormal_dispersion,
                           ivim_signal_stretched, dispersion_dstar_eff,
                           dispersion_crlb_mu_batch, crlb_dstar_batch,
                           DEFAULT_B_VALUES)
from gauge.robustness import _simulate, _draw_params, _draw_snr
from gauge.cohort import DEFAULT_SNR_GRID

B = DEFAULT_B_VALUES


def test_dispersion_biexp_limit_exact():
    """Gamma-dispersion at k=inf (CV=0) reproduces the bi-exponential exactly."""
    rng = np.random.default_rng(0)
    D = rng.uniform(0.5e-3, 3e-3, 200)
    mu = rng.uniform(10e-3, 100e-3, 200)
    f = rng.uniform(0.05, 0.40, 200)
    disp = ivim_signal_dispersion(B[None, :], D[:, None], mu[:, None], np.inf,
                                  f[:, None])
    biexp = ivim_signal(B[None, :], D[:, None], mu[:, None], f[:, None])
    assert np.max(np.abs(disp - biexp)) == 0.0


def test_stretched_biexp_limit_exact():
    """Stretched-exponential perfusion at beta=1 reproduces the bi-exponential."""
    rng = np.random.default_rng(1)
    D = rng.uniform(0.5e-3, 3e-3, 200)
    Dstar = rng.uniform(10e-3, 100e-3, 200)
    f = rng.uniform(0.05, 0.40, 200)
    st = ivim_signal_stretched(B[None, :], D[:, None], Dstar[:, None],
                               f[:, None], 1.0)
    biexp = ivim_signal(B[None, :], D[:, None], Dstar[:, None], f[:, None])
    assert np.max(np.abs(st - biexp)) == 0.0


def test_dispersion_finite_k_departs():
    """Finite k (real dispersion) is genuinely NOT bi-exponential."""
    biexp = ivim_signal(B, 1.5e-3, 40e-3, 0.2)
    disp = ivim_signal_dispersion(B, 1.5e-3, 40e-3, 4.0, 0.2)
    assert np.max(np.abs(disp - biexp)) > 1e-3


def test_dstar_eff_is_gamma_mean():
    """D*eff (surrogate A) is the gamma mean mu, independent of shape k."""
    mu = np.array([20e-3, 50e-3, 90e-3])
    assert np.allclose(dispersion_dstar_eff(mu), mu)
    # numeric check that mu IS the initial perfusion log-slope, for finite k.
    b = np.array([0.0, 0.5, 1.0])           # tiny b for a clean finite-difference
    for k in (2.0, 8.0, np.inf):
        perf = ivim_signal_dispersion(b, 0.0, 0.05, k, 1.0)   # f=1, D=0 -> perfusion only
        slope0 = -(np.log(perf[1]) - np.log(perf[0])) / (b[1] - b[0])
        assert slope0 == pytest.approx(0.05, rel=2e-2)


def test_dispersion_crlb_finite_and_grows():
    """Dispersion-Jacobian CRLB(mu) is finite, positive, and grows with mu."""
    N = 40
    rng = np.random.default_rng(2)
    D = np.full(N, 1.5e-3)
    mu = np.linspace(12e-3, 98e-3, N)
    f = np.full(N, 0.2)
    k = np.full(N, 4.0)
    snr = np.full(N, 30.0)
    sd = dispersion_crlb_mu_batch(B, D, mu, k, f, snr, fix_k=True)
    assert np.all(np.isfinite(sd)) and np.all(sd > 0)
    assert np.median(sd[-N // 3:]) > np.median(sd[:N // 3])   # balloons at high mu


def test_simulate_dispersion_cv0_matches_biexp():
    """The _simulate chokepoint at model='dispersion', disp_cv=0 == biexp clean."""
    rng = np.random.default_rng(3)
    params = _draw_params(300, rng)
    snr = _draw_snr(300, rng, DEFAULT_SNR_GRID)
    r1 = np.random.default_rng(7)
    r2 = np.random.default_rng(7)
    sig_disp = _simulate(B, params, snr, r1, model="dispersion", disp_cv=0.0)
    sig_biexp = _simulate(B, params, snr, r2, model="biexp")
    # same RNG state + exact clean-signal limit -> byte-identical noisy signals
    assert np.array_equal(sig_disp, sig_biexp)


def test_dispersion_pool_deterministic():
    """The Arm-1 dispersion cohort is a pure function of the seed."""
    from gauge.altmodel import _dispersion_pool
    p1, s1, g1 = _dispersion_pool(20260613, 0.5, B)
    p2, s2, g2 = _dispersion_pool(20260613, 0.5, B)
    assert np.array_equal(p1, p2) and np.array_equal(g1, g2) and np.array_equal(s1, s2)


# --------------------------------------------------------------------------- #
# Second dispersion kernel: log-normal (independent shape, same D*eff = mean).
# --------------------------------------------------------------------------- #
def test_lognormal_biexp_limit_exact():
    """Log-normal dispersion at CV=0 reproduces the bi-exponential exactly."""
    rng = np.random.default_rng(10)
    D = rng.uniform(0.5e-3, 3e-3, 200)
    mu = rng.uniform(10e-3, 100e-3, 200)
    f = rng.uniform(0.05, 0.40, 200)
    ln = ivim_signal_lognormal_dispersion(B[None, :], D[:, None], mu[:, None],
                                          0.0, f[:, None])
    biexp = ivim_signal(B[None, :], D[:, None], mu[:, None], f[:, None])
    assert np.max(np.abs(ln - biexp)) == 0.0


def test_lognormal_finite_cv_departs_and_differs_from_gamma():
    """Finite CV is genuinely non-bi-exp AND a distinct shape from the gamma kernel."""
    biexp = ivim_signal(B, 1.5e-3, 40e-3, 0.2)
    ln = ivim_signal_lognormal_dispersion(B, 1.5e-3, 40e-3, 0.5, 0.2)
    gam = ivim_signal_dispersion(B, 1.5e-3, 40e-3, 1.0 / 0.25, 0.2)   # CV=0.5 -> k=4
    assert np.max(np.abs(ln - biexp)) > 1e-3                          # off bi-exp
    assert np.max(np.abs(ln - gam)) > 1e-4                            # distinct from gamma


def test_lognormal_dstar_eff_is_kernel_mean():
    """The log-normal kernel's initial perfusion log-slope is the mean mu (= D*eff)."""
    b = np.array([0.0, 1e-3])                # b -> 0 for the analytic initial slope
    for cv in (0.3, 0.5, 0.8):
        perf = ivim_signal_lognormal_dispersion(b, 0.0, 0.05, cv, 1.0)  # f=1, D=0
        slope0 = -(perf[1] - perf[0]) / (b[1] - b[0])
        assert slope0 == pytest.approx(0.05, rel=2e-3)


def test_lognormal_pool_deterministic():
    """The log-normal Arm-1 cohort is a pure function of the seed."""
    from gauge.altmodel import _lognormal_pool
    p1, s1, g1 = _lognormal_pool(20260613, 0.5, B)
    p2, s2, g2 = _lognormal_pool(20260613, 0.5, B)
    assert np.array_equal(p1, p2) and np.array_equal(g1, g2) and np.array_equal(s1, s2)
