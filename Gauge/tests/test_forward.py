"""Property tests for the bi-exponential IVIM forward model and Rician noise.

These assert exact mathematical facts about the model, so they double as the
correctness spec: S(b=0)=S0, the two mono-exponential degeneracies, monotonic
decay, the high-b log-slope recovering D, and the Rician noise floor sigma*sqrt(pi/2).
"""
import numpy as np
import pytest

from gauge.forward import ivim_signal, add_rician_noise, DEFAULT_B_VALUES


# Physical-unit reference parameters (mm^2/s), liver-like.
D_REF, DSTAR_REF, F_REF = 1.5e-3, 50e-3, 0.20


def test_signal_at_b0_equals_s0():
    for S0 in (1.0, 0.5, 1234.0):
        for D, Dstar, f in [(1.0e-3, 20e-3, 0.1), (2.5e-3, 90e-3, 0.35)]:
            assert ivim_signal(0.0, D, Dstar, f, S0=S0) == pytest.approx(S0)


def test_f_zero_is_monoexponential():
    b = DEFAULT_B_VALUES
    got = ivim_signal(b, D_REF, DSTAR_REF, 0.0, S0=1.0)
    expected = np.exp(-b * D_REF)
    np.testing.assert_allclose(got, expected, rtol=1e-12, atol=0)


def test_dstar_equals_d_is_monoexponential_for_any_f():
    b = DEFAULT_B_VALUES
    for f in (0.0, 0.2, 0.5, 0.9):
        got = ivim_signal(b, D_REF, D_REF, f, S0=1.0)
        expected = np.exp(-b * D_REF)
        np.testing.assert_allclose(got, expected, rtol=1e-12, atol=0)


def test_signal_strictly_decreasing_in_b():
    b = np.linspace(0, 1000, 200)
    s = ivim_signal(b, D_REF, DSTAR_REF, F_REF, S0=1.0)
    assert np.all(np.diff(s) < 0)


def test_high_b_logslope_recovers_D():
    # With Dstar >> D the perfusion term vanishes for large b, leaving
    # log S ~= log(S0(1-f)) - b*D.  A linear fit on high-b points must recover D.
    b = np.array([300, 400, 500, 600, 700, 800], float)
    S0 = 1.0
    s = ivim_signal(b, D_REF, DSTAR_REF, F_REF, S0=S0)
    slope, intercept = np.polyfit(b, np.log(s), 1)
    assert -slope == pytest.approx(D_REF, rel=1e-3)
    assert np.exp(intercept) == pytest.approx(S0 * (1 - F_REF), rel=1e-3)


def test_rician_noise_is_deterministic_with_seed():
    clean = ivim_signal(DEFAULT_B_VALUES, D_REF, DSTAR_REF, F_REF)
    a = add_rician_noise(clean, snr=25.0, rng=np.random.default_rng(7))
    b = add_rician_noise(clean, snr=25.0, rng=np.random.default_rng(7))
    np.testing.assert_array_equal(a, b)


def test_rician_noise_floor_on_zero_signal():
    # E[|Rician(0, sigma)|] = sigma * sqrt(pi/2).
    sigma = 0.04
    snr = 1.0 / sigma
    rng = np.random.default_rng(0)
    zeros = np.zeros(200_000)
    noisy = add_rician_noise(zeros, snr=snr, rng=rng, S0=1.0)
    assert noisy.mean() == pytest.approx(sigma * np.sqrt(np.pi / 2), rel=2e-2)


def test_high_snr_approximates_clean_signal():
    clean = ivim_signal(DEFAULT_B_VALUES, D_REF, DSTAR_REF, F_REF)
    noisy = add_rician_noise(clean, snr=5000.0, rng=np.random.default_rng(1))
    np.testing.assert_allclose(noisy, clean, atol=2e-3)


def test_noise_std_scales_inversely_with_snr():
    clean = ivim_signal(np.full(50_000, 0.5), D_REF, DSTAR_REF, F_REF)  # mid-range signal
    rng = np.random.default_rng(3)
    std_lo = add_rician_noise(clean, snr=10.0, rng=rng).std()
    std_hi = add_rician_noise(clean, snr=40.0, rng=rng).std()
    # 4x SNR -> ~4x smaller noise std (Rician ~ Gaussian away from the floor).
    assert std_lo / std_hi == pytest.approx(4.0, rel=0.15)
