"""Rician noise model + the per-sample Rician Fisher-information factor."""
import numpy as np
from scipy import integrate

from levy import noise, seeding


def test_info_factor_gaussian_limit():
    # f(a) -> 1 as local SNR a -> infinity (Rician -> Gaussian)
    assert noise.rician_info_factor(100.0) > 0.99
    assert noise.rician_info_factor(50.0) > 0.98
    # and is < 1 at low SNR (information about the mean is degraded)
    assert noise.rician_info_factor(1.0) < 0.95
    assert noise.rician_info_factor(0.3) < noise.rician_info_factor(3.0)


def test_info_factor_monotone_increasing():
    a = np.array([0.1, 0.5, 1.0, 2.0, 5.0, 20.0, 100.0])
    f = noise.rician_info_factor(a)
    assert np.all(np.diff(f) >= -1e-9)


def test_logpdf_integrates_to_one():
    nu, sigma = 3.0, 1.0
    val, _ = integrate.quad(lambda M: np.exp(noise.rician_logpdf(M, nu, sigma)), 0, 40, limit=200)
    assert abs(val - 1.0) < 1e-4


def test_rician_sample_high_snr_mean():
    rng = seeding.make_rng(1)
    nu, sigma = 10.0, 1.0  # SNR 10 -> Rician mean ~ nu
    M = noise.rician_sample(np.full(20000, nu), sigma, rng)
    assert abs(M.mean() - nu) < 0.05


def test_info_factor_vs_direct_quadrature():
    # cross-check the cached interpolant against a direct E[score^2] estimate at a=4
    a = 4.0
    f_interp = float(noise.rician_info_factor(a))
    # direct: I_R = E[score^2] with sigma=1, nu=a; score=(M r(z) - nu), z=M nu
    from scipy import special
    def integrand(M):
        z = M * a
        r = special.ive(1, z) / special.ive(0, z)
        score2 = (M * r - a) ** 2
        pdf = M * np.exp(-(M - a) ** 2 / 2) * special.ive(0, M * a)
        return score2 * pdf
    direct, _ = integrate.quad(integrand, 0, a + 20, limit=200)
    assert abs(f_interp - direct) < 5e-3
