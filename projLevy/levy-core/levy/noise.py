r"""Rician measurement noise for magnitude diffusion-MRI signals, and the per-sample
Rician Fisher-information factor that the CRLB layer consumes.

SNR convention (fixed once, used everywhere)
--------------------------------------------
    sigma = S0 / SNR
i.e. SNR is the signal-to-noise ratio at b = 0 (the non-diffusion-weighted amplitude
over the Gaussian noise std of each magnitude channel). This is the standard DWI
convention.

Rician magnitude model
----------------------
A magnitude measurement is M = sqrt((nu + n1)^2 + n2^2) with n1, n2 ~ N(0, sigma^2)
independent, where nu = S(b; theta) >= 0 is the noise-free signal. Its density is
    p(M | nu, sigma) = (M/sigma^2) exp(-(M^2 + nu^2)/(2 sigma^2)) I0(M nu / sigma^2),  M >= 0.

Per-sample Fisher information about nu
--------------------------------------
The score is  d/dnu log p = (1/sigma^2)(M r(z) - nu),  z = M nu / sigma^2,  r(x)=I1(x)/I0(x).
The Fisher information I_R(nu, sigma) = E[ score^2 ] depends on (nu, sigma) only through the
local SNR  a = nu/sigma, via  I_R(nu, sigma) = f(a) / sigma^2  where

    f(a) = \int_0^inf (m r(m a) - a)^2 * m * exp(-(m-a)^2/2) * Ive(0, m a) dm        (m = M/sigma)

with Ive(0, .) the exponentially-scaled modified Bessel I0 (numerically stable). As a -> inf,
f(a) -> 1 (Rician -> Gaussian, recovering the 1/sigma^2 Gaussian information); as a -> 0 the
information about nu collapses. ``rician_info_factor`` returns f(a) (cached interpolant), and
this is what turns the high-SNR Gaussian CRLB into the honest Rician CRLB at finite SNR.
"""
from __future__ import annotations

import numpy as np
from scipy import integrate, special


def sigma_from_snr(S0: float, snr: float) -> float:
    """sigma = S0 / SNR (b=0 magnitude SNR convention)."""
    return float(S0) / float(snr)


def rician_sample(nu, sigma, rng):
    """Draw Rician magnitude samples for noise-free signal(s) ``nu`` and noise ``sigma``.

    ``rng`` is an ``np.random.Generator`` (explicit-Generator discipline; no bare np.random).
    """
    nu = np.asarray(nu, dtype=float)
    n1 = rng.normal(0.0, sigma, size=nu.shape)
    n2 = rng.normal(0.0, sigma, size=nu.shape)
    return np.sqrt((nu + n1) ** 2 + n2 ** 2)


def rician_logpdf(M, nu, sigma):
    """log p(M | nu, sigma), numerically stable via the exponentially-scaled Bessel Ive.

    log I0(z) = log Ive(0, z) + z, so the exp argument becomes -(M-nu)^2/(2 sigma^2),
    which never overflows.
    """
    M = np.asarray(M, dtype=float)
    nu = np.asarray(nu, dtype=float)
    s2 = sigma * sigma
    z = M * nu / s2
    # log I0(z) = log Ive(0,z) + z; the +z combines with -(M^2+nu^2)/2s2 into -(M-nu)^2/2s2.
    # Add ONLY log Ive(0,z) here -- the z is already folded into the (M-nu)^2 term.
    out = np.log(np.maximum(M, 1e-300)) - np.log(s2) - (M - nu) ** 2 / (2.0 * s2) + np.log(special.ive(0, z))
    return out


def _f_of_a(a: float) -> float:
    """f(a) = sigma^2 * I_R(nu, sigma), the dimensionless Rician info factor at local SNR a."""
    if a <= 0.0:
        # Rayleigh limit: information about the (zero) signal mean -> small but finite;
        # the integral below handles it, but guard a==0 to avoid log/ratio issues.
        a = 1e-6

    def integrand(m):
        za = m * a
        # r(za) = I1/I0 = Ive(1,za)/Ive(0,za); stable ratio.
        ive0 = special.ive(0, za)
        ive1 = special.ive(1, za)
        r = ive1 / ive0
        # weight = m * exp(-(m-a)^2/2) * Ive(0, m a)  [ = m * exp(-(m^2+a^2)/2) * I0(ma) ]
        w = m * np.exp(-0.5 * (m - a) ** 2) * ive0
        return (m * r - a) ** 2 * w

    # Integrand peaks near m ~ a and decays like a Gaussian; integrate generously.
    hi = a + 20.0
    val, _ = integrate.quad(integrand, 0.0, hi, limit=200)
    return float(val)


# Cached interpolant of f(a) over a log-spaced grid of local SNR a, so the wall sweep is fast.
_A_GRID = np.concatenate([[0.0], np.geomspace(1e-3, 200.0, 400)])
_F_GRID = np.array([_f_of_a(a) for a in _A_GRID])


def rician_info_factor(a):
    """f(a): dimensionless Rician Fisher-info factor at local SNR a = nu/sigma.

    Interpolated from a precomputed quadrature grid; f(a) -> 1 as a -> inf (Gaussian limit).
    Vectorized over ``a``.
    """
    a = np.asarray(a, dtype=float)
    a_clip = np.clip(a, 0.0, _A_GRID[-1])
    f = np.interp(a_clip, _A_GRID, _F_GRID)
    # For a beyond the grid, the Gaussian limit f->1 holds to <1e-3; clamp to 1.0.
    f = np.where(a > _A_GRID[-1], 1.0, f)
    return f
