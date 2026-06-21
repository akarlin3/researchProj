"""Clean-room IVIM bi-exponential forward model.

Standard intravoxel-incoherent-motion signal (Le Bihan et al., Radiology 1988),
reimplemented from the physics -- no code from Fashion or Caliper:

    S(b)/S0 = (1 - f) * exp(-b * D)  +  f * exp(-b * (D + Dstar))

Units throughout: ``D``, ``Dstar`` in mm^2/s; ``b`` in s/mm^2; ``f`` dimensionless.
Because ``Dstar >> D`` and the fast term decays within the first few b-values, the
pseudo-diffusion compartment is weakly identified at clinical-sparse sampling -- the
mechanism behind the boundary-railing (T1) and the D* under-coverage (T3) Gnomon
reproduces.

Numerics: estimators fit in **scaled** parameters ``(S0, D3, f, Ds3)`` with
``D = D3 * 1e-3`` and ``Dstar = Ds3 * 1e-3`` so every fit variable is O(1) and the
Jacobian is well-conditioned. :func:`ivim_design` and :func:`ivim_jac_scaled` are
the scaled-space model and analytic Jacobian used by NLLS, the Laplace/CRLB
covariance, and the MCMC likelihood.
"""
from __future__ import annotations

import numpy as np

# Canonical parameter order for the *physical* model and the *scaled* fit vector.
PARAM_NAMES = ("D", "Dstar", "f")
FIT_NAMES = ("S0", "D3", "f", "Ds3")  # scaled: D = D3*1e-3, Dstar = Ds3*1e-3
_SCALE = 1e-3  # mm^2/s per scaled unit for D, Dstar


def ivim(b, D, Dstar, f, s0=1.0):
    """IVIM bi-exponential signal in physical units.

    Broadcasts: ``b`` shape ``(nb,)``; ``D, Dstar, f, s0`` scalars or shape ``(n,)``.
    Returns ``(nb,)`` for scalar params or ``(n, nb)`` for vector params.
    """
    b = np.asarray(b, dtype=float)
    D = np.asarray(D, dtype=float)
    Dstar = np.asarray(Dstar, dtype=float)
    f = np.asarray(f, dtype=float)
    s0 = np.asarray(s0, dtype=float)
    if max(D.ndim, Dstar.ndim, f.ndim, s0.ndim) == 0:  # all scalar -> (nb,)
        return s0 * ((1.0 - f) * np.exp(-b * D) + f * np.exp(-b * (D + Dstar)))
    # at least one vector param -> (n, nb); reshape vectors to columns, keep scalars
    col = lambda x: x[:, None] if x.ndim == 1 else x
    bb = b[None, :]
    D, Dstar, f, s0 = col(D), col(Dstar), col(f), col(s0)
    return s0 * ((1.0 - f) * np.exp(-bb * D) + f * np.exp(-bb * (D + Dstar)))


def to_fit(D, Dstar, f, s0=1.0):
    """Physical (D, Dstar, f, S0) -> scaled fit vector (S0, D3, f, Ds3)."""
    return np.array([s0, D / _SCALE, f, Dstar / _SCALE], dtype=float)


def from_fit(p):
    """Scaled fit vector (S0, D3, f, Ds3) -> physical dict (broadcasts over rows)."""
    p = np.asarray(p, dtype=float)
    return {"S0": p[..., 0], "D": p[..., 1] * _SCALE, "f": p[..., 2],
            "Dstar": p[..., 3] * _SCALE}


def ivim_design(p, b):
    """Model signal for scaled fit vector ``p=(S0, D3, f, Ds3)`` at b-values ``b``.

    ``p`` shape ``(4,)`` -> ``(nb,)``; ``p`` shape ``(n, 4)`` -> ``(n, nb)``.
    """
    p = np.asarray(p, dtype=float)
    b = np.asarray(b, dtype=float)
    if p.ndim == 1:
        s0, D3, f, Ds3 = p
        return ivim(b, D3 * _SCALE, Ds3 * _SCALE, f, s0)
    s0, D3, f, Ds3 = p[:, 0], p[:, 1], p[:, 2], p[:, 3]
    return ivim(b, D3 * _SCALE, Ds3 * _SCALE, f, s0)


def ivim_jac_scaled(p, b):
    """Analytic Jacobian d S / d (S0, D3, f, Ds3) at scaled ``p``, shape ``(nb, 4)``.

    With ``D = D3*1e-3``, ``Dstar = Ds3*1e-3`` and chain rule for the 1e-3 scaling:
        dS/dS0  = (1-f) e^{-bD} + f e^{-b(D+Dstar)}              ( = S/S0 )
        dS/dD   = -b * S            -> dS/dD3  = -b * S * 1e-3
        dS/df   = S0 ( e^{-b(D+Dstar)} - e^{-bD} )
        dS/dDstar = -S0 f b e^{-b(D+Dstar)} -> dS/dDs3 = (that) * 1e-3
    """
    p = np.asarray(p, dtype=float)
    b = np.asarray(b, dtype=float)
    s0, D3, f, Ds3 = p
    D = D3 * _SCALE; Dstar = Ds3 * _SCALE
    e_slow = np.exp(-b * D)
    e_fast = np.exp(-b * (D + Dstar))
    S = s0 * ((1.0 - f) * e_slow + f * e_fast)
    dS0 = (1.0 - f) * e_slow + f * e_fast
    dD = -b * S
    df = s0 * (e_fast - e_slow)
    dDstar = -s0 * f * b * e_fast
    J = np.stack([dS0, dD * _SCALE, df, dDstar * _SCALE], axis=1)
    return J


def snr_to_sigma(snr, s0=1.0):
    """Noise std for a given SNR (signal-domain): sigma = S0 / SNR."""
    return s0 / float(snr)


def add_rician_noise(clean, sigma, rng):
    """Add Rician noise of std ``sigma`` to a clean magnitude signal (clean-room).

    A Rician magnitude is ``|S + n_r + i n_i|`` with ``n_r, n_i ~ N(0, sigma)``.
    """
    clean = np.asarray(clean, dtype=float)
    nr = rng.normal(0.0, sigma, size=clean.shape)
    ni = rng.normal(0.0, sigma, size=clean.shape)
    return np.sqrt((clean + nr) ** 2 + ni ** 2)
