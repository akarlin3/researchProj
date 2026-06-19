"""Clean-room IVIM forward signal generators for the Lattice DRO.

Every generator returns the *normalised* signal S(b)/S0 (set ``s0`` to scale).
All five families share one invariant: at their **continuity limit** they reduce
to the canonical bi-exponential model. Three of them reduce *exactly*
(lognormal dispersion at ``cv=0``, stretched at ``beta=1``, tri-exp at ``g=0``);
the gamma-dispersion family reduces only asymptotically as ``k -> inf``.

These are re-implementations of the standard IVIM signal equations (physics, not
borrowed code). No clinical, scanner, or third-party data is involved.

References for the equations
----------------------------
* Bi-exponential IVIM: Le Bihan et al., Radiology 1988.
* Gamma velocity dispersion: the perfusion term is the Laplace transform of a
  Gamma(shape=k, mean=mu) pseudo-diffusivity, ``(1 + b*mu/k)**(-k)``.
* Log-normal velocity dispersion: no closed form; evaluated by Gauss-Hermite
  quadrature over a LogNormal(mean=mu, cv) pseudo-diffusivity.
* Stretched-exponential (anomalous) perfusion: ``exp(-(b*Dstar)**beta)``.
* Tri-exponential: a third, fast pseudo-diffusion compartment.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "ivim_biexp",
    "ivim_dispersion_gamma",
    "ivim_dispersion_lognormal",
    "ivim_stretched",
    "ivim_triexp",
    "add_rician_noise",
    "add_gaussian_noise",
    "FAMILIES",
]


def _asarray(x) -> np.ndarray:
    return np.asarray(x, dtype=float)


def ivim_biexp(b, D, Dstar, f, s0: float = 1.0):
    """Canonical bi-exponential IVIM signal.

    ``S(b)/S0 = f * exp(-b*Dstar) + (1 - f) * exp(-b*D)``.

    All scalar parameters broadcast against ``b``; pass arrays of equal length
    for a vectorised cohort (``b`` shape ``(n_b,)``, params shape ``(n,)``
    produce an ``(n, n_b)`` signal).
    """
    b = _asarray(b)
    D = _asarray(D)[..., None] if np.ndim(D) else D
    Dstar = _asarray(Dstar)[..., None] if np.ndim(Dstar) else Dstar
    f = _asarray(f)[..., None] if np.ndim(f) else f
    return s0 * (f * np.exp(-b * Dstar) + (1.0 - f) * np.exp(-b * D))


def ivim_dispersion_gamma(b, D, mu, k, f, s0: float = 1.0):
    """Gamma velocity-dispersion IVIM.

    The pseudo-diffusivity is ``D* ~ Gamma(shape=k, mean=mu)``; its Laplace
    transform gives the perfusion term ``(1 + b*mu/k)**(-k)``. The effective
    (low-b log-slope) pseudo-diffusion is ``mu``. As ``k -> inf`` the perfusion
    term -> ``exp(-b*mu)`` and the model reduces to bi-exponential with
    ``Dstar = mu`` (asymptotic, not exact at finite k).
    """
    b = _asarray(b)
    D = _asarray(D)[..., None] if np.ndim(D) else D
    mu = _asarray(mu)[..., None] if np.ndim(mu) else mu
    k = _asarray(k)[..., None] if np.ndim(k) else k
    f = _asarray(f)[..., None] if np.ndim(f) else f
    perf = np.power(1.0 + b * mu / k, -k)
    return s0 * (f * perf + (1.0 - f) * np.exp(-b * D))


def ivim_dispersion_lognormal(b, D, mu, cv, f, s0: float = 1.0, n_quad: int = 64):
    """Log-normal velocity-dispersion IVIM.

    ``D* ~ LogNormal`` with mean ``mu`` and coefficient of variation ``cv``.
    The perfusion term ``E[exp(-b*D*)]`` has no closed form and is evaluated by
    ``n_quad``-point Gauss-Hermite quadrature. At ``cv = 0`` the distribution
    collapses to a point mass at ``mu`` and the model reduces *exactly* to
    bi-exponential with ``Dstar = mu``.
    """
    b = _asarray(b)
    D = _asarray(D)
    mu = _asarray(mu)
    cv = _asarray(cv)
    f = _asarray(f)

    # Broadcast params to a common (n,) shape, then add b on the last axis.
    D, mu, cv, f = np.broadcast_arrays(D, mu, cv, f)
    scalar = D.ndim == 0
    D = np.atleast_1d(D)[:, None]
    mu = np.atleast_1d(mu)[:, None]
    cv = np.atleast_1d(cv)[:, None]
    f = np.atleast_1d(f)[:, None]

    # LogNormal(mean=mu, cv): sigma^2 = ln(1+cv^2), m = ln(mu) - sigma^2/2.
    sigma = np.sqrt(np.log1p(cv ** 2))
    m = np.log(np.where(mu > 0, mu, 1.0)) - 0.5 * sigma ** 2

    nodes, weights = np.polynomial.hermite_e.hermegauss(n_quad)  # physicists' "probabilists" form
    w = (weights / np.sqrt(2.0 * np.pi))[None, None, :]
    z = nodes[None, None, :]
    # D*_i = exp(m + sigma * z_i); perfusion = sum_i w_i exp(-b D*_i)
    dstar = np.exp(m[..., None] + sigma[..., None] * z)  # (n, 1, n_quad)
    kernel = np.exp(-b[None, :, None] * dstar)           # (n, n_b, n_quad)
    perf = np.sum(w * kernel, axis=-1)                   # (n, n_b)

    # cv == 0 -> degenerate: use exact point mass exp(-b*mu).
    degenerate = (cv[..., 0] == 0.0)
    if np.any(degenerate):
        exact = np.exp(-b[None, :] * mu)
        perf = np.where(degenerate[:, None], exact, perf)

    sig = s0 * (f * perf + (1.0 - f) * np.exp(-b[None, :] * D))
    return sig[0] if scalar else sig


def ivim_stretched(b, D, Dstar, f, beta, s0: float = 1.0):
    """Stretched-exponential (anomalous) perfusion IVIM.

    ``S(b)/S0 = (1 - f) exp(-b*D) + f exp(-(b*Dstar)**beta)``.
    At ``beta = 1`` it reduces *exactly* to bi-exponential.
    """
    b = _asarray(b)
    D = _asarray(D)[..., None] if np.ndim(D) else D
    Dstar = _asarray(Dstar)[..., None] if np.ndim(Dstar) else Dstar
    f = _asarray(f)[..., None] if np.ndim(f) else f
    beta = _asarray(beta)[..., None] if np.ndim(beta) else beta
    perf = np.exp(-np.power(np.clip(b * Dstar, 0.0, None), beta))
    return s0 * (f * perf + (1.0 - f) * np.exp(-b * D))


def ivim_triexp(b, D, Dstar, f, Dstar2, g, s0: float = 1.0):
    """Tri-exponential IVIM: bi-exp plus a third, fast pseudo-diffusion pool.

    ``S(b)/S0 = (1-f) e^{-bD} + f(1-g) e^{-bD*} + f g e^{-bD*2}``.
    ``g`` in [0, 1] is the fraction of the perfusion pool routed to the fast
    component. At ``g = 0`` it reduces *exactly* to bi-exponential.
    """
    b = _asarray(b)
    D = _asarray(D)[..., None] if np.ndim(D) else D
    Dstar = _asarray(Dstar)[..., None] if np.ndim(Dstar) else Dstar
    f = _asarray(f)[..., None] if np.ndim(f) else f
    Dstar2 = _asarray(Dstar2)[..., None] if np.ndim(Dstar2) else Dstar2
    g = _asarray(g)[..., None] if np.ndim(g) else g
    return s0 * (
        (1.0 - f) * np.exp(-b * D)
        + f * (1.0 - g) * np.exp(-b * Dstar)
        + f * g * np.exp(-b * Dstar2)
    )


def add_gaussian_noise(signal, snr: float, rng: np.random.Generator, s0: float = 1.0):
    """Additive Gaussian noise, sigma = s0 / snr (defined at b=0)."""
    sigma = s0 / float(snr)
    return signal + rng.normal(0.0, sigma, size=signal.shape)


def add_rician_noise(signal, snr: float, rng: np.random.Generator, s0: float = 1.0):
    """Rician-distributed magnitude noise, sigma = s0 / snr (defined at b=0).

    ``mag = sqrt((S + n_re)^2 + n_im^2)`` with ``n_re, n_im ~ N(0, sigma)``.
    """
    sigma = s0 / float(snr)
    n_re = rng.normal(0.0, sigma, size=signal.shape)
    n_im = rng.normal(0.0, sigma, size=signal.shape)
    return np.sqrt((signal + n_re) ** 2 + n_im ** 2)


# Registry of generator families, keyed by the name used in ``make_cohort``.
# ``extra`` lists the family-specific parameters appended after (D, Dstar, f);
# ``continuity`` documents the exact/asymptotic reduction to bi-exponential.
FAMILIES = {
    "biexp": {
        "fn": ivim_biexp,
        "extra": (),
        "continuity": "identity",
    },
    "dispersion_gamma": {
        "fn": ivim_dispersion_gamma,
        "extra": ("k",),
        "continuity": "k -> inf (asymptotic)",
    },
    "dispersion_lognormal": {
        "fn": ivim_dispersion_lognormal,
        "extra": ("cv",),
        "continuity": "cv = 0 (exact)",
    },
    "stretched": {
        "fn": ivim_stretched,
        "extra": ("beta",),
        "continuity": "beta = 1 (exact)",
    },
    "triexp": {
        "fn": ivim_triexp,
        "extra": ("Dstar2_mult", "g"),
        "continuity": "g = 0 (exact)",
    },
}
