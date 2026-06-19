"""Cohort families, ground-truth schema, and seeded generation for Lattice.

A :class:`Cohort` bundles ground-truth parameters with the multi-b signals
generated from them. The canonical parameter order is ``(D, Dstar, f)``; any
family-specific parameters (``k``, ``cv``, ``beta``, ``g`` ...) are carried
separately in ``extra`` and named in ``extra_names``.

Every cohort is fully reproducible from its integer ``seed`` alone: base
parameters, family-specific parameters, and noise are drawn from three
independent, offset RNG streams so that the base ``(D, Dstar, f)`` draws are
*identical across families* at a fixed seed. That is what makes the continuity
gates exact: a family at its continuity limit shares the bi-exponential cohort's
ground truth bit-for-bit.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from pathlib import Path

import numpy as np

from . import generators as G

__all__ = [
    "DEFAULT_BVALUES",
    "DEFAULT_SEED",
    "PARAM_NAMES",
    "PARAM_RANGES",
    "Cohort",
    "sample_params",
    "make_cohort",
    "continuity_residual",
    "DEFAULT_EXTRA",
    "CONTINUITY_EXTRA",
]

# 22-point b-value scheme: dense at low b (perfusion), sparse at high b (tissue).
DEFAULT_BVALUES = np.array(
    [0, 5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 120,
     150, 200, 250, 300, 400, 500, 600, 700, 800, 1000],
    dtype=float,
)

DEFAULT_SEED = 20260619  # reproducible-from-seed; date-stamped, no wall-clock.

PARAM_NAMES = ("D", "Dstar", "f")

# Published physiological ranges (mm^2/s for D, Dstar; dimensionless for f).
PARAM_RANGES = {
    "D": (0.5e-3, 3.0e-3),
    "Dstar": (10e-3, 100e-3),
    "f": (0.05, 0.40),
}

# RNG stream offsets (added to the cohort seed) so base/extra/noise are
# independent and base draws are family-invariant.
_OFFSET_BASE = 0
_OFFSET_EXTRA = 4040
_OFFSET_NOISE = 7000

# Default family-specific parameters (mid-deviation, clearly non-trivial).
DEFAULT_EXTRA = {
    "biexp": {},
    "dispersion_gamma": {"k": 4.0},          # CV = 1/sqrt(k) = 0.5
    "dispersion_lognormal": {"cv": 0.5},
    "stretched": {"beta": 0.7},
    "triexp": {"Dstar2_mult": 4.0, "g": 0.3},
}

# Family-specific parameters at the continuity limit (reduces to bi-exp).
CONTINUITY_EXTRA = {
    "biexp": {},
    "dispersion_gamma": {"k": 1e8},          # asymptotic; not exact at finite k
    "dispersion_lognormal": {"cv": 0.0},     # exact
    "stretched": {"beta": 1.0},              # exact
    "triexp": {"Dstar2_mult": 4.0, "g": 0.0},  # exact (g=0 drops the 3rd pool)
}

# Latent Gaussian-copula correlation for the "realistic" prior. Mild positive
# D-f coupling and negative f-Dstar coupling, consistent with abdominal/tumour
# IVIM literature; documented, not fit to any clinical data.
_REALISTIC_CORR = np.array([
    [1.0, -0.15, 0.25],   # D
    [-0.15, 1.0, -0.30],  # Dstar
    [0.25, -0.30, 1.0],   # f
])


def _uniform_marginal(u: np.ndarray, name: str) -> np.ndarray:
    lo, hi = PARAM_RANGES[name]
    return lo + (hi - lo) * u


def sample_params(n: int, rng: np.random.Generator, prior: str = "realistic") -> np.ndarray:
    """Draw ``n`` ground-truth ``(D, Dstar, f)`` rows.

    ``prior="uniform"`` draws each parameter independently uniform over its
    physiological range. ``prior="realistic"`` applies a Gaussian copula with a
    documented latent correlation, keeping the same uniform marginals.
    """
    if prior == "uniform":
        u = rng.random((n, 3))
    elif prior == "realistic":
        L = np.linalg.cholesky(_REALISTIC_CORR)
        z = rng.standard_normal((n, 3)) @ L.T
        # Standard-normal CDF -> uniform marginals (Gaussian copula).
        from math import sqrt
        u = 0.5 * (1.0 + _erf(z / sqrt(2.0)))
    else:
        raise ValueError(f"unknown prior {prior!r}")
    return np.column_stack([
        _uniform_marginal(u[:, 0], "D"),
        _uniform_marginal(u[:, 1], "Dstar"),
        _uniform_marginal(u[:, 2], "f"),
    ])


def _erf(x: np.ndarray) -> np.ndarray:
    """Vectorised erf via numpy (avoids a scipy dependency for the core path)."""
    # Abramowitz & Stegun 7.1.26, max abs error ~1.5e-7 -- ample for a prior.
    x = np.asarray(x, dtype=float)
    sign = np.sign(x)
    ax = np.abs(x)
    t = 1.0 / (1.0 + 0.3275911 * ax)
    y = 1.0 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t
                - 0.284496736) * t + 0.254829592) * t * np.exp(-ax * ax)
    return sign * y


@dataclass
class Cohort:
    """A synthetic IVIM ground-truth cohort.

    Attributes
    ----------
    family : str
        Generator family name (see :data:`lattice.generators.FAMILIES`).
    bvalues : (n_b,) ndarray
        Acquisition b-values (s/mm^2).
    params : (n, 3) ndarray
        Ground truth, columns ``(D, Dstar, f)`` (:data:`PARAM_NAMES`).
    extra : (n, n_extra) ndarray
        Family-specific ground-truth parameters (:attr:`extra_names`).
    signals_clean : (n, n_b) ndarray
        Noise-free forward signals (normalised, S0 = 1).
    signals : (n, n_b) ndarray
        Noisy observations at :attr:`snr` under :attr:`noise`.
    snr, noise, seed, prior : metadata.
    """

    family: str
    bvalues: np.ndarray
    params: np.ndarray
    extra: np.ndarray
    extra_names: tuple
    signals_clean: np.ndarray
    signals: np.ndarray
    snr: float
    noise: str
    seed: int
    prior: str
    param_names: tuple = PARAM_NAMES

    def __len__(self) -> int:
        return self.params.shape[0]

    @property
    def D(self) -> np.ndarray:
        return self.params[:, 0]

    @property
    def Dstar(self) -> np.ndarray:
        return self.params[:, 1]

    @property
    def f(self) -> np.ndarray:
        return self.params[:, 2]

    def manifest(self) -> dict:
        """Reproducibility manifest (everything needed to regenerate)."""
        return {
            "family": self.family,
            "n": len(self),
            "bvalues": self.bvalues.tolist(),
            "snr": self.snr,
            "noise": self.noise,
            "seed": self.seed,
            "prior": self.prior,
            "param_names": list(self.param_names),
            "extra_names": list(self.extra_names),
            "param_ranges": PARAM_RANGES,
        }

    def save(self, path: str | Path) -> Path:
        """Persist arrays to ``<path>.npz`` and the manifest to ``<path>.json``."""
        path = Path(path)
        np.savez_compressed(
            path.with_suffix(".npz"),
            bvalues=self.bvalues,
            params=self.params,
            extra=self.extra,
            signals_clean=self.signals_clean,
            signals=self.signals,
        )
        path.with_suffix(".json").write_text(json.dumps(self.manifest(), indent=2))
        return path.with_suffix(".npz")


def _forward(family: str, b, params: np.ndarray, extra_kw: dict) -> np.ndarray:
    """Evaluate a family's clean forward signal for an (n,3) param block."""
    D, Dstar, f = params[:, 0], params[:, 1], params[:, 2]
    if family == "biexp":
        return G.ivim_biexp(b, D, Dstar, f)
    if family == "dispersion_gamma":
        # mu plays the role of the effective Dstar.
        return G.ivim_dispersion_gamma(b, D, Dstar, extra_kw["k"], f)
    if family == "dispersion_lognormal":
        return G.ivim_dispersion_lognormal(b, D, Dstar, extra_kw["cv"], f)
    if family == "stretched":
        return G.ivim_stretched(b, D, Dstar, f, extra_kw["beta"])
    if family == "triexp":
        Dstar2 = extra_kw["Dstar2_mult"] * Dstar
        return G.ivim_triexp(b, D, Dstar, f, Dstar2, extra_kw["g"])
    raise ValueError(f"unknown family {family!r}")


def _extra_block(family: str, extra_kw: dict, n: int) -> tuple[np.ndarray, tuple]:
    names = G.FAMILIES[family]["extra"]
    if not names:
        return np.empty((n, 0)), ()
    cols = [np.full(n, float(extra_kw[name])) for name in names]
    return np.column_stack(cols), tuple(names)


def make_cohort(
    family: str = "biexp",
    n: int = 2000,
    snr: float = 50.0,
    seed: int = DEFAULT_SEED,
    prior: str = "realistic",
    noise: str = "rician",
    bvalues: np.ndarray | None = None,
    extra: dict | None = None,
) -> Cohort:
    """Generate a reproducible cohort for ``family``.

    The base ``(D, Dstar, f)`` draws depend only on ``(seed, prior, n)`` and are
    therefore identical across families -- the property the continuity gates
    rely on. ``extra`` overrides the family defaults (:data:`DEFAULT_EXTRA`).
    """
    if family not in G.FAMILIES:
        raise ValueError(f"unknown family {family!r}; choose from {list(G.FAMILIES)}")
    b = DEFAULT_BVALUES if bvalues is None else np.asarray(bvalues, dtype=float)
    extra_kw = dict(DEFAULT_EXTRA[family])
    if extra:
        extra_kw.update(extra)

    rng_base = np.random.default_rng(seed + _OFFSET_BASE)
    params = sample_params(n, rng_base, prior=prior)

    clean = _forward(family, b, params, extra_kw)
    extra_block, extra_names = _extra_block(family, extra_kw, n)

    rng_noise = np.random.default_rng(seed + _OFFSET_NOISE)
    if noise == "rician":
        noisy = G.add_rician_noise(clean, snr, rng_noise)
    elif noise == "gaussian":
        noisy = G.add_gaussian_noise(clean, snr, rng_noise)
    elif noise == "none":
        noisy = clean.copy()
    else:
        raise ValueError(f"unknown noise {noise!r}")

    return Cohort(
        family=family,
        bvalues=b,
        params=params,
        extra=extra_block,
        extra_names=extra_names,
        signals_clean=clean,
        signals=noisy,
        snr=float(snr),
        noise=noise,
        seed=int(seed),
        prior=prior,
    )


def continuity_residual(
    family: str,
    n: int = 2000,
    seed: int = DEFAULT_SEED,
    prior: str = "realistic",
    bvalues: np.ndarray | None = None,
) -> float:
    """Max abs deviation between ``family`` at its continuity limit and bi-exp.

    Uses identical base ``(D, Dstar, f)`` draws for both, so the result is the
    pure forward-model continuity error. Exact-reduction families
    (``dispersion_lognormal`` at cv=0, ``stretched`` at beta=1, ``triexp`` at
    g=0) return 0.0; ``dispersion_gamma`` returns a small asymptotic residual.
    """
    b = DEFAULT_BVALUES if bvalues is None else np.asarray(bvalues, dtype=float)
    rng_base = np.random.default_rng(seed + _OFFSET_BASE)
    params = sample_params(n, rng_base, prior=prior)
    biexp = _forward("biexp", b, params, {})
    limit = _forward(family, b, params, CONTINUITY_EXTRA[family])
    return float(np.max(np.abs(limit - biexp)))
