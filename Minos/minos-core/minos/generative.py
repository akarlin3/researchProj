"""Synthetic latent field + measurement model, with common-random-number (CRN)
realisation so every ``(tau, delta)`` reuses one set of base variates.

The latent severity ``theta`` and the measurement map live behind a small seam so a
real IVIM parameter map + Fashion posterior can later replace the synthetic source
without touching the decision / VoI / gate core. See the marked region below.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import MinosConfig


@dataclass(frozen=True)
class BaseDraws:
    """CRN base variates, sampled once and reused across the whole sweep.

    ``theta`` is the latent severity (independent of ``tau``/``delta``); ``z_eta`` and
    ``z_w`` are standard-normal base draws transformed per shift in :func:`realise`.
    """

    theta: np.ndarray
    z_eta: np.ndarray
    z_w: np.ndarray


# ----------------------------------------------------------------------------------
# IVIM seam — Fashion integration point (deferred).
#
# ``sample_latent`` is the ONLY place the synthetic latent source is defined. To wire
# in real data later, replace its body so it returns a per-voxel IVIM parameter map
# (e.g. pseudo-diffusion fraction f) as ``theta``, and have :func:`realise` consume a
# Fashion posterior ``(mu, sigma_rep)`` instead of the synthetic Gaussian measurement.
# Nothing in decision.py / voi.py / gate.py reads the latent source directly.
# ----------------------------------------------------------------------------------
def sample_latent(cfg: MinosConfig, rng: np.random.Generator, n: int) -> np.ndarray:
    """Sample ``n`` latent severities from the 3-component Gaussian mixture prior."""
    weights = np.asarray(cfg.mix_weights, dtype=float)
    means = np.asarray(cfg.mix_means, dtype=float)
    stds = np.asarray(cfg.mix_stds, dtype=float)
    comp = rng.choice(len(weights), size=n, p=weights)
    z = rng.standard_normal(n)
    return means[comp] + stds[comp] * z


def make_population(cfg: MinosConfig, rng: np.random.Generator) -> BaseDraws:
    """Draw the CRN base variates for ``cfg.n_voxels`` voxels."""
    n = cfg.n_voxels
    theta = sample_latent(cfg, rng, n)
    z_eta = rng.standard_normal(n)
    z_w = rng.standard_normal(n)
    return BaseDraws(theta=theta, z_eta=z_eta, z_w=z_w)


def _as_shift_mask(shift, n: int) -> np.ndarray:
    if np.isscalar(shift) or (isinstance(shift, np.ndarray) and shift.ndim == 0):
        return np.full(n, bool(shift))
    shift = np.asarray(shift, dtype=bool)
    assert shift.shape == (n,), "shift mask must be scalar or length n_voxels"
    return shift


def realise(base: BaseDraws, cfg: MinosConfig, *, delta: float, shift):
    """Map base variates to observed estimate ``mu`` and acquisition feature ``w``.

    ``shift`` is a bool scalar (applies to all voxels) or a per-voxel mask. Shifted
    voxels get inflated true noise ``s*(1+alpha*delta)``, a downward bias
    ``-beta*s*delta``, and feature mean ``delta``; unshifted voxels are nominal.
    """
    n = base.theta.shape[0]
    mask = _as_shift_mask(shift, n)
    sigma_true = np.where(mask, cfg.s * (1.0 + cfg.alpha * delta), cfg.s)
    bias = np.where(mask, -cfg.beta * cfg.s * delta, 0.0)
    feat_mean = np.where(mask, delta, 0.0)
    mu = base.theta + bias + sigma_true * base.z_eta
    w = cfg.w_train_mean + feat_mean + cfg.w_train_std * base.z_w
    return mu, w
