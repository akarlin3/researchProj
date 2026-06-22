"""The twin's synthetic IVIM scanner: biexponential signal + Rician noise.

Pure-numeric, vectorised over voxels. This is the *forward* model of the synthetic
twin only — it is NOT a scanner driver and touches no real acquisition or patient data.
"""
from __future__ import annotations

import numpy as np


def ivim_signal(bvals, D, Dstar, f, S0=1.0):
    """Biexponential IVIM signal.

    Returns an array of shape ``(n_b, V)`` for b-values ``bvals`` (n_b,) and per-voxel
    parameter arrays ``D, Dstar, f`` (each length V or scalar).
    """
    b = np.asarray(bvals, float)[:, None]
    D = np.atleast_1d(np.asarray(D, float))[None, :]
    Dstar = np.atleast_1d(np.asarray(Dstar, float))[None, :]
    f = np.atleast_1d(np.asarray(f, float))[None, :]
    S0 = np.asarray(S0, float)
    return S0 * (f * np.exp(-b * (D + Dstar)) + (1.0 - f) * np.exp(-b * D))


def add_rician_noise(signal, snr, rng, S0=1.0):
    """Add Rician noise at the given SNR (relative to ``S0``) using ``rng``."""
    sigma = float(S0) / float(snr)
    real = signal + rng.normal(0.0, sigma, size=signal.shape)
    imag = rng.normal(0.0, sigma, size=signal.shape)
    return np.sqrt(real * real + imag * imag)


def simulate_scan(bvals, D, Dstar, f, snr, n_noise, rng, S0=1.0):
    """One synthetic scan: ``n_noise`` Rician realisations of the IVIM signal.

    ``snr`` may be a scalar or a per-voxel array (length V) — a per-voxel SNR map
    models a realistic abdominal acquisition-degradation (motion/susceptibility) zone.
    Returns ``(n_b, V, n_noise)``. The mean over the noise axis is the measured signal;
    the spread over the axis is what the fit turns into a raw posterior.
    """
    clean = ivim_signal(bvals, D, Dstar, f, S0=S0)              # (n_b, V)
    n_b, V = clean.shape
    reps = clean[:, :, None] + np.zeros((1, 1, n_noise))        # broadcast (n_b, V, K)
    snr = np.atleast_1d(np.asarray(snr, float))
    sig = float(S0) / (snr if snr.size == V else np.full(V, snr.reshape(-1)[0]))
    sig = sig[None, :, None]                                    # (1, V, 1)
    real = reps + rng.normal(0.0, 1.0, size=reps.shape) * sig
    imag = rng.normal(0.0, 1.0, size=reps.shape) * sig
    return np.sqrt(real * real + imag * imag)
