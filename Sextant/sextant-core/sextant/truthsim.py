"""Truth-controlled IVIM simulation for the HC2/CS2 robustness battery.

The boundary-railing primary claim (Sextant) is read off the optimiser and needs
no ground truth.  A sceptic can still ask a *different* question of it: is the
railing we see on the real abdominal scan an intrinsic pathology of the bounded
NLLS estimator under weak D* identifiability, or could it be an artefact of
*simulator-reality mismatch* (the forward / noise model not matching tissue)?

This module supplies the truth-controlled substrate that answers it.  It

  * draws ground-truth (D, D*, f) over the trained NPE BoxUniform prior,
  * renders signals under a *well-specified* bi-exponential forward model and
    under three deliberate *misspecifications* (tri-exponential tissue tail,
    diffusional kurtosis, and a non-Rician multi-coil noise model), and
  * fits every voxel with **Fashion's exact bounded bi-exponential NLLS** and
    **Fashion's exact rail thresholds**, reused read-only via
    :mod:`sextant.fashion_reuse` (no reimplementation, no silent divergence).

Because the truth is known at every voxel, the railing rate and the
parameter-recovery error can both be reported against ground truth -- the
in-silico misspecification isolation (HC2 item b) and the known-truth recovery
test (HC2 item a) the railing-first paper needs.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .fashion_reuse import load_railing

# ---- canonical analysis constants, reused read-only from Fashion ----------- #
_R = load_railing()
TARGET_BVALS: np.ndarray = np.asarray(_R["TARGET_BVALS"], float)   # 10-pt OSIPI scheme
DSTAR_LOWER_RAIL: float = float(_R["DSTAR_LOWER_RAIL"])            # 0.0033 mm^2/s
DSTAR_UPPER_RAIL: float = float(_R["DSTAR_UPPER_RAIL"])           # 0.1485 mm^2/s
_FIT = _R["fit_biexp_nlls"]                                        # bounded TRF NLLS

# Trained-NPE BoxUniform prior (display units -> mm^2/s); D* is log-uniform, as
# trained ("D* trained in log space").  Matches paper/manuscript.tex Methods.
PRIOR = {
    "D":  (0.2e-3, 3.0e-3),     # uniform
    "Ds": (3.0e-3, 150.0e-3),   # log-uniform
    "f":  (0.0, 0.5),           # uniform
}


def draw_truths(n: int, rng: np.random.Generator) -> np.ndarray:
    """Draw ``n`` ground-truth (D, D*, f) rows over the trained NPE prior."""
    D = rng.uniform(*PRIOR["D"], size=n)
    lo, hi = np.log(PRIOR["Ds"][0]), np.log(PRIOR["Ds"][1])
    Ds = np.exp(rng.uniform(lo, hi, size=n))
    f = rng.uniform(*PRIOR["f"], size=n)
    return np.column_stack([D, Ds, f])


# --------------------------------------------------------------------------- #
# Forward models: one well-specified, three misspecified.  All return the clean
# (noise-free) normalised signal S(b)/S0 at the given b-values, S0 == 1.
# --------------------------------------------------------------------------- #
def fm_biexp(b, D, Ds, f):
    """Well-specified bi-exponential IVIM (the model the NLLS assumes)."""
    return f * np.exp(-b * Ds) + (1.0 - f) * np.exp(-b * D)


def fm_triexp(b, D, Ds, f, w=0.25, slow=0.35):
    """Misspecified: tissue compartment carries a slow restricted tail.

    The single tissue exponential is replaced by a two-component mixture whose
    fast part is D and whose slow part is ``slow*D`` (weight ``w``); the bounded
    bi-exponential fit has no degree of freedom for the extra curvature.
    """
    tissue = (1.0 - w) * np.exp(-b * D) + w * np.exp(-b * slow * D)
    return f * np.exp(-b * Ds) + (1.0 - f) * tissue


def fm_kurtosis(b, D, Ds, f, K=0.8):
    """Misspecified: non-Gaussian (diffusional-kurtosis) tissue diffusion."""
    tissue = np.exp(-b * D + (1.0 / 6.0) * (b * D) ** 2 * K)
    return f * np.exp(-b * Ds) + (1.0 - f) * tissue


FORWARD_MODELS = {
    "biexp_WS": fm_biexp,
    "triexp": fm_triexp,
    "kurtosis": fm_kurtosis,
    # "noise_chi" shares fm_biexp clean signal but a non-Rician noise model
    # (handled in add_noise); registered here so callers can iterate uniformly.
    "noise_chi": fm_biexp,
}


def add_noise(clean: np.ndarray, snr: float, rng: np.random.Generator,
              model: str = "rician") -> np.ndarray:
    """Add magnitude noise at the given SNR (sigma = 1/SNR, since S0 == 1).

    ``rician`` -- standard single-coil magnitude noise.
    ``chi``    -- L=4 coil sum-of-squares (non-central chi); a *noise-model*
                  misspecification relative to the Rician assumption.
    """
    sigma = 1.0 / snr
    if model == "rician":
        re = clean + rng.normal(0.0, sigma, size=clean.shape)
        im = rng.normal(0.0, sigma, size=clean.shape)
        return np.sqrt(re ** 2 + im ** 2)
    if model == "chi":
        L = 4
        acc = np.zeros_like(clean)
        # distribute the signal power across coils, sigma per channel scaled so
        # the effective magnitude SNR matches `snr`.
        s_ch = clean / np.sqrt(L)
        sig_ch = sigma  # per-channel
        for _ in range(L):
            re = s_ch + rng.normal(0.0, sig_ch, size=clean.shape)
            im = rng.normal(0.0, sig_ch, size=clean.shape)
            acc += re ** 2 + im ** 2
        return np.sqrt(acc)
    raise ValueError(f"unknown noise model {model!r}")


def render(truths: np.ndarray, bvals: np.ndarray, snr: float,
           rng: np.random.Generator, forward: str = "biexp_WS") -> np.ndarray:
    """Render noisy normalised signals for every truth row under one forward model."""
    fm = FORWARD_MODELS[forward]
    D, Ds, f = truths[:, 0], truths[:, 1], truths[:, 2]
    clean = fm(bvals[None, :], D[:, None], Ds[:, None], f[:, None])
    noise = "chi" if forward == "noise_chi" else "rician"
    return add_noise(clean, snr, rng, model=noise)


@dataclass
class FitResult:
    dstar: np.ndarray       # (N,) fitted D*
    D: np.ndarray
    f: np.ndarray
    railed: np.ndarray      # (N,) bool
    railed_lo: np.ndarray
    railed_hi: np.ndarray


def fit_and_rail(signals: np.ndarray, bvals: np.ndarray) -> FitResult:
    """Fit every voxel with Fashion's bounded NLLS; flag rail at Fashion's thresholds."""
    out = np.array([_FIT(bvals, signals[i]) for i in range(len(signals))])
    D, Ds, f = out[:, 0], out[:, 1], out[:, 2]
    finite = np.isfinite(Ds)
    lo = finite & (Ds <= DSTAR_LOWER_RAIL)
    hi = finite & (Ds >= DSTAR_UPPER_RAIL)
    return FitResult(dstar=Ds, D=D, f=f, railed=(lo | hi),
                     railed_lo=lo, railed_hi=hi)
