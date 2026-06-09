"""
Synthetic in-distribution (ID) reference at the ACRIN-6698 acquisition scheme.

Arm 1 of Tier 2 is *estimator-free*: we build the ID distribution directly in the
acquired 4-b-value signal space (b = 0, 100, 600, 800), with the same breast-tissue
IVIM priors used in Tier 1 and Rician noise at an ACRIN-representative SNR. Real
ACRIN units are then scored against this synthetic ID distribution.

Normalisation note (stated explicitly, as the directive asks for the Family-1
asymmetry): in vivo the DWI signal is only known up to the unknown S0 scaling, so
every unit is normalised by its own measured b=0 signal, S(b)/S(b=0). After this
scale normalisation the b=0 entry is identically 1 and carries no information, so
the operative feature space is the 3 informative decay ratios
[S(100)/S(0), S(600)/S(0), S(800)/S(0)]. We keep the "4-b normalized-signal space"
framing but flag that it is 3-D once S0 is removed -- there is no honest way to get
a 4th independent dimension out of a single-S0 acquisition.
"""

import numpy as np
import torch

from sibyl.data.synthetic import sample_prior, PRIOR_BOUNDS
from sibyl.data.shift import apply_noise
from sibyl.forward_model.ivim import ivim_biexponential, ACRIN_B_SCHEME

# ACRIN DWI is acquired at 1.5T/3T with NEX 2-5; there is *no* published numeric
# SNR (it is only a categorical QA rating). We treat the ID-reference SNR as a
# free parameter to be matched to / swept around the empirical b0 SNR of the real
# data. This default is a plausible mid-range value, used for synthetic validation.
ACRIN_REF_SNR = 40.0


def normalize_by_b0(signal_4b: np.ndarray, b0_index: int = 0, eps: float = 1e-6) -> np.ndarray:
    """
    Divide each 4-b signal vector by its own b=0 entry, reproducing the in-vivo
    situation where only S(b)/S0 ratios are observable.

    Parameters
    ----------
    signal_4b : np.ndarray
        Raw signal (shape: [N, 4] or [4]).

    Returns
    -------
    np.ndarray
        b0-normalised signal (shape: [N, 4]); column ``b0_index`` is ~1.
    """
    signal_4b = np.atleast_2d(np.asarray(signal_4b, dtype=np.float64))
    b0 = signal_4b[:, b0_index:b0_index + 1]
    return signal_4b / np.clip(b0, eps, None)


def signal_features(signal_4b_norm: np.ndarray, b0_index: int = 0) -> np.ndarray:
    """
    Drop the (constant ~1) b=0 column to give the 3 informative decay ratios that
    the estimator-free detector actually operates on.

    Parameters
    ----------
    signal_4b_norm : np.ndarray
        b0-normalised signal (shape: [N, 4]).

    Returns
    -------
    np.ndarray
        Decay-ratio features (shape: [N, 3]).
    """
    signal_4b_norm = np.atleast_2d(np.asarray(signal_4b_norm, dtype=np.float64))
    keep = [j for j in range(signal_4b_norm.shape[1]) if j != b0_index]
    return signal_4b_norm[:, keep]


def generate_acrin_id_reference(
    n_samples: int,
    snr: float = ACRIN_REF_SNR,
    seed: int = None,
    rician: bool = True,
):
    """
    Generate the synthetic ID reference at the ACRIN 4-b scheme.

    Returns
    -------
    theta : torch.Tensor
        IVIM parameters [f, D, Dstar] (shape: [N, 3]).
    raw_4b : np.ndarray
        Noisy raw signal at b=[0,100,600,800] with S0=1 (shape: [N, 4]).
    features : np.ndarray
        b0-normalised decay-ratio features (shape: [N, 3]).
    """
    theta = sample_prior(n_samples, seed=seed)
    f, D, Dstar = theta[:, 0], theta[:, 1], theta[:, 2]
    clean = ivim_biexponential(ACRIN_B_SCHEME, f, D, Dstar)  # [N, 4], S0=1
    noisy = apply_noise(clean, snr=snr, rician=rician, seed=seed)
    raw_4b = noisy.cpu().numpy().astype(np.float64)
    feats = signal_features(normalize_by_b0(raw_4b))
    return theta, raw_4b, feats
