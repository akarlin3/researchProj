"""
Arm 2 imputation: lift a 4-b ACRIN unit onto the dense 10-b grid so the dense
uGUIDE estimator (trained in Tier 1) can be applied.

Per unit we fit a quick segmented IVIM estimate from the 4 real points, then
forward-project that IVIM model onto the 6 missing dense-grid b-values and splice
the 4 real (normalised) values back in at their grid positions.

Stated circularity (this is exactly why Arm 2 is SECONDARY): the 6 imputed points
are generated *by the IVIM forward model*. They therefore carry the model
assumption, so a low OOD score at those positions can be the imputer fitting the
model rather than the unit being genuinely in-distribution. Agreement between the
Arm-1 (estimator-free) and Arm-2 scores is the check that this circularity is not
driving the result; divergence flags imputation-injected IVIM artifacts.
"""

import numpy as np
import torch

from sibyl.forward_model.ivim import ivim_biexponential, DEFAULT_B_SCHEME, ACRIN_B_SCHEME, ACRIN_TO_DENSE_IDX
from sibyl.data.synthetic import PRIOR_BOUNDS


def segmented_ivim_fit(sig4_norm: np.ndarray, eps: float = 1e-6):
    """
    Fast segmented IVIM fit from the 4 b0-normalised ACRIN points.

    Step 1: high-b (b>=100) log-linear fit -> D (slope) and (1-f) (intercept),
            since the perfusion term is ~0 there.
    Step 2: low-b residual at b=100 -> D* (poorly constrained from one point;
            falls back to the prior-mean D* where the residual is non-physical).

    Parameters
    ----------
    sig4_norm : np.ndarray
        b0-normalised 4-b signal [N, 4] at b = [0, 100, 600, 800] (col 0 ~ 1).

    Returns
    -------
    f, D, Dstar : np.ndarray, each [N]
    """
    sig4_norm = np.atleast_2d(np.asarray(sig4_norm, dtype=np.float64))
    bvals = ACRIN_B_SCHEME.numpy().astype(np.float64)

    # --- Step 1: high-b log-linear fit over b = 100, 600, 800 ---
    hi = bvals >= 100
    b_hi = bvals[hi]                       # [3]
    y = np.log(np.clip(sig4_norm[:, hi], eps, None))  # [N, 3]
    bbar = b_hi.mean()
    bc = b_hi - bbar                       # [3]
    denom = np.sum(bc ** 2)
    yc = y - y.mean(axis=1, keepdims=True)
    slope = np.sum(yc * bc, axis=1) / denom  # [N]; slope = -D
    D = np.clip(-slope, 1e-4, 5e-3)
    intercept = y.mean(axis=1) - slope * bbar
    one_minus_f = np.clip(np.exp(intercept), 1e-3, 1.0)
    f = np.clip(1.0 - one_minus_f, 0.0, 0.6)

    # --- Step 2: D* from the b=100 perfusion residual ---
    s100 = sig4_norm[:, 1]
    tissue_100 = one_minus_f * np.exp(-100.0 * D)
    perf_100 = s100 - tissue_100                       # = f * exp(-100 (D + D*))
    dstar_prior = float(np.mean(PRIOR_BOUNDS["Dstar"]))
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = perf_100 / np.clip(f, 1e-3, None)      # = exp(-100 (D + D*))
        Dstar = -np.log(np.clip(ratio, 1e-6, 1.0)) / 100.0 - D
    bad = ~np.isfinite(Dstar) | (perf_100 <= 0) | (f < 1e-2)
    Dstar = np.where(bad, dstar_prior, Dstar)
    Dstar = np.clip(Dstar, 3e-3, 100e-3)
    return f, D, Dstar


def impute_dense(sig4_norm: np.ndarray) -> np.ndarray:
    """
    Build the length-10 dense-grid signal for each unit: forward-project the
    segmented IVIM fit onto DEFAULT_B_SCHEME, then splice the 4 real normalised
    values back in at their dense-grid indices.

    Returns
    -------
    np.ndarray
        Dense normalised signal [N, 10] suitable as a uGUIDE input.
    """
    sig4_norm = np.atleast_2d(np.asarray(sig4_norm, dtype=np.float64))
    f, D, Dstar = segmented_ivim_fit(sig4_norm)
    dense = ivim_biexponential(
        DEFAULT_B_SCHEME,
        torch.tensor(f, dtype=torch.float32),
        torch.tensor(D, dtype=torch.float32),
        torch.tensor(Dstar, dtype=torch.float32),
    ).cpu().numpy().astype(np.float64)  # [N, 10], S0=1

    # Splice the real measured values back at the ACRIN positions.
    for acrin_j, dense_idx in enumerate(ACRIN_TO_DENSE_IDX):
        dense[:, dense_idx] = sig4_norm[:, acrin_j]
    return dense
