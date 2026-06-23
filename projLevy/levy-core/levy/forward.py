"""The diffusion-MRI signal-decay forward model S(b; theta) and its Jacobian.

LEAD LANE -- stretched-exponential (single fractional order alpha):

    S(b; S0, D, alpha) = S0 * exp( -(b * D)^alpha )                          (Bennett 2003)

with alpha in (0, 1] the stretching / intra-voxel-heterogeneity exponent, D the
diffusion coefficient (units mm^2/s), b the diffusion weighting (units s/mm^2), S0
the non-diffusion-weighted amplitude. alpha = 1 recovers the mono-exponential
S0 * exp(-b D).

This is the object whose identifiability Levy bounds. Crucially, the parameter enters
the likelihood ONLY through the b-indexed *signal attenuation* S(b; theta) -- NOT through
a trajectory / mean-squared-displacement / increment process. That is the forward-model
difference from the fBm/Hurst CRB of Coeurjolly-Istas (2001): same "fractional exponent"
word, different statistical experiment.

The joint CTRW / fractional Bloch-Torrey two-exponent model

    S(b; S0, D, alpha, beta) = S0 * exp( -(b * D)^alpha )   [time-alpha lane; beta enters
                                                             the spatial/space-fractional
                                                             attenuation -- Phase 3]

is stubbed via ``forward_joint`` for the Phase-3 (alpha, beta) degeneracy work and is
NOT exercised by CP0.

PARAMETER ORDER convention (used by the Jacobian, the Fisher matrix, and CRLB indexing):
    theta = (S0, D, alpha)        indices 0, 1, 2
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

PARAM_NAMES = ("S0", "D", "alpha")
IDX = {name: i for i, name in enumerate(PARAM_NAMES)}


@dataclass(frozen=True)
class StretchedExp:
    """Stretched-exponential ground-truth parameters and a fixed b-value design.

    Parameters
    ----------
    S0, D, alpha:
        Ground-truth parameters. D in mm^2/s (e.g. 1.5e-3 for tissue), alpha in (0, 1].
    """

    S0: float = 1.0
    D: float = 1.5e-3
    alpha: float = 0.75

    @property
    def theta(self) -> np.ndarray:
        return np.array([self.S0, self.D, self.alpha], dtype=float)


def signal(b, theta):
    """S(b; theta) = S0 * exp(-(b D)^alpha).  ``b`` array-like, ``theta`` = (S0, D, alpha)."""
    b = np.asarray(b, dtype=float)
    S0, D, alpha = float(theta[0]), float(theta[1]), float(theta[2])
    u = np.power(b * D, alpha)  # (b D)^alpha; (b=0) -> 0
    return S0 * np.exp(-u)


def jacobian(b, theta):
    """Closed-form dS/dtheta at each b. Returns array of shape (len(b), 3).

    With u = (b D)^alpha and S = S0 exp(-u):
        dS/dS0    =  S / S0
        dS/dD     = -S * alpha * u / D
        dS/dalpha = -S * u * ln(b D)
    At b = 0: u = 0 so dS/dD = 0 and dS/dalpha = 0 (only S0 is informed).
    """
    b = np.asarray(b, dtype=float)
    S0, D, alpha = float(theta[0]), float(theta[1]), float(theta[2])
    bD = b * D
    u = np.power(bD, alpha)
    S = S0 * np.exp(-u)

    dS0 = S / S0
    dD = -S * alpha * u / D
    # ln(bD) is -inf at b=0, but u=0 there so u*ln(bD) -> 0. Evaluate log only where bD>0
    # to avoid 0*(-inf)=nan; the derivative is exactly 0 at b=0 (only S0 is informed).
    safe_bD = np.where(bD > 0.0, bD, 1.0)
    ln_bD = np.log(safe_bD)
    dalpha = np.where(bD > 0.0, -S * u * ln_bD, 0.0)

    return np.stack([dS0, dD, dalpha], axis=-1)


def forward_joint(b, theta_joint):  # pragma: no cover - Phase 3 stub
    """Joint CTRW / fractional Bloch-Torrey two-exponent forward model (Phase 3 only).

    theta_joint = (S0, D, alpha, beta). For CP0 only the time-alpha lane is active and
    this reduces to the stretched-exponential. The space-fractional beta coupling that
    creates the (alpha, beta) degeneracy is deferred to Phase 3 (only built if Gate C
    stands), per the protocol.
    """
    raise NotImplementedError(
        "forward_joint (alpha, beta) is the Phase-3 object; CP0 uses the stretched-exp lead lane."
    )
