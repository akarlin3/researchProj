"""Fisher information and Cramer-Rao lower bounds for the fractional order under the
diffusion-MRI signal-decay forward model -- THIS IS LEVY'S NET-NEW CONTRIBUTION.

The upstream Ouroboros tooling has no Fisher/CRLB/identifiability layer (audited, GATE A);
this module builds it. The estimand is theta = (S0, D, alpha) estimated *jointly* from a
finite set of b-values under Rician noise.

Two Fisher matrices (parameter enters likelihood only through nu_i = S(b_i; theta)):

    FIM_G(theta) = (1/sigma^2) * sum_i g_i g_i^T                     (high-SNR Gaussian limit)
    FIM_R(theta) = sum_i I_R(nu_i, sigma) g_i g_i^T
                 = (1/sigma^2) * sum_i f(nu_i/sigma) g_i g_i^T       (honest finite-SNR Rician)

with g_i = dS(b_i)/dtheta (closed form, forward.jacobian) and f the Rician info factor
(noise.rician_info_factor). FIM_R -> FIM_G as SNR -> inf. CRLB(theta_k) = [FIM^{-1}]_kk;
the standard error of any unbiased estimator obeys SE(theta_k) >= sqrt(CRLB_k).

Identifiability diagnostics for the WALL:
  * crlb_alpha and the *relative* CRLB  cv_alpha = sqrt(CRLB_alpha)/alpha  (info-limit on alpha)
  * rho_alpha_D: the (alpha, D) correlation from the inverse FIM -> +-1 signals degeneracy
  * cond: FIM condition number -> blows up where the experiment cannot separate the parameters

LICENSE/SCOPE NOTE: a CRLB is an *identifiability / information* statement, never an
impossibility one. Every number here is scoped to its regime (forward model, SNR, b-design).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import forward, noise


@dataclass(frozen=True)
class CRLBResult:
    """CRLB diagnostics at one (truth, b-design, SNR) cell."""

    theta: np.ndarray          # (S0, D, alpha) ground truth
    snr: float
    sigma: float
    fim: np.ndarray            # 3x3 Fisher information matrix (Rician or Gaussian)
    cov: np.ndarray            # FIM^{-1} = CRLB covariance
    model: str                 # "rician" or "gaussian"

    @property
    def crlb(self) -> np.ndarray:
        return np.diag(self.cov)

    @property
    def crlb_alpha(self) -> float:
        return float(self.cov[forward.IDX["alpha"], forward.IDX["alpha"]])

    @property
    def crlb_D(self) -> float:
        return float(self.cov[forward.IDX["D"], forward.IDX["D"]])

    @property
    def se_alpha(self) -> float:
        return float(np.sqrt(max(self.crlb_alpha, 0.0)))

    @property
    def cv_alpha(self) -> float:
        """Relative CRLB on alpha: sqrt(CRLB_alpha)/alpha. The info-limit recovery metric."""
        return self.se_alpha / float(self.theta[forward.IDX["alpha"]])

    @property
    def rho_alpha_D(self) -> float:
        """(alpha, D) correlation implied by the inverse FIM; +-1 = degeneracy."""
        ia, iD = forward.IDX["alpha"], forward.IDX["D"]
        denom = np.sqrt(self.cov[ia, ia] * self.cov[iD, iD])
        if denom <= 0:
            return float("nan")
        return float(self.cov[ia, iD] / denom)

    @property
    def cond(self) -> float:
        """FIM condition number (2-norm). Large => parameters not separable by this design."""
        return float(np.linalg.cond(self.fim))


def fisher_matrix(b, theta, snr, model: str = "rician"):
    """Fisher information matrix for theta = (S0, D, alpha) at b-design ``b`` and ``snr``.

    model="gaussian": (1/sigma^2) J^T J  (high-SNR reference)
    model="rician":   J^T diag(I_R(nu_i, sigma)) J  (honest finite-SNR)
    """
    b = np.asarray(b, dtype=float)
    theta = np.asarray(theta, dtype=float)
    S0 = float(theta[0])
    sigma = noise.sigma_from_snr(S0, snr)

    J = forward.jacobian(b, theta)                 # (m, 3)
    nu = forward.signal(b, theta)                  # (m,)

    if model == "gaussian":
        weights = np.full(nu.shape, 1.0 / (sigma * sigma))
    elif model == "rician":
        a = nu / sigma                             # local SNR per b
        weights = noise.rician_info_factor(a) / (sigma * sigma)
    else:
        raise ValueError(f"unknown model {model!r} (use 'gaussian' or 'rician')")

    # FIM = sum_i w_i g_i g_i^T  =  J^T diag(w) J
    fim = J.T @ (weights[:, None] * J)
    return fim


def crlb(b, theta, snr, model: str = "rician") -> CRLBResult:
    """Compute the CRLB diagnostics at one cell. Inverts the FIM (pinv fallback if singular)."""
    b = np.asarray(b, dtype=float)
    theta = np.asarray(theta, dtype=float)
    sigma = noise.sigma_from_snr(float(theta[0]), snr)
    fim = fisher_matrix(b, theta, snr, model=model)
    try:
        cov = np.linalg.inv(fim)
    except np.linalg.LinAlgError:
        cov = np.linalg.pinv(fim)
    return CRLBResult(theta=theta, snr=float(snr), sigma=float(sigma), fim=fim, cov=cov, model=model)
