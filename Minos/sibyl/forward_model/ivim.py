import torch
from typing import Optional, Tuple

# Default dense synthetic 10-value b-scheme for Tier 1
DEFAULT_B_SCHEME = torch.tensor([0, 10, 20, 50, 100, 200, 400, 600, 800, 1000], dtype=torch.float32)

# ACRIN-6698 / I-SPY2 acquired b-scheme (s/mm^2). Verified against the TCIA data
# descriptor and Partridge 2018 / Newitt 2019: fat-suppressed SS-EPI at exactly
# these four b-values, 3-direction trace.
ACRIN_B_SCHEME = torch.tensor([0, 100, 600, 800], dtype=torch.float32)

# Indices of the four ACRIN b-values inside DEFAULT_B_SCHEME (the dense grid):
# 0->0, 100->4, 600->7, 800->8. The other six positions are imputed in Arm 2.
ACRIN_TO_DENSE_IDX = (0, 4, 7, 8)

def ivim_biexponential(
    bvals: torch.Tensor,
    f: torch.Tensor,
    D: torch.Tensor,
    Dstar: torch.Tensor,
    S0: Optional[torch.Tensor] = None
) -> torch.Tensor:
    """
    IVIM biexponential forward model.

    Parameters
    ----------
    bvals : torch.Tensor
        1D tensor of b-values (shape: [B])
    f : torch.Tensor
        Perfusion fraction (shape: [N, 1] or [N])
    D : torch.Tensor
        Tissue diffusion coefficient (shape: [N, 1] or [N])
    Dstar : torch.Tensor
        Pseudo-diffusion coefficient (shape: [N, 1] or [N])
    S0 : torch.Tensor, optional
        Signal at b=0. If None, S0=1.0 is used. (shape: [N, 1] or [N])

    Returns
    -------
    torch.Tensor
        Simulated signal (shape: [N, B])
    """
    if f.ndim == 1:
        f = f.unsqueeze(-1)
    if D.ndim == 1:
        D = D.unsqueeze(-1)
    if Dstar.ndim == 1:
        Dstar = Dstar.unsqueeze(-1)
    if S0 is not None and S0.ndim == 1:
        S0 = S0.unsqueeze(-1)
        
    bvals = bvals.view(1, -1)  # shape [1, B]
    
    signal = (f * torch.exp(-bvals * (D + Dstar))) + ((1 - f) * torch.exp(-bvals * D))
    
    if S0 is not None:
        signal = S0 * signal

    return signal


def adc_monoexp_fit(
    signal: torch.Tensor,
    bvals: torch.Tensor,
    eps: float = 1e-8,
) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Mono-exponential ADC by log-linear least squares over ALL supplied b-values.

    This is the canonical ACRIN-6698 ADC definition (Partridge 2018, Newitt 2019):
    a linear least-squares fit of log(signal) vs b across b = 0, 100, 600, 800,
    under S(b) = S0 * exp(-b * ADC). It is the project's *external* trust reference
    and is deliberately independent of uGUIDE.

    Parameters
    ----------
    signal : torch.Tensor
        Signal intensities (shape: [N, B] or [B]). Must be positive.
    bvals : torch.Tensor
        b-values (shape: [B]).
    eps : float
        Floor applied before taking the log, to keep noisy/near-zero signal finite.

    Returns
    -------
    S0 : torch.Tensor
        Fitted signal at b=0 (shape: [N]).
    ADC : torch.Tensor
        Apparent diffusion coefficient, same length-units as 1/b (shape: [N]).
    """
    if signal.ndim == 1:
        signal = signal.unsqueeze(0)
    signal = signal.to(torch.float64)
    bvals = bvals.to(torch.float64).view(-1)
    B = bvals.shape[0]
    assert signal.shape[1] == B, f"signal has {signal.shape[1]} b-values, bvals has {B}"

    # Design matrix A = [1, -b]; model log(S) = log(S0) - ADC * b.
    A = torch.stack([torch.ones(B, dtype=torch.float64), -bvals], dim=1)  # [B, 2]
    pinv = torch.linalg.pinv(A)  # [2, B], shared across all samples

    y = torch.log(signal.clamp_min(eps))  # [N, B]
    beta = y @ pinv.T  # [N, 2] -> columns [log S0, ADC]
    S0 = torch.exp(beta[:, 0])
    ADC = beta[:, 1]
    return S0.to(torch.float32), ADC.to(torch.float32)
