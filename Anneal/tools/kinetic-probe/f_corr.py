"""CP2 plug-in point — finite-N kinetic-theory correction to the collective flow.

THIS FILE IS THE HUMAN ANALYTICAL CHECKPOINT (CP2).

The body of `f_corr` below is intentionally the ZERO correction. It is the
placeholder for the human-supplied system-size-expansion correction
(Buice-Chow / Hildebrand-Buice-Chow class) to the deterministic reduced flow
of tools/reduced-ode/reduced_core.py (rhs_3d). The verification harness
(harness.py / run_probe.py) integrates rhs_3d PLUS whatever this function
returns and scores the result against four pre-committed conditions, untuned.

The harness authors (CC) must NOT fill this in. The expansion, the closure,
the boundary treatment near the homoclinic ghost, and every coefficient are
the human's to supply. If the derivation does not close, or produces no
N-independent term, that is itself the result and is recorded as a clean
negative.

Conventions the human body must match (see HARNESS.md for the full list):
  - State is (rho1, rho2, psi) in the conventions of reduced_core.rhs_3d:
    mu = (1+A)/2, nu = (1-A)/2, alpha = pi/2 - beta, omega = 0; at the
    operating corner A = 0.5, beta = 0.05.
  - 1:1 clock: one model time unit = one second; rates in s^-1.
  - drift_vec is an ADDITIVE deterministic correction to rhs_3d (s^-1).
  - diff_matrix is the noise-AMPLITUDE matrix B in
        dX = [rhs_3d(X) + drift_vec(X, N)] dt + B(X, N) dW,
    with W a 3-vector of independent standard Wiener processes, so the
    Fokker-Planck diffusion matrix is D = B B^T. If the derivation yields D,
    supply B = any matrix square root of D (e.g. Cholesky).
  - N-dependence lives entirely inside this function; the harness sweeps
    N in {8, 16, 32, 64} and never rescales the returned values.
  - Coefficients must be theory-fixed. No free parameters tuned toward the
    3.2x target. If a coefficient is genuinely undetermined by the
    derivation, state its range in F_CORR_META["undetermined_coefficients"]
    and the harness will report the range of outcomes, not a chosen value.
"""
import numpy as np

F_CORR_META = {
    "name": "zero",
    "description": "placeholder: no correction (deterministic reduced flow)",
    "supplied_by": None,  # set to "human" (with a derivation reference) at CP2
    "derivation_reference": None,
    "theory_fixed_coefficients": None,
    "undetermined_coefficients": None,  # {"name": [lo, hi]} if any
    # CP2 basis decision (human, 2026-06-10): the derivation will follow the
    # Tyulkina-Goldobin-Klimenko-Pikovsky circular-cumulant route
    # (PRL 120, 264101 (2018)). Basis chosen != body supplied: the
    # system-specific closure choices and coefficients are still pending
    # human input. See CP2_TGKP_WORKSHEET.md for the open slots.
    "basis_chosen": "TGKP circular-cumulant (human decision 2026-06-10)",
}

_ZERO3 = np.zeros(3)
_ZERO33 = np.zeros((3, 3))


def f_corr(rho1, rho2, psi, N):
    """Finite-N kinetic-theory correction to the reduced collective flow.

    Parameters
    ----------
    rho1, rho2, psi : float
        Current collective state (order parameters and phase difference).
    N : int
        Oscillators per population.

    Returns
    -------
    drift_vec : ndarray, shape (3,)
        Additional drift on (rho1, rho2, psi), units s^-1.
    diff_matrix : ndarray, shape (3, 3)
        Noise-amplitude matrix B (see module docstring); zero matrix means
        no stochastic term. Returned arrays must not be mutated by callers.
    """
    # >>> HUMAN ANALYTICAL INPUT GOES HERE (CP2). Do not auto-fill. <<<
    return _ZERO3, _ZERO33
