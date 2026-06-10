"""AGENT-DERIVED, EXPLORATORY — TGKP circular-cumulant correction (slaved kappa).

STATUS / PROVENANCE (read before trusting anything in this file):
The human (Avery) authorized agent completion of the TGKP route on
2026-06-10 ("Do it as you think is best, but feel free to complete it").
Under the original v6 protocol this is a DOWNGRADE: the closure below is
agent-assembled, so its score CANNOT carry the evidentiary weight of the
pre-committed bet (the agent knows the 3.2x target). Outputs are therefore
labelled exploratory; f_corr.py remains zero, F_CORR_META['supplied_by']
remains None, and run_probe.py still gates on a human-supplied derivation.
This file answers a cheaper question: is the TGKP route FEASIBLE at
data-anchored magnitudes?

Derivation basis (published): Tyulkina, Goldobin, Klimenko & Pikovsky,
PRL 120, 264101 (2018), Eq. (13), two-cumulant truncation, identical
oscillators (gamma = 0), rotating frame (Omega_0 = 0):

    Z'     = h - conj(h) Z^2 - sigma^2 Z - conj(h) kappa
    kappa' = -4 conj(h) Z kappa - sigma^2 (4 kappa + 2 Z^2)

AGENT CHOICES (worksheet slots), each declared:
  Slot 2 (forcing, EXACT + verified): matching Z_s' = h_s - conj(h_s) z_s^2
    to reduced_core.rhs_3d gives h_1 = (1/2) e^{-i alpha} (mu z_1 + nu z_2),
    h_2 = (1/2) e^{-i alpha} (mu z_2 + nu z_1), z_s = rho_s e^{i phi_s},
    psi = phi_1 - phi_2. verify_oa_match() asserts machine-precision
    agreement with rhs_3d; this part involves no judgment.
  Slot 3 (cumulant treatment): kappa adiabatically SLAVED (kappa' = 0):
    kappa_s = -sigma^2 z_s^2 / (2 (conj(h_s) z_s + sigma^2)).
    Justification: TGKP's own scaling kappa ~ sigma^2; truncation error
    O(sigma^4). Dynamical-kappa augmentation deferred (would need an
    augmented-state driver); this is a known judgment call near the ghost.
  Slot 4 (character): drift-only (B = 0) primary — the TGKP correction is
    deterministic at the ensemble level, and the stochastic 1/sqrt(N) CLT
    part was separately excluded by Appendix B. Variant V5 adds the
    Appendix-B-calibrated additive noise (c = 0.05/sqrt(N)) on top.
  Slot 1 (sigma^2_eff — the LOAD-BEARING, genuinely underived slot):
    The measured per-population order-parameter fluctuation satisfies
    sigma_HF * sqrt(N) ~ 0.040-0.054 ~ const (results_mech.json,
    physical_estimate), supporting an N-FREE effective per-oscillator
    phase-noise intensity. Converting the measured AMPLITUDE to a RATE
    requires a relaxation rate gamma_R (OU inversion var ~ D/(2 gamma_R)
    => per-oscillator sigma^2 ~ 2 gamma_R (sigma_HF sqrt(N))^2). gamma_R
    is NOT derived here; it is declared as the undetermined coefficient
    with bracket gamma_R in [0.1, 1] s^-1 (relaxation times 1-10 s, i.e.
    sub-breath-period), scanned and never selected.

PRE-DECLARED VARIANTS (frozen in this file BEFORE any scoring run; the
no-tuning rule applies — no variant is added, removed, or re-parameterized
after results are seen):
  V1 sigma^2 = 2*0.3*0.047^2 = 1.325e-3 s^-1 (pooled c, gamma_R = 0.3), B=0
  V2 sigma^2 = 2*0.1*0.047^2 = 4.418e-4 s^-1 (bracket low), B=0
  V3 sigma^2 = 2*1.0*0.047^2 = 4.418e-3 s^-1 (bracket high), B=0
  V4 per-N: sigma^2(N) = 2*0.3*(sigma_HF(N)*sqrt(N))^2 (per-N measured
     c(N) = 0.0402, 0.0454, 0.0473, 0.0544 -> carries the data's mild
     N-trend honestly), B=0
  V5 = V1 drift + Appendix-B additive noise B = (0.05/sqrt(N)) I

Direction-of-effect note (pre-run): the -sigma^2 rho drag pulls the order
parameters away from the capture boundary, so prolongation (not shortening)
is the expected sign; whether the magnitude reaches ~3.2x is what the scan
measures.
"""
import cmath
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "tools/reduced-ode"))

A = 0.5
BETA = 0.05
MU = (1.0 + A) / 2.0
NU = (1.0 - A) / 2.0
ALPHA = np.pi / 2.0 - BETA
_EIA = cmath.exp(-1j * ALPHA)
_EYE3 = np.eye(3)
_ZERO33 = np.zeros((3, 3))

# Slot-1 anchors (measured, results_mech.json physical_estimate)
C_POOLED = 0.047                       # pooled sigma_HF*sqrt(N)
C_PER_N = {8: 0.0402, 16: 0.0454, 32: 0.0473, 64: 0.0544}
GAMMA_R_CENTRAL, GAMMA_R_LO, GAMMA_R_HI = 0.3, 0.1, 1.0

SIGMA2_V1 = 2.0 * GAMMA_R_CENTRAL * C_POOLED**2
SIGMA2_V2 = 2.0 * GAMMA_R_LO * C_POOLED**2
SIGMA2_V3 = 2.0 * GAMMA_R_HI * C_POOLED**2
SIGMA2_V4 = {N: 2.0 * GAMMA_R_CENTRAL * c**2 for N, c in C_PER_N.items()}
NOISE_C_V5 = 0.05                      # Appendix-B calibrated amplitude


def _u(rho1, rho2, psi):
    """Gauge-invariant forcings u_s = h_s e^{-i phi_s} (Slot 2, exact)."""
    u1 = 0.5 * _EIA * (MU * rho1 + NU * rho2 * cmath.exp(-1j * psi))
    u2 = 0.5 * _EIA * (MU * rho2 + NU * rho1 * cmath.exp(1j * psi))
    return u1, u2


def tgkp_drift(rho1, rho2, psi, sigma2):
    """Slaved-kappa TGKP drift correction on (rho1, rho2, psi).

    Gauge-invariant transcription: with u_s = h_s e^{-i phi_s} and
    kt_s = kappa_s e^{-2 i phi_s} (slaved):
        kt_s    = -sigma^2 rho_s^2 / (2 (conj(u_s) rho_s + sigma^2))
        Delta_s = -sigma^2 rho_s - conj(u_s) kt_s
        d rho_s = Re(Delta_s);  d phi_s = Im(Delta_s)/rho_s;
        d psi   = d phi_1 - d phi_2.
    Im(Delta_s)/rho_s -> 0 as rho_s -> 0 (Delta_s imaginary part is
    O(rho^2)); the epsilon guard only avoids 0/0 at exactly rho = 0.
    """
    u1, u2 = _u(rho1, rho2, psi)
    kt1 = -sigma2 * rho1 * rho1 / (2.0 * (u1.conjugate() * rho1 + sigma2))
    kt2 = -sigma2 * rho2 * rho2 / (2.0 * (u2.conjugate() * rho2 + sigma2))
    d1 = -sigma2 * rho1 - u1.conjugate() * kt1
    d2 = -sigma2 * rho2 - u2.conjugate() * kt2
    r1 = rho1 if rho1 > 1e-12 else 1e-12
    r2 = rho2 if rho2 > 1e-12 else 1e-12
    return np.array([d1.real, d2.real, d1.imag / r1 - d2.imag / r2])


# ---- picklable f_corr factories (module-level, used via functools.partial) --
def fcorr_tgkp(rho1, rho2, psi, N, sigma2):
    return tgkp_drift(rho1, rho2, psi, sigma2), _ZERO33


def fcorr_tgkp_perN(rho1, rho2, psi, N, sigma2_per_N):
    return tgkp_drift(rho1, rho2, psi, sigma2_per_N[N]), _ZERO33


_B_CACHE = {}


def fcorr_tgkp_noise(rho1, rho2, psi, N, sigma2, noise_c):
    B = _B_CACHE.get(N)
    if B is None:
        B = _B_CACHE.setdefault(N, (noise_c / np.sqrt(N)) * _EYE3)
    return tgkp_drift(rho1, rho2, psi, sigma2), B


# ---- self-verification gates (run before any scoring) -----------------------
def verify_oa_match(n_pts=2000, seed=0):
    """The deterministic OA form built from u_s must reproduce rhs_3d to
    machine precision (Slot-2 exactness gate):
      rho1' = (1-rho1^2) Re(u1); rho2' = (1-rho2^2) Re(u2);
      psi'  = (1+rho1^2) Im(u1)/rho1 - (1+rho2^2) Im(u2)/rho2.
    """
    import reduced_core as rc
    p = rc.Params(A=A, beta=BETA)
    rng = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(n_pts):
        r1, r2 = rng.uniform(0.05, 0.999, 2)
        psi = rng.uniform(-np.pi, np.pi)
        u1, u2 = _u(r1, r2, psi)
        mine = np.array([
            (1 - r1 * r1) * u1.real,
            (1 - r2 * r2) * u2.real,
            (1 + r1 * r1) * u1.imag / r1 - (1 + r2 * r2) * u2.imag / r2,
        ])
        ref = rc.rhs_3d([r1, r2, psi], p)
        worst = max(worst, float(np.max(np.abs(mine - ref))))
    assert worst < 1e-12, f"OA forcing does not match rhs_3d: {worst}"
    return worst


def verify_zero_limit(n_pts=200, seed=1):
    """sigma^2 -> 0 must kill the correction (and scale ~ sigma^2)."""
    rng = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(n_pts):
        r1, r2 = rng.uniform(0.05, 0.999, 2)
        psi = rng.uniform(-np.pi, np.pi)
        d_small = tgkp_drift(r1, r2, psi, 1e-12)
        worst = max(worst, float(np.max(np.abs(d_small))))
    assert worst < 1e-10, f"zero-noise limit violated: {worst}"
    return worst


VARIANTS = [
    ("V1_central", dict(kind="plain", sigma2=SIGMA2_V1)),
    ("V2_gammaR_low", dict(kind="plain", sigma2=SIGMA2_V2)),
    ("V3_gammaR_high", dict(kind="plain", sigma2=SIGMA2_V3)),
    ("V4_perN_measured", dict(kind="perN", sigma2_per_N=SIGMA2_V4)),
    ("V5_central_plus_CLTnoise",
     dict(kind="noise", sigma2=SIGMA2_V1, noise_c=NOISE_C_V5)),
]
