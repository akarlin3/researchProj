"""The fractionated-session axis — projSentinel's Phase-1 enabler.

Matrix's drift is *within one run* (dose-response, no session axis). The enabler
wraps Matrix to build a **course**: the same patient is monitored across
``n_sessions`` fractions, and a hidden measurement drift **accumulates
session-to-session**, concentrated in the decision band (the marginal, boosted
voxels). Matrix supplies the patient (ground-truth perfusion field + a
Matrix-grounded reporting scale), read **read-only**; the session axis,
measurement noise, and accumulating drift are projSentinel's own objects.

Per session ``k`` the wrapper exposes everything the three competing stopping
rules consume: reported points ``mu_k`` (for the regret-targeted monitor),
per-session miscoverage of the 1-alpha interval (for ACI/PID), conformal
nonconformity scores (for WATCH), and the oracle decision regret (validation).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import SentinelConfig
from .matrix_bridge import load_matrix
from .seeding import make_rng


# --------------------------------------------------------------------------------------
# patient: ground-truth perfusion field, pooled from the Matrix twin (read-only)
# --------------------------------------------------------------------------------------
def matrix_patient(cfg: SentinelConfig) -> tuple[np.ndarray, float]:
    """Pool Matrix twin runs into a patient: per-voxel true perfusion ``f_true`` of
    length ``cfg.n_voxels`` plus the Matrix-empirical reported-perfusion residual sd.

    Read-only: runs Matrix's closed loop (byte-identity asserted in ``load_matrix``)
    and reads ``state.truth['f']`` (ground truth) and ``state.mu['f']`` (reported) at
    iteration ``cfg.matrix_iter``. Never edits Matrix.
    """
    matrix = load_matrix()
    run_loop = matrix.run_loop
    MatrixConfig = matrix.MatrixConfig
    Interfaces = matrix.Interfaces

    f_true_parts: list[np.ndarray] = []
    resid_parts: list[np.ndarray] = []
    sub = 0
    while sum(len(p) for p in f_true_parts) < cfg.n_voxels:
        mcfg = MatrixConfig(seed=cfg.seed + 7919 * sub)
        _, states = run_loop(mcfg, Interfaces.placeholders())
        st = states[cfg.matrix_iter]
        f_true_parts.append(np.asarray(st.truth["f"], dtype=float))
        resid_parts.append(np.asarray(st.mu["f"], dtype=float) - np.asarray(st.truth["f"], dtype=float))
        sub += 1

    f_true = np.concatenate(f_true_parts)[: cfg.n_voxels]
    matrix_resid_sd = float(np.std(np.concatenate(resid_parts)))
    return f_true, matrix_resid_sd


# --------------------------------------------------------------------------------------
# the accumulating, threshold-concentrated hidden drift (the enabler's mechanism)
# --------------------------------------------------------------------------------------
def drift_bias(z_true: np.ndarray, session: int, cfg: SentinelConfig) -> np.ndarray:
    """Accumulated reported-perfusion bias at session ``k`` for voxels at true
    distance-to-threshold ``z_true`` (in reported-sd units).

    ``b_k(z) = sign * drift_rate * k * exp(-0.5 (z / drift_band)^2)`` — grows linearly
    with session and is Gaussian-localised to the decision band. Far-from-threshold
    voxels (|z| >> band) are essentially untouched: that co-location of the drift with
    the decision is what splits regret from coverage.
    """
    kernel = np.exp(-0.5 * (z_true / cfg.drift_band) ** 2)
    return cfg.drift_sign * cfg.drift_rate * float(session) * kernel


@dataclass(frozen=True)
class Session:
    """Per-session observables consumed by the competing stopping rules."""

    k: int
    mu: np.ndarray            # reported perfusion points (drifted), length n_voxels
    f_true: np.ndarray        # ground-truth perfusion (fixed patient)
    score: np.ndarray         # conformal nonconformity scores |mu - f_true|
    miscover: float           # empirical miscoverage of the session-0 (1-alpha) interval
    regret: float             # oracle decision regret (cost-weighted), validation only


@dataclass(frozen=True)
class Course:
    sessions: list[Session]
    f_treat: float
    s_f: float
    matrix_resid_sd: float
    cfg: SentinelConfig

    @property
    def regrets(self) -> np.ndarray:
        return np.array([s.regret for s in self.sessions])


# --------------------------------------------------------------------------------------
# decision regret (oracle, validation only)
# --------------------------------------------------------------------------------------
def _treat(mu: np.ndarray, f_treat: float) -> np.ndarray:
    """TREAT (1) low-perfusion voxels (reported below threshold), else SPARE (0)."""
    return (mu < f_treat).astype(int)


def decision_regret(mu: np.ndarray, f_true: np.ndarray, cfg: SentinelConfig) -> float:
    """Mean cost-weighted regret of acting on reported ``mu`` vs the oracle on truth.

    Oracle treats iff ``f_true < f_treat``. Under-treatment (oracle treats, we spare)
    costs ``k_under``; over-treatment (oracle spares, we treat) costs ``k_over``. The
    asymmetry ``k_under > k_over`` is the decision-value stakes the monitor weights.
    """
    a = _treat(mu, cfg.f_treat)
    oracle = _treat(f_true, cfg.f_treat)
    under = (oracle == 1) & (a == 0)
    over = (oracle == 0) & (a == 1)
    return float(cfg.k_under * np.mean(under) + cfg.k_over * np.mean(over))


# --------------------------------------------------------------------------------------
# build the course
# --------------------------------------------------------------------------------------
def build_course(cfg: SentinelConfig, *, f_true: np.ndarray | None = None,
                 matrix_resid_sd: float | None = None) -> Course:
    """Run the fractionated-session axis and return per-session observables.

    The patient (``f_true``) comes from the Matrix twin unless injected (tests inject a
    fixed field to stay hermetic). Each session re-measures the patient with fresh
    noise and the accumulated near-threshold drift.
    """
    if f_true is None:
        f_true, matrix_resid_sd = matrix_patient(cfg)
    if matrix_resid_sd is None:
        matrix_resid_sd = cfg.s_f

    rng = make_rng(cfg.seed + 31)
    f_treat, s_f = cfg.f_treat, cfg.s_f
    z_true = (f_true - f_treat) / s_f

    # session-0 reference interval half-width (the calibration the baselines start from)
    half = float(np.quantile(np.abs(rng.normal(0.0, s_f, size=20000)), 1.0 - cfg.aci_alpha))

    sessions: list[Session] = []
    for k in range(cfg.n_sessions):
        noise = rng.normal(0.0, s_f, size=f_true.shape)
        mu = f_true + noise + drift_bias(z_true, k, cfg)
        score = np.abs(mu - f_true)
        miscover = float(np.mean(np.abs(mu - f_true) > half))
        regret = decision_regret(mu, f_true, cfg)
        sessions.append(Session(k=k, mu=mu, f_true=f_true, score=score,
                                miscover=miscover, regret=regret))

    return Course(sessions=sessions, f_treat=f_treat, s_f=s_f,
                  matrix_resid_sd=matrix_resid_sd, cfg=cfg)
