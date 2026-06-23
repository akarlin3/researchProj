"""Frozen configuration for projSentinel.

One dataclass holds (i) the course / fractionated-session structure, (ii) the
decision problem (perfusion threshold + cost asymmetry), (iii) the
**fractionated-session hidden-drift model** (the Phase-1 enabler's scientific
object), (iv) the regret-targeted monitor's stakes kernel, (v) the two real
baselines' knobs (ACI/PID step sizes; WATCH false-alarm budget), and (vi) the
pre-registered gate thresholds for the CP0 separation test.

All thresholds that decide the verdict live here so the gate is auditable.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

from .seeding import GLOBAL_SEED


@dataclass(frozen=True)
class SentinelConfig:
    # ---- course / fractionated-session structure -----------------------------
    n_sessions: int = 20          # length of the RT course (sessions monitored)
    n_voxels: int = 4000          # reported points per session (decision units)
    matrix_iter: int = 0          # which Matrix loop iteration to read as the "scan"

    # ---- decision problem (perfusion treat/spare at a threshold) -------------
    f_treat: float = 0.16         # treat/spare perfusion threshold (Matrix-aligned)
    s_f: float = 0.04             # reported sd of the perfusion estimate
    k_under: float = 4.0          # cost of UNDER-treating (miss a low-perfusion voxel)
    k_over: float = 1.0           # cost of OVER-treating (asymmetry: under > over)

    # ---- fractionated-session hidden-drift model (the enabler) ---------------
    # A measurement drift on the reported perfusion, accumulating session-to-session
    # and CONCENTRATED near the decision band (the boosted, marginal voxels): the
    # cumulative-RT signal perturbation is co-located with the voxels actually
    # receiving dose. drift_rate is per-session bias amplitude (in perfusion units);
    # drift_band is the kernel half-width in reported-sd units around the threshold.
    drift_rate: float = 0.010     # per-session near-threshold reported-bias amplitude
    drift_band: float = 0.80      # Gaussian kernel width (reported-sd) around threshold
    drift_sign: float = +1.0      # +1: inflate reported perfusion of marginal voxels

    # ---- regret-targeted monitor M (ported from Minos) -----------------------
    z_range: float = 5.0          # reported-coordinate histogram half-range (sd)
    z_bins: int = 50
    mon_alpha: float = 0.05       # false-alarm level for m* (no-drift null quantile)
    mon_null_seeds: int = 80      # null draws to calibrate m*

    # ---- baseline 1: ACI / conformal-PID (recalibrate forever) ---------------
    aci_alpha: float = 0.10       # target miscoverage (90% intervals)
    aci_gamma: float = 0.05       # ACI step size (P term)
    pid_k_i: float = 0.02         # PID integral gain (0 => pure ACI)
    pid_sat: float = 50.0         # integral saturation bound
    aci_width_cap_mult: float = 6.0  # "can no longer hold coverage" only if width must
                                     # exceed this multiple of the nominal width

    # ---- baseline 2: WATCH-style conformal test martingale -------------------
    watch_delta: float = 0.05     # false-alarm budget; Ville threshold c = 1/delta
    watch_bins_per_session: int = 200  # conformal p-values drawn per session

    # ---- CP0 separation gate thresholds (pre-registered) ---------------------
    min_stop_gap: float = 1.0     # |t_watch - t_regret| must exceed this (sessions)
    boot_n: int = 1000            # bootstrap resamples for the stop-time-gap CI
    boot_ci: float = 0.95
    # "decision value can no longer be held": regret at t_regret must have left the
    # no-drift null band (statistical) AND risen a practically-meaningful amount.
    regret_null_z: float = 3.0    # >= null_mean + z*null_sd (statistical significance)
    regret_rise_min: float = 0.10 # AND >= (1 + this) x no-drift regret floor (practical)
    coverage_ok_tol: float = 0.03 # "ACI still holds coverage" = |cov - target| <= tol

    # ---- reproducibility -----------------------------------------------------
    seed: int = GLOBAL_SEED

    def replace(self, **kw) -> "SentinelConfig":
        return replace(self, **kw)


DEFAULT = SentinelConfig()
