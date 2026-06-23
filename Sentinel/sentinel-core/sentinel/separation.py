"""CP0 — the separation gate that decides the paper.

Run the fractionated course, fire all three rules, and test whether the
regret-targeted stop halts at a **different session than both baselines**, in the
named regime: ACI still holds coverage (by widening) while decision value is dead.

Pre-registered REFUTE (any one kills the wedge -> RED):
  R-ACI   : ACI's recalibration reproduces the regret-stop's halt (it stops, or it
            fails to hold coverage at/ before t_regret) -> no decision/coverage gap.
  R-WATCH : the WATCH coverage-validity alarm fires at the same session as the
            regret-stop (gap CI includes 0) -> the regret-stop is just a coverage
            changepoint detector.
  R-REGIME: at t_regret, decision value is NOT yet dead, OR coverage is NOT still
            held -> the separation is not in the claimed regime.

The separation (t_watch - t_regret > 0, CI excludes 0, with ACI never stopping while
holding coverage, in the named regime) IS the paper.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .baselines import conformal_pvalues, run_aci, run_watch
from .config import SentinelConfig
from .course import build_course, decision_regret, Course
from .monitor import monitor_M, reference_density
from .seeding import make_rng
from .stopping import calibrate_regret_stop, regret_stop, regret_stop_seq, RegretStopCal, _cusum_path


@dataclass(frozen=True)
class SeparationResult:
    t_regret: int | None
    t_watch: int | None
    aci_stops: bool             # ACI ever stops? (should be False)
    aci_holds_at_regret: bool   # ACI still holds coverage at t_regret?
    regret_dead_at_regret: bool # decision value dead at t_regret?
    gap: float | None           # t_watch - t_regret (None if WATCH never fires)
    gap_ci: tuple[float, float] | None
    gap_boot: np.ndarray | None
    M_seq: np.ndarray
    cusum_b: float
    cusum_h: float
    regrets: np.ndarray
    aci_coverage: np.ndarray
    aci_width_mult: np.ndarray
    watch_logmart: np.ndarray
    watch_thr: float
    refute: dict

    @property
    def separated(self) -> bool:
        """Wedge survives: WATCH fires strictly after regret-stop with CI>min_gap,
        ACI never stops, and the named regime holds at t_regret."""
        if self.gap is None or self.gap_ci is None or self.t_regret is None:
            return False
        return (not self.aci_stops and self.aci_holds_at_regret
                and self.regret_dead_at_regret
                and self.gap_ci[0] > 0.0)


def _watch_stop_from_scores(score_mat: np.ndarray, cfg: SentinelConfig,
                            rng: np.random.Generator, eps: np.ndarray,
                            log_eps: np.ndarray) -> int | None:
    """WATCH alarm session from a (n_sessions x n_voxels) score matrix (bootstrap-fast)."""
    cal = score_mat[0]
    cum = np.zeros(eps.size)
    log_thr = float(np.log(1.0 / cfg.watch_delta))
    nb = min(cfg.watch_bins_per_session, score_mat.shape[1])
    for k in range(1, score_mat.shape[0]):
        idx = rng.choice(score_mat.shape[1], size=nb, replace=False)
        p = np.clip(conformal_pvalues(cal, score_mat[k][idx], rng), 1e-12, 1.0)
        cum = cum + (log_eps[:, None] + np.outer(eps - 1.0, np.log(p))).sum(axis=1)
        m = float(np.max(cum)); lm = m + float(np.log(np.sum(np.exp(cum - m)))) - np.log(eps.size)
        if lm >= log_thr:
            return k
    return None


def _regret_stop_from_mu(mu_mat: np.ndarray, p_ref: np.ndarray, cal: RegretStopCal,
                         cfg: SentinelConfig) -> int | None:
    M = np.array([monitor_M(mu_mat[k], p_ref, cfg) for k in range(mu_mat.shape[0])])
    cr = np.where(_cusum_path(M, cal.b) >= cal.h)[0]
    return int(cr[0]) if cr.size else None


def run_separation(cfg: SentinelConfig, *, f_true: np.ndarray | None = None,
                   matrix_resid_sd: float | None = None,
                   boot_n: int | None = None) -> SeparationResult:
    course = build_course(cfg, f_true=f_true, matrix_resid_sd=matrix_resid_sd)
    f_true = course.sessions[0].f_true

    # --- point estimates -------------------------------------------------------
    p_ref = reference_density(course.sessions[0].mu, cfg)
    cal = calibrate_regret_stop(f_true, p_ref, cfg)
    t_regret, _, M_seq = regret_stop(course, cfg, cal=cal)
    aci = run_aci(course, cfg)
    watch = run_watch(course, cfg)
    t_watch = watch.stop_session

    aci_stops = aci.stop_session is not None
    regrets = course.regrets

    # no-drift regret null: same patient, fresh noise, drift switched off, over the
    # course length -> the irreducible measurement-noise regret floor + its spread.
    null_course = build_course(cfg.replace(drift_rate=0.0), f_true=f_true,
                               matrix_resid_sd=course.matrix_resid_sd)
    null_regret = null_course.regrets
    null_mean, null_sd = float(np.mean(null_regret)), float(np.std(null_regret))
    regret_floor = null_mean

    def regime_at(t: int | None) -> tuple[bool, bool]:
        if t is None:
            return False, False
        holds = bool(aci.holds_coverage[t])
        left_band = regrets[t] >= null_mean + cfg.regret_null_z * max(null_sd, 1e-12)
        practical = regrets[t] >= (1.0 + cfg.regret_rise_min) * max(regret_floor, 1e-12)
        return holds, bool(left_band and practical)

    aci_holds_at_regret, regret_dead_at_regret = regime_at(t_regret)
    gap = None if (t_watch is None or t_regret is None) else float(t_watch - t_regret)

    # --- voxel bootstrap on the stop-time gap ---------------------------------
    gap_boot = None
    gap_ci = None
    if t_regret is not None:
        boot_n = cfg.boot_n if boot_n is None else boot_n
        mu_mat = np.stack([s.mu for s in course.sessions])      # (K, V)
        score_mat = np.stack([s.score for s in course.sessions])
        rng = make_rng(cfg.seed + 4242)
        eps = np.linspace(1e-3, 1.0 - 1e-3, 199)
        log_eps = np.log(eps)
        V = f_true.size
        gaps = []
        for _ in range(boot_n):
            bidx = rng.integers(0, V, size=V)                    # resample patient voxels
            mub = mu_mat[:, bidx]
            scb = score_mat[:, bidx]
            p_ref_b = reference_density(mub[0], cfg)
            tr = _regret_stop_from_mu(mub, p_ref_b, cal, cfg)    # fixed rule (b, h)
            tw = _watch_stop_from_scores(scb, cfg, rng, eps, log_eps)
            if tr is not None and tw is not None:
                gaps.append(tw - tr)
        gap_boot = np.array(gaps, dtype=float)
        if gap_boot.size:
            lo = float(np.quantile(gap_boot, (1 - cfg.boot_ci) / 2))
            hi = float(np.quantile(gap_boot, 1 - (1 - cfg.boot_ci) / 2))
            gap_ci = (lo, hi)

    refute = {
        "R-ACI (ACI stops or fails coverage at t_regret)": bool(
            aci_stops or (t_regret is not None and not aci_holds_at_regret)),
        "R-WATCH (gap CI includes 0)": bool(
            gap_ci is None or gap_ci[0] <= 0.0),
        "R-REGIME (decision value not dead at t_regret)": bool(
            t_regret is not None and not regret_dead_at_regret),
    }

    return SeparationResult(
        t_regret=t_regret, t_watch=t_watch, aci_stops=aci_stops,
        aci_holds_at_regret=aci_holds_at_regret,
        regret_dead_at_regret=regret_dead_at_regret,
        gap=gap, gap_ci=gap_ci, gap_boot=gap_boot, M_seq=M_seq,
        cusum_b=cal.b, cusum_h=cal.h,
        regrets=regrets, aci_coverage=aci.coverage, aci_width_mult=aci.width_mult,
        watch_logmart=watch.log_mart, watch_thr=watch.threshold, refute=refute)
