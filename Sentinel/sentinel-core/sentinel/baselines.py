"""The two real baselines the wedge must out-separate.

(1) **ACI / conformal-PID** (Gibbs & Candes, NeurIPS 2021, arXiv:2106.00170;
    Angelopoulos et al., NeurIPS 2023, arXiv:2307.16895): online recalibration of
    the interval level. ``alpha_{t+1} = alpha_t + gamma (alpha - err_t)`` (PID adds an
    integral term on the running coverage error). It **recalibrates forever — it
    never stops**; it "holds coverage by widening". We record, per session, the
    marginal coverage it achieves and the width multiple it needs.

(2) **WATCH-style conformal test martingale** (Prinster et al., ICML 2025,
    arXiv:2505.04608): stream conformal p-values of each session's nonconformity
    scores against the exchangeable session-0 reference; bet with a mixture
    (power) martingale; **alarm at the Ville stopping time** ``M_t >= 1/delta``.
    We use the *unweighted* martingale — the earliest-firing, strongest form — so a
    separation in favour of the regret-stop is conservative (the paper's covariate
    weighting only delays the alarm further).

Both are implemented from scratch (no faithful in-repo prior art existed; Gate A).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import SentinelConfig
from .course import Course
from .seeding import make_rng


# ======================================================================================
# (1) ACI / conformal-PID — recalibrate forever
# ======================================================================================
@dataclass(frozen=True)
class ACITrace:
    alpha: np.ndarray          # adaptive level alpha_t per session
    half_width: np.ndarray     # interval half-width per session
    coverage: np.ndarray       # marginal coverage achieved per session
    width_mult: np.ndarray     # half_width / session-0 half_width
    holds_coverage: np.ndarray # bool: |coverage - target| <= tol AND within width cap
    stop_session: int | None   # ACI never stops -> None


def run_aci(course: Course, cfg: SentinelConfig, *, use_pid: bool = True) -> ACITrace:
    """Batched ACI/PID. Half-width at session t is the ``(1 - alpha_t)`` quantile of the
    frozen session-0 calibration scores; ``err_t`` is the session-t empirical
    miscoverage at that width; ``alpha`` recalibrates toward the target.
    """
    cal = course.sessions[0].score                      # exchangeable calibration scores
    target = cfg.aci_alpha
    w0 = float(np.quantile(cal, 1.0 - target))

    alpha_t = target
    err_sum = 0.0
    alphas, halfs, covs, mults, holds = [], [], [], [], []
    for s in course.sessions:
        a = min(max(alpha_t, 1e-4), 0.999)
        w = float(np.quantile(cal, 1.0 - a))            # widen as alpha_t drops
        err_t = float(np.mean(s.score > w))             # batched miscoverage
        cov = 1.0 - err_t
        within_cap = w <= cfg.aci_width_cap_mult * w0
        holds_t = (abs(cov - (1.0 - target)) <= cfg.coverage_ok_tol) and within_cap

        alphas.append(a); halfs.append(w); covs.append(cov)
        mults.append(w / w0); holds.append(holds_t)

        # PID update on the coverage error g_t = err_t - target
        g = err_t - target
        err_sum += g
        integral = cfg.pid_k_i * float(np.clip(err_sum, -cfg.pid_sat, cfg.pid_sat)) if use_pid else 0.0
        alpha_t = alpha_t + cfg.aci_gamma * (target - err_t) - integral

    return ACITrace(alpha=np.array(alphas), half_width=np.array(halfs),
                    coverage=np.array(covs), width_mult=np.array(mults),
                    holds_coverage=np.array(holds), stop_session=None)


# ======================================================================================
# (2) WATCH-style conformal test martingale
# ======================================================================================
def conformal_pvalues(cal_scores: np.ndarray, test_scores: np.ndarray,
                      rng: np.random.Generator) -> np.ndarray:
    """Smoothed conformal p-values of ``test_scores`` against ``cal_scores``.

    ``p = (#{cal > v} + u #{cal == v}) / (n_cal + 1)``, ``u ~ Unif[0,1]``. Under
    exchangeability these are iid Uniform[0,1] — the null WATCH monitors. Large
    nonconformity (miscoverage) yields small p (mass piles near 0) and grows the bet.
    """
    cal = np.sort(np.asarray(cal_scores, float))
    n = cal.size
    gt = n - np.searchsorted(cal, test_scores, side="right")   # #{cal > v}
    eq = (np.searchsorted(cal, test_scores, side="right")
          - np.searchsorted(cal, test_scores, side="left"))     # #{cal == v}
    u = rng.uniform(size=test_scores.shape)
    return (gt + u * eq) / (n + 1.0)


@dataclass(frozen=True)
class WatchTrace:
    log_mart: np.ndarray       # log martingale after each session
    threshold: float           # log(1/delta) Ville threshold
    stop_session: int | None   # first session crossing the threshold (the alarm)


def run_watch(course: Course, cfg: SentinelConfig, *, eps_grid: int = 199) -> WatchTrace:
    """Mixture (power) conformal test martingale; Ville alarm at ``M_t >= 1/delta``.

    ``M_t = integral_0^1 prod_i eps p_i^{eps-1} d eps`` approximated on an ``eps`` grid,
    tracked in log space (cumulative log-bets per grid point, combined by logsumexp).
    Fires when ``log M_t >= log(1/delta)``.
    """
    rng = make_rng(cfg.seed + 9001)
    cal = course.sessions[0].score
    eps = np.linspace(1e-3, 1.0 - 1e-3, eps_grid)
    log_eps = np.log(eps)
    cum = np.zeros(eps_grid)                       # per-grid cumulative sum of log-bets
    log_thr = float(np.log(1.0 / cfg.watch_delta))

    log_mart, stop = [], None
    for s in course.sessions[1:]:                  # session 0 is the calibration reference
        idx = rng.choice(s.score.size, size=min(cfg.watch_bins_per_session, s.score.size),
                         replace=False)
        p = np.clip(conformal_pvalues(cal, s.score[idx], rng), 1e-12, 1.0)
        # log power-bet f_eps(p) = log eps + (eps-1) log p, summed over this session's p's
        cum = cum + (log_eps[:, None] + np.outer(eps - 1.0, np.log(p))).sum(axis=1)
        lm = float(_logsumexp(cum) - np.log(eps_grid))
        log_mart.append(lm)
        if stop is None and lm >= log_thr:
            stop = s.k
    return WatchTrace(log_mart=np.array(log_mart), threshold=log_thr, stop_session=stop)


def _logsumexp(x: np.ndarray) -> float:
    m = float(np.max(x))
    return m + float(np.log(np.sum(np.exp(x - m))))
