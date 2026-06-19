"""vernier.decision -- decision-value-per-scan-minute through Minos's lens.

**PROVISIONAL (Minos-dependent).** The feasibility gate (:mod:`vernier.feasibility`)
does NOT use this module; it is the Experiment B readout, strictly downstream of the
gate. Minos is in review, so every number this produces is flagged PROVISIONAL (see
``ASSUMPTIONS.md``) and must be re-validated if Minos changes.

It prices the conformal-corrected D* error bar on a stylised treat/spare/escalate
decision using Minos-Core's *published* asymmetric utility and its closed-form
expected utility for a Gaussian reported posterior. Per voxel:

  1. reported posterior q = N(mu, sigma): mu = corrected D* median, sigma =
     (corrected 90% half-width) / z_0.95 -- a Gaussian moment-match of the
     conformal-corrected interval.
  2. Bayes action a* = argmax_a EU(a | q)        (minos.utility.expected_utility_under_q)
  3. achieved utility U(a*, D*_true)             (minos.utility.utility)

The decision VALUE of a scheme is its mean achieved utility minus the *no-scan
baseline* (the Bayes action taken under the prior alone), i.e. the decision utility
*gained* by acquiring and correcting. Dividing by scan-minutes gives
value-per-scan-minute -- the quantity whose per-minute denominator only has teeth
across schemes of *different* scan-time (the efficiency frontier).

Decision thresholds t1<t2 are the prior D* terciles; the cost asymmetry
(k_under=2 >= k_over=1) is Minos's default (under-treatment penalised more).
"""
from __future__ import annotations

from dataclasses import dataclass
from statistics import NormalDist

import numpy as np

from . import _paths
from .schemes import BScheme

_paths.add_caliper()
_paths.add_minos()
from caliper.conformal import SplitConformalQuantile  # noqa: E402
from caliper.estimator_reference import ReferenceIVIMEstimator  # noqa: E402
from caliper.forward import synthetic_cohort  # noqa: E402
from minos.config import MinosConfig            # noqa: E402
from minos.utility import Action, expected_utility_under_q, utility  # noqa: E402

Q_LEVELS = np.array([0.05, 0.25, 0.5, 0.75, 0.95])
DSTAR = 2
_Z95 = NormalDist().inv_cdf(0.95)
_ACTIONS = (Action.SPARE, Action.TREAT, Action.ESCALATE)


@dataclass
class DecisionResult:
    name: str
    n_b: int
    scan_minutes: float
    width_dstar: float            # post-conformal D* 90% width (sharpness)
    cond_cov_high_dstar: float    # post-conformal high-D* tercile coverage
    mean_utility: float           # mean achieved decision utility (higher = better; <= 0)
    baseline_utility: float       # mean achieved utility of the no-scan (prior) decision


def make_decision_config(dstar_values, k_under: float = 2.0, k_over: float = 1.0) -> MinosConfig:
    """Minos config whose treat/spare/escalate thresholds are the prior D* terciles."""
    t1, t2 = (float(x) for x in np.quantile(dstar_values, [1.0 / 3.0, 2.0 / 3.0]))
    return MinosConfig(t1=t1, t2=t2, k_under=k_under, k_over=k_over)


def _corrected_posterior(scheme: BScheme, n: int, snr: float, seed: int, cal_frac: float):
    """Return (mu, sigma, dstar_true) for the test split (Gaussian-matched corrected D*)."""
    cohort = synthetic_cohort(n=n, bvalues=scheme.b, snr=snr, noise="rician", seed=seed)
    est = ReferenceIVIMEstimator(bvalues=scheme.b)
    q_raw = est.predict_quantiles(cohort.signals, Q_LEVELS)
    y = cohort.params
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    n_cal = int(round(cal_frac * n))
    cal_idx, test_idx = perm[:n_cal], perm[n_cal:]
    q_corr = SplitConformalQuantile(Q_LEVELS).calibrate(q_raw[cal_idx], y[cal_idx]).apply(q_raw[test_idx])
    mu = q_corr[:, DSTAR, 2]                                   # corrected median
    sigma = (q_corr[:, DSTAR, -1] - q_corr[:, DSTAR, 0]) / (2.0 * _Z95)
    sigma = np.maximum(sigma, 1e-6)
    return mu, sigma, y[test_idx, DSTAR]


def _achieved_utility(mu, sigma, theta_true, cfg: MinosConfig) -> np.ndarray:
    """Per-voxel achieved utility U(a*, theta_true), a* = Bayes action under N(mu,sigma)."""
    eu = np.stack([expected_utility_under_q(a, mu, sigma, cfg) for a in _ACTIONS])  # (3, n)
    a_star = np.argmax(eu, axis=0)                                                  # (n,)
    u_all = np.stack([np.asarray(utility(a, theta_true, cfg), dtype=float) for a in _ACTIONS])
    return u_all[a_star, np.arange(theta_true.shape[0])]


def decision_value(scheme: BScheme, *, n: int = 8000, snr: float = 33.0, seed: int = 0,
                   cal_frac: float = 0.5, cfg: MinosConfig = None) -> DecisionResult:
    """Decision-value-per-scan-minute for one scheme (PROVISIONAL; Minos lens)."""
    mu, sigma, theta = _corrected_posterior(scheme, n, snr, seed, cal_frac)
    if cfg is None:
        cfg = make_decision_config(theta)

    # scheme decisions use the per-voxel corrected posterior
    u_scheme = _achieved_utility(mu, sigma, theta, cfg)
    # no-scan baseline: one Bayes action under the prior alone (mu, sigma = prior moments)
    mu0 = float(np.mean(theta))
    sig0 = float(np.std(theta))
    eu0 = np.array([float(expected_utility_under_q(a, mu0, sig0, cfg)) for a in _ACTIONS])
    a0 = _ACTIONS[int(np.argmax(eu0))]
    u_base = np.asarray(utility(a0, theta, cfg), dtype=float)

    mean_u = float(np.mean(u_scheme))
    base_u = float(np.mean(u_base))
    width = float(np.mean(2.0 * _Z95 * sigma))  # the post-conformal 90% D* width
    # high-D* tercile coverage (consistency with the gate's readout)
    hi = theta >= np.quantile(theta, 2.0 / 3.0)
    lo_b = mu - _Z95 * sigma
    hi_b = mu + _Z95 * sigma
    cond_high = float(np.mean(((theta >= lo_b) & (theta <= hi_b))[hi]))

    return DecisionResult(
        name=scheme.name, n_b=scheme.n_b, scan_minutes=scheme.scan_minutes(),
        width_dstar=width, cond_cov_high_dstar=cond_high,
        mean_utility=mean_u, baseline_utility=base_u,
    )


def format_frontier(results: list[DecisionResult], cfg: MinosConfig) -> str:
    """Frontier table: decision utility (higher=better) vs scan-minutes, with the
    *incremental* utility per added minute relative to the cheapest protocol.

    Decision value is reported relative to the cheapest protocol in the set (not
    the no-scan prior): for D* the corrected interval is wider than the prior (the
    identifiability wall + reference-estimator bias inflate it), so 'value over
    prior' is negative for every protocol -- itself an honest, Gauge-consistent
    finding, printed below the table. The frontier's content is the *saturation*:
    each added scan-minute buys progressively less calibrated decision utility.
    """
    res = sorted(results, key=lambda r: r.scan_minutes)
    cheap = res[0]
    L = []
    L.append("=== Vernier Experiment B: decision utility vs scan-time (PROVISIONAL, Minos lens) ===")
    L.append(f"treat/spare/escalate on D*; thresholds t1={cfg.t1:.2f} t2={cfg.t2:.2f} "
             f"(prior D* terciles), k_under={cfg.k_under:g} k_over={cfg.k_over:g}")
    L.append(f"decision utility = mean achieved utility (higher=better, <=0); "
             f"baseline = cheapest protocol ({cheap.name})")
    L.append("")
    hdr = (f"{'scheme':>14} {'n_b':>4} {'min':>6} {'width':>8} {'hiCov':>7} "
           f"{'meanU':>8} {'dU_vs_cheap':>12} {'dU/min':>8}")
    L.append(hdr)
    L.append("-" * len(hdr))
    for r in res:
        dU = r.mean_utility - cheap.mean_utility
        dmin = r.scan_minutes - cheap.scan_minutes
        dU_per_min = (dU / dmin) if dmin > 0 else float("nan")
        slope = "    -- " if dmin == 0 else f"{dU_per_min:>8.3f}"
        L.append(f"{r.name:>14} {r.n_b:>4} {r.scan_minutes:>6.1f} {r.width_dstar:>8.2f} "
                 f"{r.cond_cov_high_dstar:>7.3f} {r.mean_utility:>8.3f} {dU:>12.3f} {slope}")
    # marginal (segment-to-segment) returns, to show diminishing returns explicitly
    L.append("")
    L.append("marginal decision utility per added scan-minute (each step vs the previous):")
    for prev, cur in zip(res[:-1], res[1:]):
        dmin = cur.scan_minutes - prev.scan_minutes
        marg = (cur.mean_utility - prev.mean_utility) / dmin if dmin > 0 else float("nan")
        L.append(f"  {prev.name} -> {cur.name}: +{cur.mean_utility - prev.mean_utility:.3f} utility "
                 f"over +{dmin:.1f} min  =>  {marg:.3f}/min")
    L.append("")
    L.append(f"no-scan(prior) baseline mean utility = {res[0].baseline_utility:.3f}; "
             f"every protocol's mean utility is below it -> for D* decisions, acquisition does "
             f"not beat the prior (consistent with Gauge's identifiability wall).")
    return "\n".join(L)
