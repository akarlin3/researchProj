"""Synthetic test--retest harness -- Echo's method self-test (and the Reverb fallback spec).

This module generates controlled scan--rescan data with KNOWN truth, KNOWN measurement
noise, KNOWN systematic bias, and a deployed interval of KNOWN scale. It exists for two
reasons:

  1. METHOD SELF-TEST (run at CP1, fully synthetic/open): prove that
     ``echo_repeat.statistic`` recovers the analytic test--retest coverage, detects width
     mis-scaling, and is INVARIANT TO BIAS (the precision-not-accuracy guarantee).
  2. REVERB FALLBACK: if the public repeatability data turns out unsuitable for Echo's
     regime (CP2 data gate), this synthetic harness *is* the deliverable -- a constrained
     validation on controlled truth. See PROMOTION.md.

Everything here is seeded and synthetic; it touches no external or private data.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import statistic as st


@dataclass
class TestRetest:
    """One synthetic scan--rescan cohort with a deployed (1-level)-truth-coverage interval."""
    est_a: np.ndarray
    est_b: np.ndarray
    lo_a: np.ndarray
    hi_a: np.ndarray
    lo_b: np.ndarray
    hi_b: np.ndarray
    truth: np.ndarray
    sigma_meas: np.ndarray
    level: float
    bias: float
    width_scale: float
    model_to_meas_ratio: float


def simulate(n: int = 76, *, level: float = 0.10, bias: float = 0.0,
             width_scale: float = 1.0, model_to_meas_ratio: float = 0.0,
             heteroscedastic: bool = True, mu: float = 1.2e-3, cv_truth: float = 0.35,
             snr_floor: float = 8.0, seed: int = 0) -> TestRetest:
    """Generate ``n`` units, each with two repeat estimates and a deployed interval.

    Model: ``est_repeat = truth + bias + eps``, ``eps ~ N(0, sigma_meas^2)`` independent
    across repeats. The deployed interval half-width is
    ``hw = width_scale * z * sqrt(sigma_meas^2 + sigma_model^2)`` with
    ``sigma_model = model_to_meas_ratio * sigma_meas`` (non-repeating spread). The interval
    is centred on each repeat's own estimate -- exactly the deployment posture.

    Defaults evoke a tissue-diffusion D regime (mu ~ 1.2e-3 mm^2/s) at n=76 (ACRIN-6698).
    """
    rng = np.random.default_rng(seed)
    truth = rng.lognormal(mean=np.log(mu), sigma=cv_truth, size=n)
    # per-unit measurement SD: a CV floor plus heteroscedastic spread tied to magnitude
    base_cv = 1.0 / snr_floor
    if heteroscedastic:
        cv = base_cv * (1.0 + 0.8 * rng.random(n))
    else:
        cv = np.full(n, base_cv)
    sigma_meas = cv * truth
    sigma_model = model_to_meas_ratio * sigma_meas

    eps_a = rng.normal(0.0, 1.0, n) * sigma_meas
    eps_b = rng.normal(0.0, 1.0, n) * sigma_meas
    est_a = truth + bias + eps_a
    est_b = truth + bias + eps_b

    z = st.norm_ppf(1.0 - level / 2.0)
    hw = width_scale * z * np.sqrt(sigma_meas ** 2 + sigma_model ** 2)
    return TestRetest(
        est_a=est_a, est_b=est_b,
        lo_a=est_a - hw, hi_a=est_a + hw,
        lo_b=est_b - hw, hi_b=est_b + hw,
        truth=truth, sigma_meas=sigma_meas, level=level, bias=bias,
        width_scale=width_scale, model_to_meas_ratio=model_to_meas_ratio,
    )


def measure(tr: TestRetest, *, n_boot: int = 2000, seed: int = 0) -> st.CoverageResult:
    """Run Echo's scorecard on a synthetic cohort."""
    return st.evaluate("synthetic", tr.est_a, tr.est_b, tr.lo_a, tr.hi_a, tr.lo_b, tr.hi_b,
                       level=tr.level, n_boot=n_boot, seed=seed)


def self_test(n: int = 20000, level: float = 0.10, seed: int = 0) -> dict:
    """The locked CP1 method self-test.

    Returns a dict of checks (each ``{value, target, pass}``). Large ``n`` keeps Monte
    Carlo error small so we can assert against the analytic targets.
    """
    out: dict = {"n": n, "level": level}

    # (1) correctly-scaled interval recovers the analytic ~0.755 repeat-coverage
    tr1 = simulate(n, level=level, width_scale=1.0, model_to_meas_ratio=0.0,
                   heteroscedastic=True, seed=seed)
    cov1 = st.test_retest_coverage(tr1.est_a, tr1.est_b, tr1.lo_a, tr1.hi_a,
                                   tr1.lo_b, tr1.hi_b)
    tgt1 = st.analytic_repeat_coverage(level, scale=1.0)
    out["scaled_recovers_analytic"] = {
        "value": cov1, "target": tgt1, "pass": abs(cov1 - tgt1) < 0.01}

    # (2) BIAS INVARIANCE -- the precision-not-accuracy guarantee
    tr2 = simulate(n, level=level, width_scale=1.0, bias=5.0e-3, seed=seed)  # huge bias
    cov2 = st.test_retest_coverage(tr2.est_a, tr2.est_b, tr2.lo_a, tr2.hi_a,
                                   tr2.lo_b, tr2.hi_b)
    out["bias_invariance"] = {
        "value": cov2, "target": tgt1, "pass": abs(cov2 - tgt1) < 0.01}

    # (3) SCALE SENSITIVITY -- coverage tracks width mis-scaling per the analytic law
    scale_checks = []
    for s in (0.5, 0.75, 1.0, 1.5, 2.0):
        tr = simulate(n, level=level, width_scale=s, seed=seed + 1)
        cov = st.test_retest_coverage(tr.est_a, tr.est_b, tr.lo_a, tr.hi_a, tr.lo_b, tr.hi_b)
        tgt = st.analytic_repeat_coverage(level, scale=s)
        scale_checks.append({"scale": s, "value": cov, "target": tgt,
                             "pass": abs(cov - tgt) < 0.012})
    out["scale_sensitivity"] = scale_checks

    # (4) DISTINCTNESS FROM GAUGE -- a pure rescale leaves Spearman fixed, moves coverage
    base = simulate(n, level=level, width_scale=1.0, seed=seed + 2)
    width = 0.5 * ((base.hi_a - base.lo_a) + (base.hi_b - base.lo_b))
    dabs = np.abs(base.est_b - base.est_a)
    rho_base = st.spearman(width, dabs)
    # rescale all widths by 0.5 -> Spearman identical, coverage collapses
    hw2 = 0.5 * (base.hi_a - base.est_a)
    cov_base = st.test_retest_coverage(base.est_a, base.est_b, base.lo_a, base.hi_a,
                                       base.lo_b, base.hi_b)
    cov_resc = st.test_retest_coverage(base.est_a, base.est_b, base.est_a - hw2,
                                       base.est_a + hw2, base.est_b - hw2, base.est_b + hw2)
    rho_resc = st.spearman(0.5 * width, dabs)
    out["distinct_from_gauge"] = {
        "spearman_before": rho_base, "spearman_after_rescale": rho_resc,
        "coverage_before": cov_base, "coverage_after_rescale": cov_resc,
        "pass": abs(rho_base - rho_resc) < 1e-9 and abs(cov_base - cov_resc) > 0.05}

    checks = [out["scaled_recovers_analytic"]["pass"], out["bias_invariance"]["pass"],
              all(c["pass"] for c in out["scale_sensitivity"]),
              out["distinct_from_gauge"]["pass"]]
    out["ALL_PASS"] = bool(all(checks))
    return out
