"""Tests for the split-conformal / CQR wrapper (numpy-only)."""
from statistics import NormalDist

import numpy as np
import pytest

from caliper import conformal as C
from caliper import metrics as M
from caliper.estimator_reference import ReferenceIVIMEstimator
from caliper.forward import PARAM_NAMES, synthetic_cohort

LEVELS = np.array([0.05, 0.25, 0.5, 0.75, 0.95])
_Z = np.array([NormalDist().inv_cdf(p) for p in LEVELS])
ALPHA = 0.10
DSTAR = 2  # index of D* in (D, f, D*)


def _ivim_split(snr=40.0, n_cal=4000, n_test=9000):
    """Real reference-estimator run on synthetic IVIM: (cal, test, q_cal, q_test)."""
    cal = synthetic_cohort(n=n_cal, snr=snr, seed=1)
    test = synthetic_cohort(n=n_test, snr=snr, seed=2)
    est = ReferenceIVIMEstimator()
    return cal, test, est


def _overconfident(n, rng, spread=0.5):
    mu = rng.normal(size=(n, 1))
    y = mu + rng.normal(size=(n, 1))               # true sigma = 1
    q = mu[:, :, None] + spread * _Z[None, None, :]  # under-dispersed
    return y, q


def test_conformity_scores_definition():
    s = C.conformity_scores(np.array([0.0]), np.array([2.0]), np.array([3.0]))
    assert s[0] == pytest.approx(1.0)  # y above hi by 1
    s = C.conformity_scores(np.array([0.0]), np.array([2.0]), np.array([-1.0]))
    assert s[0] == pytest.approx(1.0)  # y below lo by 1
    s = C.conformity_scores(np.array([0.0]), np.array([2.0]), np.array([1.0]))
    assert s[0] == pytest.approx(-1.0)  # inside -> negative


def test_offset_increases_with_lower_alpha():
    rng = np.random.default_rng(0)
    scores = rng.normal(size=1000)
    q90 = C.conformal_offset(scores, alpha=0.1)
    q50 = C.conformal_offset(scores, alpha=0.5)
    assert q90 > q50


def test_conformal_restores_marginal_coverage():
    rng = np.random.default_rng(0)
    y_cal, q_cal = _overconfident(3000, rng)
    y_te, q_te = _overconfident(5000, rng)
    raw = M.score_quantiles(y_te, q_te, LEVELS, alpha=0.1)[0]
    cq = C.SplitConformalQuantile(LEVELS).calibrate(q_cal, y_cal)
    q_corr = cq.apply(q_te)
    cor = M.score_quantiles(y_te, q_corr, LEVELS, alpha=0.1)[0]
    assert raw.coverage < 0.75              # raw is over-confident
    assert abs(cor.coverage_gap) < 0.03     # conformal restores nominal


def test_apply_preserves_monotonic_quantiles():
    rng = np.random.default_rng(1)
    y_cal, q_cal = _overconfident(2000, rng)
    _, q_te = _overconfident(1000, rng)
    cq = C.SplitConformalQuantile(LEVELS).calibrate(q_cal, y_cal)
    q_corr = cq.apply(q_te)
    assert np.all(np.diff(q_corr, axis=2) >= -1e-9)


def test_apply_before_calibrate_raises():
    cq = C.SplitConformalQuantile(LEVELS)
    with pytest.raises(RuntimeError):
        cq.apply(np.zeros((4, 1, 5)))


def test_calibrate_shape_validation():
    cq = C.SplitConformalQuantile(LEVELS)
    with pytest.raises(ValueError):
        cq.calibrate(np.zeros((4, 2, 5)), np.zeros((4, 3)))


# --------------------------------------------------------------------------- #
# Real IVIM runs: the headline marginal / conditional results
# --------------------------------------------------------------------------- #
def test_cqr_restores_marginal_coverage_on_ivim():
    """CQR restores marginal coverage to within ~2 points, per IVIM parameter."""
    cal, test, est = _ivim_split()
    q_cal = est.predict_quantiles(cal.signals, LEVELS)
    q_test = est.predict_quantiles(test.signals, LEVELS)
    raw = M.score_quantiles(test.params, q_test, LEVELS, alpha=ALPHA,
                            param_names=PARAM_NAMES)
    cq = C.SplitConformalQuantile(LEVELS).calibrate(q_cal, cal.params)
    cor = M.score_quantiles(test.params, cq.apply(q_test), LEVELS, alpha=ALPHA,
                            param_names=PARAM_NAMES)
    for r, c in zip(raw, cor):
        assert r.coverage < 0.80, f"{r.name} raw not over-confident ({r.coverage})"
        assert abs(c.coverage_gap) <= 0.025, f"{c.name} post-CQR gap {c.coverage_gap}"


def test_split_residual_restores_marginal_coverage_on_ivim():
    cal, test, est = _ivim_split()
    pt_cal = est.predict_point(cal.signals)
    pt_test = est.predict_point(test.signals)
    sr = C.SplitConformalResidual(alpha=ALPHA).calibrate(pt_cal, cal.params)
    lo, hi = sr.apply(pt_test)
    for j, name in enumerate(PARAM_NAMES):
        cov = M.empirical_coverage(test.params[:, j], lo[:, j], hi[:, j])
        assert abs(cov - 0.90) <= 0.025, f"{name} split-residual coverage {cov}"


def test_marginal_cqr_leaves_conditional_gap_high_dstar():
    """Marginal CQR fixes the pool but not the conditional D* coverage:
    low-D* ends over-covered, high-D* under-covered."""
    cal, test, est = _ivim_split()
    q_cal = est.predict_quantiles(cal.signals, LEVELS)
    q_test = est.predict_quantiles(test.signals, LEVELS)
    cq = C.SplitConformalQuantile(LEVELS).calibrate(q_cal, cal.params)
    q_corr = cq.apply(q_test)
    lo, hi = M.central_interval(q_corr[:, DSTAR, :], LEVELS, ALPHA)
    strata = M.tercile_groups(test.params[:, DSTAR])
    by = C.conditional_coverage_by_strata(test.params[:, DSTAR], lo, hi, strata)
    # pooled coverage is restored ...
    pooled = M.empirical_coverage(test.params[:, DSTAR], lo, hi)
    assert abs(pooled - 0.90) <= 0.03
    # ... but conditional coverage is not uniform: low over-covers, high under-covers
    assert by[0].coverage > 0.92, f"low-D* not over-covered: {by[0].coverage}"
    assert by[2].coverage < 0.90, f"high-D* not under-covered: {by[2].coverage}"
    assert by[0].coverage - by[2].coverage >= 0.03


def test_mondrian_restores_conditional_coverage_at_width_cost():
    """Mondrian CQR equalizes per-tercile coverage but inflates high-D* width."""
    cal, test, est = _ivim_split()
    q_cal = est.predict_quantiles(cal.signals, LEVELS)
    q_test = est.predict_quantiles(test.signals, LEVELS)
    groups_cal = M.tercile_groups(cal.params[:, DSTAR])
    strata = M.tercile_groups(test.params[:, DSTAR])
    mq = C.MondrianConformalQuantile(LEVELS).calibrate(q_cal, cal.params, groups_cal)
    q_corr = mq.apply(q_test, strata)
    lo, hi = M.central_interval(q_corr[:, DSTAR, :], LEVELS, ALPHA)
    by = C.conditional_coverage_by_strata(test.params[:, DSTAR], lo, hi, strata)
    # every tercile is near nominal (conditional coverage delivered)
    for g in (0, 1, 2):
        assert abs(by[g].coverage - 0.90) <= 0.04, f"g{g} coverage {by[g].coverage}"
    # but at a real width cost: high-D* intervals far wider than low-D*
    ratio = by[2].mean_width / by[0].mean_width
    assert ratio > 1.5, f"expected high/low width ratio > 1.5, got {ratio}"


def test_conditional_coverage_by_strata_values():
    y = np.array([0.0, 0.0, 5.0, 5.0])
    lo = np.array([-1.0, -1.0, 0.0, 0.0])      # stratum 0 covers, stratum 1 doesn't
    hi = np.array([1.0, 1.0, 4.0, 4.0])
    strata = np.array([0, 0, 1, 1])
    out = C.conditional_coverage_by_strata(y, lo, hi, strata)
    assert out[0].coverage == 1.0 and out[0].mean_width == 2.0 and out[0].n == 2
    assert out[1].coverage == 0.0 and out[1].mean_width == 4.0 and out[1].n == 2


def test_mondrian_falls_back_for_unseen_group():
    cal, test, est = _ivim_split(n_cal=600, n_test=600)
    q_cal = est.predict_quantiles(cal.signals, LEVELS)
    q_test = est.predict_quantiles(test.signals, LEVELS)
    groups_cal = np.zeros(len(cal), dtype=int)   # only group 0 seen in calibration
    groups_test = np.ones(len(test), dtype=int)  # only group 1 at test time
    mq = C.MondrianConformalQuantile(LEVELS).calibrate(q_cal, cal.params, groups_cal)
    q_corr = mq.apply(q_test, groups_test)       # must not raise; uses global fallback
    assert q_corr.shape == q_test.shape
