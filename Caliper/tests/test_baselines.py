"""Tests for caliper.baselines -- the box-constrained NLLS IVIM baseline.

These pin the two qualitative Fashion findings the baseline exists to reproduce
on synthetic data: (1) D* rails against the box bounds at a clearly nonzero rate,
and (2) the calibration ruler flags the NLLS D* intervals as overconfident
(coverage below nominal). All data is in-repo synthetic; no clinical numbers.
"""
import numpy as np
import pytest

from caliper import metrics as M
from caliper.baselines import NLLSIVIMEstimator
from caliper.forward import PARAM_NAMES, synthetic_cohort

LEVELS = np.array([0.05, 0.25, 0.5, 0.75, 0.95])


def _score(cohort, est):
    q = est.predict_quantiles(cohort.signals, LEVELS)
    scores = M.score_quantiles(
        cohort.params, q, LEVELS, alpha=0.1,
        param_names=PARAM_NAMES, conditioning=cohort.params,
    )
    return {s.name: s for s in scores}


# --------------------------------------------------------------------------- #
# Contract / shape / monotonicity
# --------------------------------------------------------------------------- #
def test_predict_quantiles_shape_and_monotone():
    c = synthetic_cohort(n=120, snr=30.0, seed=0)
    est = NLLSIVIMEstimator()
    q = est.predict_quantiles(c.signals, LEVELS)
    assert q.shape == (120, 3, LEVELS.size)
    # quantiles strictly non-decreasing along the level axis
    assert np.all(np.diff(q, axis=2) >= -1e-9)


def test_predict_point_shape_matches_truth_columns():
    c = synthetic_cohort(n=80, snr=40.0, seed=1)
    point = NLLSIVIMEstimator().predict_point(c.signals)
    assert point.shape == (80, 3)  # (D, f, D*)


def test_quantiles_reject_out_of_range_levels():
    c = synthetic_cohort(n=10, snr=40.0, seed=0)
    est = NLLSIVIMEstimator()
    with pytest.raises(ValueError):
        est.predict_quantiles(c.signals, np.array([0.0, 0.5, 1.0]))


def test_solve_rejects_wrong_bvalue_count():
    est = NLLSIVIMEstimator()
    with pytest.raises(ValueError):
        est.solve(np.zeros((5, 3)))  # 3 != 11 b-values


def test_bad_bounds_rejected():
    with pytest.raises(ValueError):
        NLLSIVIMEstimator(lower=(0.1, 0.1, 0.0, 3.0), upper=(2.0, 4.0, 1.0, 3.0))


# --------------------------------------------------------------------------- #
# Point-estimate sanity: well-identified params recover at high SNR
# --------------------------------------------------------------------------- #
def test_high_snr_recovers_D_and_f():
    c = synthetic_cohort(n=400, snr=80.0, seed=3)
    point = NLLSIVIMEstimator().predict_point(c.signals)
    med_abs_D = np.median(np.abs(point[:, 0] - c.params[:, 0]))
    med_abs_f = np.median(np.abs(point[:, 1] - c.params[:, 1]))
    assert med_abs_D < 0.15, med_abs_D     # D within ~0.15e-3 mm^2/s
    assert med_abs_f < 0.05, med_abs_f     # f within ~0.05


# --------------------------------------------------------------------------- #
# CLAIM 2a -- D* boundary railing at a clearly nonzero (but non-clinical) rate
# --------------------------------------------------------------------------- #
def test_dstar_rails_at_nonzero_rate():
    c = synthetic_cohort(n=600, snr=20.0, seed=0)
    est = NLLSIVIMEstimator()
    fit = est.solve(c.signals)
    rate = est.railing_rate_from_fit(fit, "Dstar")
    # clearly nonzero, yet nowhere near the clinical figure (reproduce the
    # phenomenon, not the number): a wide, defensible band.
    assert 0.03 < rate < 0.40, rate


def test_railed_mask_consistent_with_bounds():
    c = synthetic_cohort(n=300, snr=15.0, seed=2)
    est = NLLSIVIMEstimator()
    fit = est.solve(c.signals)
    lo, hi = est._lo[3], est._hi[3]
    span = hi - lo
    dstar = fit.params[:, 2]
    on_bound = (dstar - lo <= est.rail_tol * span) | (hi - dstar <= est.rail_tol * span)
    assert np.array_equal(on_bound, fit.dstar_railed)


def test_railing_increases_as_snr_drops():
    est = NLLSIVIMEstimator()
    hi_snr = est.boundary_railing_rate(
        synthetic_cohort(n=500, snr=60.0, seed=1).signals, "Dstar")
    lo_snr = est.boundary_railing_rate(
        synthetic_cohort(n=500, snr=12.0, seed=1).signals, "Dstar")
    assert lo_snr > hi_snr, (lo_snr, hi_snr)


def test_railing_rate_validates_param():
    c = synthetic_cohort(n=20, snr=40.0, seed=0)
    est = NLLSIVIMEstimator()
    with pytest.raises(ValueError):
        est.boundary_railing_rate(c.signals, "nope")


# --------------------------------------------------------------------------- #
# CLAIM 2b -- the ruler flags NLLS D* intervals as overconfident (under-covers)
# --------------------------------------------------------------------------- #
def test_ruler_flags_dstar_overconfident():
    c = synthetic_cohort(n=600, snr=20.0, seed=0)
    est = NLLSIVIMEstimator()
    s = _score(c, est)
    # D* central 90% interval covers well below nominal -> negative gap.
    assert s["Dstar"].coverage < 0.90
    assert s["Dstar"].coverage_gap < -0.03, s["Dstar"].coverage_gap
    # the well-identified params are far closer to nominal than D* is.
    assert s["Dstar"].coverage_gap < s["D"].coverage_gap
    assert s["Dstar"].coverage_gap < s["f"].coverage_gap


# --------------------------------------------------------------------------- #
# Retooled SD convention -- honest (default) vs the floored illustration
# --------------------------------------------------------------------------- #
def test_sd_convention_defaults_to_honest():
    # The retooled Fashion reports under the honest CRLB; honest is the default.
    assert NLLSIVIMEstimator().sd_convention == "honest"


def test_invalid_sd_convention_rejected():
    with pytest.raises(ValueError):
        NLLSIVIMEstimator(sd_convention="overconfident")
    with pytest.raises(ValueError):
        NLLSIVIMEstimator(railed_sd_floor=0.0)


def test_floored_convention_narrows_railed_dstar_and_lowers_coverage():
    # Honest CRLB (default) reports a *wide* SD for a railed/unidentified D*; the
    # floored "overconfident" convention overwrites it with the narrow floor. The
    # honest path is left untouched (so downstream numbers don't move); the floored
    # path under-covers D* harder -- the choice that reconstructs the now-dropped
    # marginal severity, kept only as a labelled illustration.
    c = synthetic_cohort(n=600, snr=20.0, seed=0)
    honest = NLLSIVIMEstimator()                       # default sd_convention
    floored = NLLSIVIMEstimator(sd_convention="floored")

    fit_h = honest.solve(c.signals)
    fit_f = floored.solve(c.signals)
    railed = fit_h.dstar_railed
    assert railed.any(), "need some railed D* voxels for this test"

    # point estimates and railing flags are identical -- only the SD differs.
    assert np.array_equal(fit_h.params, fit_f.params)
    assert np.array_equal(fit_h.dstar_railed, fit_f.dstar_railed)
    # honest railed-D* SD is wide; floored clamps it to the narrow floor (3.0).
    assert np.all(fit_f.sigma[railed, 2] == floored.railed_sd_floor)
    assert np.median(fit_h.sigma[railed, 2]) > floored.railed_sd_floor

    # consequence on the ruler: floored under-covers D* at least as hard as honest,
    # while honest itself still under-covers (the kept, milder finding).
    s_h = _score(c, honest)
    s_f = _score(c, floored)
    assert s_f["Dstar"].coverage <= s_h["Dstar"].coverage + 1e-9
    assert s_h["Dstar"].coverage_gap < 0.0
