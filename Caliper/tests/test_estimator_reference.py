"""Tests for the numpy-only over-confident reference IVIM estimator."""
import numpy as np

from caliper import metrics as M
from caliper.estimator_reference import ReferenceIVIMEstimator
from caliper.forward import PARAM_NAMES, synthetic_cohort

LEVELS = np.array([0.05, 0.25, 0.5, 0.75, 0.95])


def test_predict_quantiles_shape_and_monotone():
    cohort = synthetic_cohort(n=200, snr=40.0, seed=0)
    est = ReferenceIVIMEstimator()
    q = est.predict_quantiles(cohort.signals, LEVELS)
    assert q.shape == (200, len(PARAM_NAMES), len(LEVELS))
    # quantiles strictly increasing along the level axis (non-degenerate width)
    assert np.all(np.diff(q, axis=2) > 0)


def test_point_recovers_D_and_f_at_high_snr():
    # At high SNR the segmented fit should recover the well-identified params.
    cohort = synthetic_cohort(n=1500, snr=200.0, seed=3)
    est = ReferenceIVIMEstimator()
    pt = est.predict_point(cohort.signals)
    err_D = np.abs(pt[:, 0] - cohort.params[:, 0])
    err_f = np.abs(pt[:, 1] - cohort.params[:, 1])
    assert np.median(err_D) < 0.15      # D within ~0.15 (1e-3 mm^2/s)
    assert np.median(err_f) < 0.05      # f within ~0.05


def test_raw_estimator_is_overconfident():
    cohort = synthetic_cohort(n=4000, snr=40.0, seed=2)
    est = ReferenceIVIMEstimator()
    q = est.predict_quantiles(cohort.signals, LEVELS)
    scores = M.score_quantiles(cohort.params, q, LEVELS, alpha=0.1,
                               param_names=PARAM_NAMES)
    # every parameter's 90% interval under-covers (reported sigma too narrow)
    for s in scores:
        assert s.coverage < 0.80, f"{s.name} coverage {s.coverage} not over-confident"


def test_dstar_error_grows_with_true_dstar():
    # The identifiability wall: D* point error is larger in the high-D* tercile.
    cohort = synthetic_cohort(n=6000, snr=40.0, seed=4)
    est = ReferenceIVIMEstimator()
    pt = est.predict_point(cohort.signals)
    true_dstar = cohort.params[:, 2]
    err = np.abs(pt[:, 2] - true_dstar)
    strata = M.tercile_groups(true_dstar)
    med_low = np.median(err[strata == 0])
    med_high = np.median(err[strata == 2])
    assert med_high > med_low, f"high-D* err {med_high} should exceed low-D* err {med_low}"


def test_rejects_wrong_signal_width():
    est = ReferenceIVIMEstimator()
    try:
        est.predict_point(np.zeros((10, 3)))  # wrong number of b-values
    except ValueError:
        return
    raise AssertionError("expected ValueError on mismatched b-value count")
