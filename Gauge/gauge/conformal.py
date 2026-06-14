"""Split-conformal and conformalized quantile regression (CQR).

Both methods turn a held-out *calibration* set into a finite-sample, distribution-
free prediction interval: under exchangeability of calibration and test points,
the marginal coverage satisfies

    1 - alpha  <=  P(Y in C(X))  <=  1 - alpha + 1/(n_cal + 1).

The guarantee holds for ANY base predictor -- a wrong/biased base predictor only
widens the interval, it does not break coverage. That is what makes conformal a
correctness tool rather than a heuristic.

References
----------
Vovk et al. (2005); Lei et al. (2018) split conformal; Romano, Patterson & Candes
(2019) "Conformalized Quantile Regression".
"""
import numpy as np


def conformal_quantile(scores, alpha):
    """The conformal quantile of nonconformity ``scores`` at level ``alpha``.

    Returns the k-th smallest score with k = ceil((n + 1) * (1 - alpha)). When
    k > n (calibration set too small for the requested alpha) there is no finite
    quantile and the interval must be infinite, so ``np.inf`` is returned. This
    +1 correction is what delivers the finite-sample (not merely asymptotic)
    coverage guarantee.
    """
    scores = np.asarray(scores, dtype=float)
    n = scores.size
    k = int(np.ceil((n + 1) * (1.0 - alpha)))
    if k > n:
        return np.inf
    return float(np.sort(scores)[k - 1])


def split_conformal(cal_pred, cal_true, test_pred, alpha):
    """Symmetric absolute-residual split conformal.

    Nonconformity score on calibration: s_i = |y_i - yhat_i|. The calibrated
    radius q is the conformal quantile of those scores; the test interval is
    ``[test_pred - q, test_pred + q]``.

    Returns ``(lower, upper, q)``.
    """
    cal_pred = np.asarray(cal_pred, dtype=float)
    cal_true = np.asarray(cal_true, dtype=float)
    test_pred = np.asarray(test_pred, dtype=float)
    scores = np.abs(cal_true - cal_pred)
    q = conformal_quantile(scores, alpha)
    return test_pred - q, test_pred + q, q


def cqr(cal_lo, cal_hi, cal_true, test_lo, test_hi, alpha):
    """Conformalized Quantile Regression (Romano et al. 2019).

    Given base lower/upper conditional-quantile predictions, the nonconformity
    score is the signed distance outside the predicted band:
    ``E_i = max(qlo_i - y_i, y_i - qhi_i)``. The calibrated offset q (conformal
    quantile of E) inflates (or, if negative, shrinks) the band:
    ``[test_lo - q, test_hi + q]``.

    Returns ``(lower, upper, q)``.
    """
    cal_lo = np.asarray(cal_lo, dtype=float)
    cal_hi = np.asarray(cal_hi, dtype=float)
    cal_true = np.asarray(cal_true, dtype=float)
    test_lo = np.asarray(test_lo, dtype=float)
    test_hi = np.asarray(test_hi, dtype=float)
    scores = np.maximum(cal_lo - cal_true, cal_true - cal_hi)
    q = conformal_quantile(scores, alpha)
    return test_lo - q, test_hi + q, q


def empirical_coverage(lower, upper, true):
    """Fraction of points whose truth lies within [lower, upper] (inclusive)."""
    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    true = np.asarray(true, dtype=float)
    return float(np.mean((true >= lower) & (true <= upper)))


def interval_width(lower, upper):
    """Per-point interval width upper - lower."""
    return np.asarray(upper, dtype=float) - np.asarray(lower, dtype=float)


def interval_score(lower, upper, true, alpha):
    """Interval score for a central (1 - alpha) prediction interval.

    Gneiting & Raftery (2007): a proper scoring rule for the (alpha/2, 1-alpha/2)
    interval. Lower is better -- it rewards sharp intervals but adds a 2/alpha
    penalty proportional to how far the truth falls outside the band:

        IS = (u - l) + (2/alpha)(l - y) 1{y<l} + (2/alpha)(y - u) 1{y>u}.

    This lets the benchmark score raw, conformal, and conformalized-model-based
    methods on the same footing (coverage and sharpness combined).
    """
    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    true = np.asarray(true, dtype=float)
    width = upper - lower
    below = np.maximum(lower - true, 0.0)
    above = np.maximum(true - upper, 0.0)
    return width + (2.0 / alpha) * (below + above)
