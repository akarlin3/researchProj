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


# --------------------------------------------------------------------------- #
# Feature-conditional conformal toolkit (Gauge 03).
#
# Split/Mondrian/CQR above give MARGINAL (or group-marginal) coverage. The two
# methods below target *approximate conditional* coverage -- coverage that holds
# locally in feature space -- which is what the Gauge 03 high-D* attack needs.
# --------------------------------------------------------------------------- #
def weighted_conformal_quantile(scores, weights, alpha, self_weight=None):
    """Test-point-augmented weighted conformal quantile.

    Given calibration ``scores`` with nonnegative localization ``weights`` and
    the test point's own ``self_weight`` (its score is unknown, so it sits at
    +inf), return the smallest threshold t at which the normalized weight of
    ``{s_i <= t}`` reaches ``1 - alpha``. The +inf test-point mass is what keeps
    the localized interval valid under weighted exchangeability (Tibshirani et
    al. 2019; Guan 2023). With uniform weights this reduces exactly to
    :func:`conformal_quantile`. Returns +inf when 1-alpha is unreachable.
    """
    scores = np.asarray(scores, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if self_weight is None:
        self_weight = float(weights.max()) if weights.size else 1.0
    order = np.argsort(scores, kind="mergesort")
    s_sorted = scores[order]
    total = weights.sum() + self_weight
    if total <= 0:
        return np.inf
    cum = np.cumsum(weights[order]) / total
    idx = int(np.searchsorted(cum, 1.0 - alpha, side="left"))
    if idx >= s_sorted.size:
        return np.inf
    return float(s_sorted[idx])


def localized_conformal(cal_scores, cal_feat, test_feat, alpha, bandwidth,
                        kernel="gaussian"):
    """Localized conformal prediction (Guan 2023, LCP).

    For each test point, calibration nonconformity scores are reweighted by a
    kernel on standardized feature distance, then a test-point-augmented weighted
    quantile gives a *locally* calibrated threshold. Coverage adapts to feature
    space without ever using the label -- the label-free way to condition on a
    regime through features that correlate with it. Features are standardized by
    the calibration mean/std (deterministic). Returns an (n_test,) threshold
    array.
    """
    cal_scores = np.asarray(cal_scores, dtype=float)
    cal_feat = np.atleast_2d(np.asarray(cal_feat, dtype=float))
    test_feat = np.atleast_2d(np.asarray(test_feat, dtype=float))
    mu = cal_feat.mean(0)
    sd = cal_feat.std(0) + 1e-12
    C = (cal_feat - mu) / sd
    T = (test_feat - mu) / sd
    h = float(bandwidth)
    thr = np.empty(T.shape[0])
    for i in range(T.shape[0]):
        d2 = np.sum((C - T[i]) ** 2, axis=1)
        if kernel == "gaussian":
            w = np.exp(-0.5 * d2 / (h * h))
        elif kernel == "boxcar":
            w = (np.sqrt(d2) <= h).astype(float)
        else:
            raise ValueError(f"unknown kernel {kernel!r}")
        thr[i] = weighted_conformal_quantile(cal_scores, w, alpha,
                                             self_weight=1.0)
    return thr


def conditional_conformal(cal_scores, cal_Phi, test_Phi, alpha, ridge=0.0):
    """Conditional conformal via quantile regression of scores on features.

    The Gibbs, Cherian & Candes (2023) construction: over the linear class
    ``{x -> Phi(x)^T beta}``, the threshold that gives coverage conditional on
    that class is the level-(1-alpha) quantile regression of the nonconformity
    score on ``Phi``. We fit it by pinball loss (sklearn ``QuantileRegressor``,
    LP solver) and return ``t(x) = Phi(x)^T beta_hat`` per test point. This is
    the asymptotic-guarantee form (the finite-sample per-point augmentation is
    omitted, as is standard for moderate n). Coverage is then controlled
    conditional on the span of ``Phi`` -- e.g. on a plug-in D-hat* feature --
    but NOT on quantities outside that span.

    For symmetric absolute-residual scores the caller should clip the threshold
    at 0; for CQR scores a negative threshold legitimately shrinks the band.
    """
    from sklearn.linear_model import QuantileRegressor
    cal_scores = np.asarray(cal_scores, dtype=float)
    cal_Phi = np.atleast_2d(np.asarray(cal_Phi, dtype=float))
    test_Phi = np.atleast_2d(np.asarray(test_Phi, dtype=float))
    qr = QuantileRegressor(quantile=1.0 - alpha, alpha=ridge, solver="highs")
    qr.fit(cal_Phi, cal_scores)
    return qr.predict(test_Phi)
