"""Echo's repeatability-based validation statistics.

Echo asks a *scale-calibration* question that Gauge's published width-tracks-repeatability
check (Gauge paper Sec. 4.2.2 -- the "Sec. 3.7" of the protocol) does **not** ask:

    Gauge:  does the conformal interval *width* RANK-track scan--rescan variability?
            (Spearman correlation of width vs |Delta|; ordering only.)
    Echo:   is the conformal interval the right *SIZE* to capture that variability?
            (a coverage rate -- does one repeat's estimate fall inside the other
            repeat's deployed interval -- reported with a bootstrap CI.)

These are mathematically independent: rescaling every width by a constant leaves the
Spearman correlation invariant but moves Echo's coverage arbitrarily. See
``test_retest_coverage`` (the headline) and ``standardized_residual_dispersion`` (the
companion scale statistic).

THE LEGITIMACY RAZOR -- precision, not accuracy. Write each estimate as
``est = theta_true + bias + eps`` (eps = measurement noise). The test--retest
discrepancy ``Delta = est_B - est_A = eps_B - eps_A`` cancels ``bias`` exactly. So
every statistic in this module is *invariant to any systematic error common to both
repeats* (proven empirically in ``echo_repeat.harness`` and ``tests/``). Echo therefore
certifies that an interval is correctly **sized to measurement irreproducibility**
(precision); it is provably **blind to accuracy/bias** and makes **no ground-truth
(calibration) coverage claim**.

numpy-only by design (mirrors Caliper's posture); no scipy dependency.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

__all__ = [
    "norm_cdf",
    "norm_ppf",
    "analytic_repeat_coverage",
    "test_retest_coverage",
    "standardized_residual_dispersion",
    "spearman",
    "bca_ci",
    "CoverageResult",
    "evaluate",
]


# --------------------------------------------------------------------------------------
# numpy-only normal CDF / inverse-CDF (no scipy)
# --------------------------------------------------------------------------------------
def norm_cdf(x: np.ndarray | float) -> np.ndarray | float:
    """Standard normal CDF via ``math.erf`` (vectorised)."""
    x = np.asarray(x, dtype=float)
    out = 0.5 * (1.0 + np.vectorize(math.erf)(x / math.sqrt(2.0)))
    return float(out) if out.ndim == 0 else out


def norm_ppf(p: float) -> float:
    """Standard normal inverse CDF (Acklam's rational approximation; |abs err| < 1.2e-9)."""
    if not (0.0 < p < 1.0):
        if p == 0.0:
            return -math.inf
        if p == 1.0:
            return math.inf
        raise ValueError(f"norm_ppf requires 0<p<1, got {p}")
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    plow, phigh = 0.02425, 1.0 - 0.02425
    if p < plow:
        q = math.sqrt(-2.0 * math.log(p))
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    if p > phigh:
        q = math.sqrt(-2.0 * math.log(1.0 - p))
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
                ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    q = p - 0.5
    r = q * q
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / \
           (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)


# --------------------------------------------------------------------------------------
# Analytic reference: what test-retest coverage SHOULD a correctly-scaled interval show?
# --------------------------------------------------------------------------------------
def analytic_repeat_coverage(level: float, scale: float = 1.0,
                             model_to_meas_ratio: float = 0.0) -> float:
    """Predicted test--retest interval coverage for a (1-``level``)-truth-coverage interval.

    A conformal interval calibrated to cover the *truth* at level ``1-level`` (e.g. 0.90)
    spans ONE measurement-noise draw. Test--retest coverage asks whether it spans the
    difference of TWO independent draws (``Delta = eps_B - eps_A``, variance 2*sigma^2).
    In the Gaussian reference with half-width ``hw = scale * z * sqrt(sigma_m^2+sigma_mod^2)``,

        repeat_coverage = 2*Phi( scale * z * sqrt(1+r^2) / sqrt(2) ) - 1,   z = Phi^{-1}(1-level/2)

    with ``r = sigma_mod/sigma_m`` (``model_to_meas_ratio``). The key facts this encodes:

      * ``scale=1, r=0`` -> ``2*Phi(z/sqrt(2))-1`` ~ 0.755 for level=0.10: a perfectly
        measurement-scaled 90% interval is EXPECTED to show ~76% test-retest coverage,
        NOT 90% -- the derivable gap between accuracy-coverage and repeat-coverage.
      * ``r>0`` (non-repeating model/posterior spread) pushes coverage UP toward 1.0
        (interval conservative for repeatability).
      * ``scale<1`` (width too narrow) pushes coverage DOWN.
    """
    z = norm_ppf(1.0 - level / 2.0)
    arg = scale * z * math.sqrt(1.0 + model_to_meas_ratio ** 2) / math.sqrt(2.0)
    return float(2.0 * norm_cdf(arg) - 1.0)


# --------------------------------------------------------------------------------------
# The headline statistic and its scale companion
# --------------------------------------------------------------------------------------
def test_retest_coverage(est_a, est_b, lo_a, hi_a, lo_b=None, hi_b=None,
                         symmetrize: bool = True) -> float:
    """Fraction of units where one repeat's estimate falls in the other repeat's interval.

    Direction A->B: ``est_b in [lo_a, hi_a]``. If ``symmetrize`` and B's interval is
    supplied, also score B->A (``est_a in [lo_b, hi_b]``) and pool both directions.
    """
    est_a = np.asarray(est_a, float); est_b = np.asarray(est_b, float)
    lo_a = np.asarray(lo_a, float); hi_a = np.asarray(hi_a, float)
    hits = [(est_b >= lo_a) & (est_b <= hi_a)]
    if symmetrize and lo_b is not None and hi_b is not None:
        lo_b = np.asarray(lo_b, float); hi_b = np.asarray(hi_b, float)
        hits.append((est_a >= lo_b) & (est_a <= hi_b))
    return float(np.mean(np.concatenate(hits)))


def standardized_residual_dispersion(est_a, est_b, w_a, w_b, level: float,
                                     robust: bool = True) -> float:
    """Dispersion of the standardized test--retest residual ``z_i``.

    With half-width ``hw = w/2 = scale * z * sigma`` and ``Delta = est_b-est_a ~ N(0,2sigma^2)``,
    ``z_i = z_level * Delta_i / (sqrt(2) * hw_i)`` has unit dispersion iff the width is
    correctly measurement-scaled. >1 means widths too narrow; <1 means too wide.
    ``robust`` uses the MAD-based estimate (1.4826*MAD); else RMS.
    """
    est_a = np.asarray(est_a, float); est_b = np.asarray(est_b, float)
    hw = 0.5 * (np.asarray(w_a, float) + np.asarray(w_b, float)) / 2.0  # mean half-width
    z_level = norm_ppf(1.0 - level / 2.0)
    zi = z_level * (est_b - est_a) / (math.sqrt(2.0) * hw)
    zi = zi[np.isfinite(zi)]
    if robust:
        return float(1.4826 * np.median(np.abs(zi - np.median(zi))))
    return float(np.sqrt(np.mean(zi ** 2)))


def spearman(x, y) -> float:
    """Spearman rank correlation (numpy-only) -- the GAUGE-style statistic, for contrast."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if x.size < 3:
        return float("nan")
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    rx -= rx.mean(); ry -= ry.mean()
    denom = math.sqrt(float(np.sum(rx ** 2)) * float(np.sum(ry ** 2)))
    return float(np.sum(rx * ry) / denom) if denom > 0 else float("nan")


# --------------------------------------------------------------------------------------
# BCa bootstrap (numpy-only) over per-unit resampling
# --------------------------------------------------------------------------------------
def bca_ci(stat_fn, n: int, *, n_boot: int = 4000, alpha: float = 0.05,
           seed: int = 0) -> tuple[float, float, float]:
    """Bias-corrected-and-accelerated bootstrap CI for a statistic over ``n`` units.

    ``stat_fn(idx)`` takes an array of unit indices and returns a scalar. Returns
    ``(point, lo, hi)`` for the central ``1-alpha`` interval. Mirrors the BCa interval
    Gauge reports for its Spearman check, so Echo's CIs are method-comparable.
    """
    idx_all = np.arange(n)
    theta_hat = float(stat_fn(idx_all))
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot, float)
    for b in range(n_boot):
        boots[b] = stat_fn(rng.integers(0, n, size=n))
    boots = boots[np.isfinite(boots)]
    if boots.size < 10:
        return theta_hat, float("nan"), float("nan")
    # bias correction z0
    prop = float(np.mean(boots < theta_hat))
    prop = min(max(prop, 1.0 / (2 * boots.size)), 1.0 - 1.0 / (2 * boots.size))
    z0 = norm_ppf(prop)
    # acceleration via jackknife
    jack = np.array([stat_fn(np.delete(idx_all, i)) for i in range(n)], float)
    jack = jack[np.isfinite(jack)]
    jbar = jack.mean()
    num = float(np.sum((jbar - jack) ** 3))
    den = 6.0 * (float(np.sum((jbar - jack) ** 2)) ** 1.5)
    a = num / den if den != 0 else 0.0
    zl, zu = norm_ppf(alpha / 2.0), norm_ppf(1.0 - alpha / 2.0)
    a1 = norm_cdf(z0 + (z0 + zl) / (1.0 - a * (z0 + zl)))
    a2 = norm_cdf(z0 + (z0 + zu) / (1.0 - a * (z0 + zu)))
    lo = float(np.percentile(boots, 100.0 * a1))
    hi = float(np.percentile(boots, 100.0 * a2))
    return theta_hat, lo, hi


# --------------------------------------------------------------------------------------
# Bundled evaluation
# --------------------------------------------------------------------------------------
@dataclass
class CoverageResult:
    param: str
    n: int
    level: float
    coverage: float
    coverage_ci: tuple[float, float]
    z_dispersion: float
    z_dispersion_ci: tuple[float, float]
    spearman_width_vs_repeat: float
    analytic_target_meas_only: float
    extras: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        d = self.__dict__.copy()
        d["coverage_ci"] = list(self.coverage_ci)
        d["z_dispersion_ci"] = list(self.z_dispersion_ci)
        return d


def evaluate(param: str, est_a, est_b, lo_a, hi_a, lo_b, hi_b, *, level: float,
             n_boot: int = 4000, alpha: float = 0.05, seed: int = 0) -> CoverageResult:
    """Compute Echo's full repeatability scorecard for one parameter, with BCa CIs."""
    est_a = np.asarray(est_a, float); est_b = np.asarray(est_b, float)
    lo_a = np.asarray(lo_a, float); hi_a = np.asarray(hi_a, float)
    lo_b = np.asarray(lo_b, float); hi_b = np.asarray(hi_b, float)
    w_a = hi_a - lo_a; w_b = hi_b - lo_b
    n = est_a.size

    def cov_fn(idx):
        return test_retest_coverage(est_a[idx], est_b[idx], lo_a[idx], hi_a[idx],
                                    lo_b[idx], hi_b[idx], symmetrize=True)

    def disp_fn(idx):
        return standardized_residual_dispersion(est_a[idx], est_b[idx], w_a[idx], w_b[idx],
                                                 level=level)

    cov, cov_lo, cov_hi = bca_ci(cov_fn, n, n_boot=n_boot, alpha=alpha, seed=seed)
    disp, disp_lo, disp_hi = bca_ci(disp_fn, n, n_boot=n_boot, alpha=alpha, seed=seed + 1)
    # Gauge-style contrast: width vs measured scan-rescan |Delta|
    width = 0.5 * (w_a + w_b)
    rho = spearman(width, np.abs(est_b - est_a))
    return CoverageResult(
        param=param, n=int(n), level=float(level),
        coverage=cov, coverage_ci=(cov_lo, cov_hi),
        z_dispersion=disp, z_dispersion_ci=(disp_lo, disp_hi),
        spearman_width_vs_repeat=rho,
        analytic_target_meas_only=analytic_repeat_coverage(level),
    )
