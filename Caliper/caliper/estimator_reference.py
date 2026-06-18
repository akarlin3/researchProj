"""caliper.estimator_reference -- a numpy-only, deliberately over-confident IVIM
reference estimator.

This is the *device under test* for the conformal wrapper. It stands in for the
torch MAF (``caliper.estimator_maf``) so the conformal layer and its tests run
without torch, while exposing the identical contract::

    predict_quantiles(signals, q_levels) -> ndarray (n, n_params, n_levels)

The point estimate is a genuine **segmented least-squares IVIM fit** (the classic
two-step method): a high-b log-linear fit recovers the tissue diffusivity ``D``
and amplitude, from which the perfusion fraction ``f`` follows; a low-b
perfusion-residual log-linear fit recovers the pseudo-diffusion ``D*``. ``D`` and
``f`` come out close to truth; ``D*`` is poorly identified -- and *systematically*
worse at high ``D*``, where the perfusion term has decayed away by the second or
third b-value. That heteroscedastic-in-D* error is not a bug: it is the
identifiability wall this toolkit exists to measure.

The *reported* predictive quantiles, by contrast, are homoscedastic and
deliberately **too narrow** (a fixed believed sigma per parameter). The mismatch
between honest (heteroscedastic, large) error and the over-confident (constant,
small) reported spread is what makes the raw estimator miscalibrated -- exactly
the model-based-UQ failure mode that conformal prediction is meant to repair.

Numpy only. No torch, no scipy, no fitting loop -- the segmented fit is closed
form and vectorised over voxels.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from statistics import NormalDist

import numpy as np

from .forward import DEFAULT_BVALUES, PARAM_NAMES

__all__ = ["ReferenceIVIMEstimator"]


def _loglinfit(x: np.ndarray, y: np.ndarray):
    """Vectorised ordinary least squares of ``y`` (n, k) on a shared ``x`` (k,).

    Returns ``(slope, intercept)`` each shape (n,). ``x`` is identical for every
    row, so the normal equations reduce to closed-form moments.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    xbar = x.mean()
    xc = x - xbar
    denom = float(np.sum(xc * xc))
    ybar = y.mean(axis=1)
    slope = (y - ybar[:, None]) @ xc / denom
    intercept = ybar - slope * xbar
    return slope, intercept


@dataclass
class ReferenceIVIMEstimator:
    """Over-confident segmented-fit IVIM estimator (numpy only).

    Parameters
    ----------
    bvalues : the acquisition b-schedule the signals were sampled on.
    sigma_D, sigma_f, sigma_Dstar : the (deliberately narrow) believed standard
        deviations the estimator reports for each parameter. These are constant
        across voxels -- the estimator advertises the same precision everywhere,
        which is precisely why it is over-confident where the true error is large.
    high_b_min : b-values at or above this define the tissue (D) fit window.
    low_b_max : positive b-values at or below this define the perfusion (D*) fit
        window.
    eps : positive floor applied before taking logs (guards noise-driven
        non-positive signals / perfusion residuals).

    The defaults are tuned so the raw marginal coverage of a 90% interval lands
    well below nominal (~0.5-0.65) on the standard synthetic cohort -- comparable
    to the raw MAF -- without peeking at the conformal result.
    """

    bvalues: np.ndarray = field(default_factory=lambda: DEFAULT_BVALUES.copy())
    sigma_D: float = 0.11
    sigma_f: float = 0.022
    sigma_Dstar: float = 6.0
    high_b_min: float = 200.0
    low_b_max: float = 50.0
    eps: float = 1e-4

    param_names: tuple[str, ...] = PARAM_NAMES

    def __post_init__(self) -> None:
        self.bvalues = np.asarray(self.bvalues, dtype=float)
        self._high = self.bvalues >= self.high_b_min
        self._low = (self.bvalues <= self.low_b_max) & (self.bvalues > 0.0)
        if self._high.sum() < 2:
            raise ValueError("need >=2 b-values at/above high_b_min for the D fit")
        if self._low.sum() < 2:
            raise ValueError("need >=2 positive b-values at/below low_b_max for D*")
        self._reported_sigma = np.array(
            [self.sigma_D, self.sigma_f, self.sigma_Dstar], dtype=float
        )

    # API symmetry with MAFPosterior -- the segmented fit needs no training.
    def fit(self, signals=None, params=None) -> "ReferenceIVIMEstimator":
        """No-op (the segmented fit is closed-form); returns ``self`` for API parity."""
        return self

    def predict_point(self, signals: np.ndarray) -> np.ndarray:
        """Segmented least-squares IVIM fit. Returns (n, 3) = (D, f, D*)."""
        S = np.asarray(signals, dtype=float)
        if S.ndim != 2 or S.shape[1] != self.bvalues.shape[0]:
            raise ValueError("signals must be (n, n_bvalues) matching the schedule")

        # S0 from the lowest-b acquisition (b=0 if present, else the smallest b).
        s0 = np.maximum(S[:, int(np.argmin(self.bvalues))], self.eps)

        # --- step 1: high-b tissue fit -> D and amplitude A = (1-f) * S0 ------
        b_hi = self.bvalues[self._high]
        y_hi = np.log(np.maximum(S[:, self._high], self.eps))
        slope_hi, intercept_hi = _loglinfit(b_hi, y_hi)
        D = np.clip(-slope_hi / 1e-3, 0.1, 4.0)
        A = np.clip(np.exp(intercept_hi), self.eps, None)

        # --- step 2: perfusion fraction f = 1 - A / S0 -----------------------
        f = np.clip(1.0 - A / s0, 0.001, 0.95)

        # --- step 3: low-b perfusion-residual fit -> D* ----------------------
        b_lo = self.bvalues[self._low]
        tissue_lo = A[:, None] * np.exp(-b_lo[None, :] * D[:, None] * 1e-3)
        perf = S[:, self._low] - tissue_lo
        y_lo = np.log(np.maximum(perf, self.eps))
        slope_lo, _ = _loglinfit(b_lo, y_lo)
        Dstar = np.clip(-slope_lo / 1e-3, 1.0, 300.0)

        return np.stack([D, f, Dstar], axis=1)

    def predict_quantiles(self, signals: np.ndarray, q_levels) -> np.ndarray:
        """Over-confident Gaussian quantiles about the segmented-fit point.

        Returns (n, n_params, n_levels). Quantiles are point + sigma_p * z(level),
        with z the standard-normal inverse CDF, so they are monotone in level and
        symmetric about the median by construction.
        """
        q_levels = np.asarray(q_levels, dtype=float)
        if np.any((q_levels <= 0.0) | (q_levels >= 1.0)):
            raise ValueError("q_levels must lie strictly in (0, 1)")
        point = self.predict_point(signals)  # (n, P)
        z = np.array([NormalDist().inv_cdf(float(p)) for p in q_levels])  # (L,)
        q = point[:, :, None] + self._reported_sigma[None, :, None] * z[None, None, :]
        return q


if __name__ == "__main__":
    from caliper import metrics as M
    from caliper.forward import synthetic_cohort

    levels = np.array([0.05, 0.25, 0.5, 0.75, 0.95])
    test = synthetic_cohort(n=4000, snr=40.0, seed=2)
    est = ReferenceIVIMEstimator()
    q = est.predict_quantiles(test.signals, levels)
    scores = M.score_quantiles(test.params, q, levels, alpha=0.1,
                               param_names=PARAM_NAMES, conditioning=test.params)
    print(M.format_scorecard(scores, title="reference estimator (raw, over-confident)"))
