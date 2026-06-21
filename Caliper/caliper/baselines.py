"""caliper.baselines -- a box-constrained NLLS bi-exponential IVIM baseline.

The classical reference estimator: a full four-parameter ``(S0, D, f, D*)``
non-linear least-squares fit of the bi-exponential IVIM model
(``caliper.forward.ivim_signal``) to each voxel's multi-b decay, with the
optimiser confined to a physiological box via ``scipy.optimize.least_squares``
(trust-region reflective). Uncertainty is the textbook asymptotic Gaussian:
``cov = sigma^2 * (J^T J)^-1`` from the Jacobian at the solution.

This module exists to reproduce -- on *synthetic phantoms only* -- the two
findings the (retooled) Fashion manuscript builds on for constrained NLLS, in
their retooled priority order:

1. **Boundary railing (the assumption-free primary).** The pseudo-diffusion
   ``D*`` is weakly identified from a sparse, noisy acquisition, so the
   constrained optimiser frequently pins it against a box bound.
   ``boundary_railing_rate`` measures the fraction of voxels whose ``D*``
   estimate lands on a bound -- a per-voxel identifiability signature that needs
   no ground truth. (The manuscript's real-data rates -- ~55% on the open OSIPI
   abdomen, replicated 47.8-73.4% across cohorts including an independent liver
   -- live in the paper and are NOT reproduced here.)
2. **Conditional interval under-coverage (a scoped, ground-truth-only readout).**
   The asymptotic-Gaussian error bars are reported under the **honest CRLB**
   convention (wide where ``D*`` is unidentified; see ``_asymptotic_sigma`` and
   ``sd_convention``). Scored by the calibration ruler (``caliper.metrics``) on
   synthetic ground truth, they under-cover ``D*`` *conditionally* in the
   high-``D*`` tercile -- not the dramatic marginal severity an earlier framing
   reported. The ruler is a scoped secondary instrument: it needs ground truth
   and cannot be applied to the real scan.

SD convention (the reviewer-flagged choice; see Gnomon ``docs/METHODS.md`` §5b).
A railed/unidentified ``D*`` carries no local Fisher information, so its SD is a
modelling *choice*. ``sd_convention="honest"`` (default) widens it -- the
statistically honest admission that ``D*`` is undetermined; ``"floored"`` is an
illustrative reconstruction that instead reports a narrow floored SD
("overconfident by design") and is shown only to explain how the now-dropped
marginal 0.30 / 0.67 D* coverage severity arose. Default and recommended: honest.

NOTE (in review): the Fashion manuscript is under peer review at **NMR in
Biomedicine**. This module reproduces only the *phenomenon* on in-repo synthetic
data; the clinical/real-data percentages live in the paper and are deliberately
NOT reproduced here. Keep this module private until the paper clears.

scipy + numpy. The fit is a per-voxel optimisation loop (no torch).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from statistics import NormalDist

import numpy as np
from scipy.optimize import least_squares

from .forward import _SCALE, DEFAULT_BVALUES, PARAM_NAMES

__all__ = ["NLLSFit", "NLLSIVIMEstimator"]

# Order of the four fitted parameters inside the optimiser.
_FIT_NAMES = ("S0", "D", "f", "Dstar")
# Indices of the three reported IVIM parameters (D, f, D*) within the 4-vector.
_REPORT_IDX = np.array([1, 2, 3])
# D* lives at fit-vector index 3; it is the weakly-identified, railing parameter.
_DSTAR_FIT_IDX = 3


def _ivim_residuals(theta: np.ndarray, b: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Residual ``model(theta) - y`` for one voxel; theta = (S0, D, f, D*)."""
    s0, D, f, Ds = theta
    model = s0 * (f * np.exp(-b * Ds * _SCALE) + (1.0 - f) * np.exp(-b * D * _SCALE))
    return model - y


@dataclass
class NLLSFit:
    """Result of fitting a cohort with :class:`NLLSIVIMEstimator`.

    Attributes
    ----------
    params : (n, 3) reported IVIM estimates, columns ``(D, f, D*)``.
    s0 : (n,) fitted baseline amplitude.
    sigma : (n, 3) asymptotic Gaussian SD for ``(D, f, D*)``.
    railed : (n, 3) bool; True where the estimate sits on a box bound.
    success : (n,) bool; optimiser convergence flag.
    """

    params: np.ndarray
    s0: np.ndarray
    sigma: np.ndarray
    railed: np.ndarray
    success: np.ndarray

    @property
    def dstar_railed(self) -> np.ndarray:
        """Per-voxel boolean: is the ``D*`` estimate railed on a bound?"""
        return self.railed[:, 2]

    @property
    def any_railed(self) -> np.ndarray:
        """Per-voxel boolean: is *any* reported parameter railed on a bound?"""
        return np.any(self.railed, axis=1)


@dataclass
class NLLSIVIMEstimator:
    """Box-constrained four-parameter NLLS bi-exponential IVIM fitter.

    Exposes the same prediction contract as the other Caliper estimators::

        predict_quantiles(signals, q_levels) -> ndarray (n, 3, n_levels)

    so it can be scored by ``caliper.metrics`` side by side with the MAF flow
    posterior. The reported quantiles are the asymptotic-Gaussian model intervals
    ``point + sigma * z(level)`` -- deliberately *not* re-calibrated; their
    over-confidence is the thing under test.

    Parameters
    ----------
    bvalues : the acquisition b-schedule the signals were sampled on.
    lower, upper : the box bounds for ``(S0, D, f, D*)``. Defaults are standard
        physiological IVIM bounds (D, D* in 1e-3 mm^2/s); ``D*`` in [3, 300]
        brackets the synthetic prior [10, 100] with the head-room a clinical fit
        would allow, leaving room for the weakly-identified fits to rail.
    init : the initial guess ``(D, f, D*)`` for every voxel; ``S0`` is initialised
        per-voxel from the b=0 signal. A shared, deliberately neutral start --
        railing is driven by identifiability, not by a peeked-at warm start.
    rail_tol : a parameter counts as railed when it lies within ``rail_tol`` of a
        bound, measured as a fraction of that parameter's bound span.
    sigma_floor : a small positive floor on the reported SD (guards the degenerate
        Jacobian at a railed solution from producing a zero-width interval).
    sd_convention : how a railed/unidentified ``D*``'s SD is reported. ``"honest"``
        (default) keeps the wide, span-clipped asymptotic SD -- the statistically
        honest admission that ``D*`` is undetermined where it rails. ``"floored"``
        is an illustrative reconstruction that overwrites a railed ``D*``'s SD with
        the narrow ``railed_sd_floor`` (the "overconfident by design" convention
        that manufactures the now-dropped marginal severity). See Gnomon
        ``docs/METHODS.md`` §5b. Default and recommended: honest.
    railed_sd_floor : the floored-convention SD for a railed ``D*``, in Caliper's
        ``D*`` units (Ds, where 1 unit == 1e-3 mm^2/s); default 3.0 == 0.003
        mm^2/s, matching Gnomon's ``RAILED_SD_FLOOR``. Ignored under the honest
        convention.
    max_nfev : optimiser iteration budget per voxel.
    """

    bvalues: np.ndarray = field(default_factory=lambda: DEFAULT_BVALUES.copy())
    lower: tuple[float, float, float, float] = (0.1, 0.1, 0.0, 3.0)
    upper: tuple[float, float, float, float] = (2.0, 4.0, 1.0, 300.0)
    init: tuple[float, float, float] = (1.0, 0.1, 20.0)
    rail_tol: float = 1e-3
    sigma_floor: float = 1e-6
    sd_convention: str = "honest"
    railed_sd_floor: float = 3.0
    max_nfev: int = 400
    eps: float = 1e-4

    param_names: tuple[str, ...] = PARAM_NAMES

    def __post_init__(self) -> None:
        self.bvalues = np.asarray(self.bvalues, dtype=float)
        self._lo = np.asarray(self.lower, dtype=float)
        self._hi = np.asarray(self.upper, dtype=float)
        if self._lo.shape != (4,) or self._hi.shape != (4,):
            raise ValueError("lower/upper must each have 4 entries (S0, D, f, D*)")
        if np.any(self._hi <= self._lo):
            raise ValueError("each upper bound must exceed its lower bound")
        if self.sd_convention not in ("honest", "floored"):
            raise ValueError("sd_convention must be 'honest' or 'floored'")
        if self.railed_sd_floor <= 0.0:
            raise ValueError("railed_sd_floor must be positive")
        self._span = self._hi - self._lo
        self._b0_idx = int(np.argmin(self.bvalues))

    # API symmetry with the other estimators -- NLLS needs no training.
    def fit(self, signals=None, params=None) -> "NLLSIVIMEstimator":
        """No-op (NLLS fits per voxel at predict time); returns ``self``."""
        return self

    # ----- core per-cohort solve ----------------------------------------- #
    def solve(self, signals: np.ndarray) -> NLLSFit:
        """Fit every voxel; return point estimates, SDs and railing flags.

        ``signals`` is ``(n, n_bvalues)`` matching the b-schedule.
        """
        S = np.asarray(signals, dtype=float)
        if S.ndim != 2 or S.shape[1] != self.bvalues.shape[0]:
            raise ValueError("signals must be (n, n_bvalues) matching the schedule")
        n = S.shape[0]
        b = self.bvalues

        point4 = np.empty((n, 4))
        sigma4 = np.empty((n, 4))
        success = np.empty(n, dtype=bool)
        for i in range(n):
            y = S[i]
            s0_init = float(np.clip(y[self._b0_idx], self.eps, None))
            x0 = np.array([s0_init, self.init[0], self.init[1], self.init[2]])
            x0 = np.clip(x0, self._lo, self._hi)
            res = least_squares(
                _ivim_residuals, x0, bounds=(self._lo, self._hi),
                args=(b, y), method="trf", max_nfev=self.max_nfev,
            )
            point4[i] = res.x
            success[i] = bool(res.success)
            sigma4[i] = self._asymptotic_sigma(res)

        railed4 = (point4 - self._lo[None, :] <= self.rail_tol * self._span[None, :]) | (
            self._hi[None, :] - point4 <= self.rail_tol * self._span[None, :]
        )
        if self.sd_convention == "floored":
            # Illustrative "overconfident floor" convention (Fashion-implied): a
            # railed/unidentified D* is reported with a *narrow* floored SD rather
            # than the honest wide one, so its interval under-covers. This is the
            # choice that reconstructs the now-dropped marginal D* severity; it is
            # shown only as a labelled illustration (see Gnomon docs/METHODS.md
            # §5b) and is never the default. We floor the railed D* SD only.
            dstar_railed = railed4[:, _DSTAR_FIT_IDX]
            sigma4[dstar_railed, _DSTAR_FIT_IDX] = self.railed_sd_floor
        return NLLSFit(
            params=point4[:, _REPORT_IDX],
            s0=point4[:, 0],
            sigma=sigma4[:, _REPORT_IDX],
            railed=railed4[:, _REPORT_IDX],
            success=success,
        )

    def _asymptotic_sigma(self, res) -> np.ndarray:
        """Textbook NLLS SD: sqrt(diag(sigma^2 (J^T J)^-1)), per parameter.

        The variance is clipped to each parameter's bound span: on a weakly
        identified (near-singular Jacobian) fit the raw asymptotic variance can
        diverge to values far exceeding the feasible range, which no real fit
        pipeline would report as an error bar. Clipping to the span keeps the SD
        physical without hand-tuning, and never tightens an interval.
        """
        J = np.asarray(res.jac, dtype=float)  # (m, 4)
        m, p = J.shape
        dof = max(m - p, 1)
        chi2 = float(np.sum(res.fun ** 2))
        resid_var = chi2 / dof
        JtJ = J.T @ J
        cov = resid_var * np.linalg.pinv(JtJ)
        var = np.clip(np.diag(cov), 0.0, self._span ** 2)
        return np.sqrt(var) + self.sigma_floor

    # ----- prediction contract ------------------------------------------- #
    def predict_point(self, signals: np.ndarray) -> np.ndarray:
        """Box-constrained NLLS point estimate. Returns (n, 3) = (D, f, D*)."""
        return self.solve(signals).params

    def predict_quantiles(self, signals: np.ndarray, q_levels) -> np.ndarray:
        """Asymptotic-Gaussian model quantiles. Returns (n, 3, n_levels).

        ``q = point + sigma * z(level)`` -- monotone in level and symmetric about
        the point estimate by construction. These intervals are *not* conformally
        corrected: their over-confidence is exactly the failure mode under test.
        """
        q_levels = np.asarray(q_levels, dtype=float)
        if np.any((q_levels <= 0.0) | (q_levels >= 1.0)):
            raise ValueError("q_levels must lie strictly in (0, 1)")
        fit = self.solve(signals)
        z = np.array([NormalDist().inv_cdf(float(p)) for p in q_levels])  # (L,)
        q = fit.params[:, :, None] + fit.sigma[:, :, None] * z[None, None, :]
        return q

    # ----- railing diagnostics ------------------------------------------- #
    def boundary_railing_rate(self, signals: np.ndarray, param: str = "Dstar") -> float:
        """Fraction of voxels whose ``param`` estimate is railed on a box bound.

        ``param`` is one of ``("D", "f", "Dstar")`` or ``"any"``. Default
        ``"Dstar"`` -- the weakly-identified parameter the manuscript highlights.
        """
        return self.railing_rate_from_fit(self.solve(signals), param)

    @staticmethod
    def railing_rate_from_fit(fit: NLLSFit, param: str = "Dstar") -> float:
        """Railing rate from an already-computed :class:`NLLSFit` (avoids refit)."""
        if param == "any":
            return float(np.mean(fit.any_railed))
        try:
            j = PARAM_NAMES.index(param)
        except ValueError as exc:
            raise ValueError(
                f"param must be one of {PARAM_NAMES + ('any',)}, got {param!r}"
            ) from exc
        return float(np.mean(fit.railed[:, j]))


if __name__ == "__main__":
    from caliper import metrics as M
    from caliper.forward import synthetic_cohort

    levels = np.array([0.05, 0.25, 0.5, 0.75, 0.95])
    test = synthetic_cohort(n=400, snr=20.0, seed=2)
    est = NLLSIVIMEstimator()
    fit = est.solve(test.signals)
    q = est.predict_quantiles(test.signals, levels)
    scores = M.score_quantiles(test.params, q, levels, alpha=0.1,
                               param_names=PARAM_NAMES, conditioning=test.params)
    print(M.format_scorecard(scores, title="constrained NLLS baseline (raw)"))
    print()
    for p in ("D", "f", "Dstar"):
        print(f"boundary-railing rate [{p:>5}] = "
              f"{est.railing_rate_from_fit(fit, p):.3f}")
