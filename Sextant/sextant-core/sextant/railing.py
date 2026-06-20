"""Boundary-railing — the primary, assumption-free diagnostic.

A conventional bi-exponential IVIM fit is a bounded NLLS problem. When the
pseudo-diffusion coefficient D* is weakly identifiable (the IVIM perfusion
regime), the optimiser frequently terminates *on* a parameter bound rather than
at an interior minimum. The fraction of high-SNR voxels whose NLLS D* estimate
rails to a bound is a property of the optimiser + data, requiring no ground truth
and no calibration-trust argument — it cannot be "overextended".

All fitting and thresholds are reused read-only from Fashion
(:mod:`sextant.fashion_reuse`); this module only orchestrates and characterises.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .fashion_reuse import load_railing, load_wide

# Default SNR strata for the "when does it rail" characterisation.
SNR_STRATA = [(8.0, 15.0), (15.0, 30.0), (30.0, 60.0), (60.0, np.inf)]


def fit_dstar(fit_signals: np.ndarray, wide: bool = False) -> np.ndarray:
    """Return per-voxel NLLS D* estimates (mm^2/s) for the 10-pt TARGET_BVALS.

    ``wide=False`` uses Fashion's tight prior-box bounds (the original analysis);
    ``wide=True`` uses the generous wide bounds (sensitivity check that railing is
    not merely an artefact of tight bounds).
    """
    R = load_railing()
    b = R["TARGET_BVALS"]
    fit_signals = np.asarray(fit_signals, float)
    if wide:
        W = load_wide()
        fit = W["fit_biexp_wide"]
    else:
        fit = R["fit_biexp_nlls"]
    return np.array([fit(b, fit_signals[i])[1] for i in range(len(fit_signals))])


def rail_mask(dstar: np.ndarray, lower: float, upper: float):
    """Boolean rail mask plus separate lower/upper rail masks (NaNs -> not railed)."""
    d = np.asarray(dstar, float)
    finite = np.isfinite(d)
    lo = finite & (d <= lower)
    hi = finite & (d >= upper)
    return (lo | hi), lo, hi


@dataclass
class RailingResult:
    cohort: str
    n_voxels: int
    n_high_snr: int
    snr_floor: float
    bounds: str               # "tight" or "wide"
    lower_rail: float
    upper_rail: float
    frac_railed: float
    frac_lower: float
    frac_upper: float
    n_railed: int
    by_snr: list = field(default_factory=list)   # [{lo,hi,n,frac_railed,frac_lower,frac_upper}]

    def to_dict(self) -> dict:
        return {k: getattr(self, k) for k in self.__dataclass_fields__}


def _strata(snr_hi, railed, lower, upper, dstar_hi):
    out = []
    for lo, hi in SNR_STRATA:
        sel = (snr_hi >= lo) & (snr_hi < hi)
        n = int(sel.sum())
        if n == 0:
            out.append({"lo": lo, "hi": (None if np.isinf(hi) else hi),
                        "n": 0, "frac_railed": None, "frac_lower": None, "frac_upper": None})
            continue
        d = dstar_hi[sel]
        m, ml, mu = rail_mask(d, lower, upper)
        out.append({
            "lo": lo, "hi": (None if np.isinf(hi) else hi), "n": n,
            "frac_railed": float(np.mean(m)),
            "frac_lower": float(np.mean(ml)), "frac_upper": float(np.mean(mu)),
        })
    return out


def analyze_cohort(cohort, *, bounds: str = "tight") -> RailingResult:
    """Compute the railing diagnostic for one cohort (see :class:`RailingResult`)."""
    R = load_railing()
    snr_floor = float(R["SNR_FLOOR"])
    snrs = np.asarray(cohort.snrs, float)
    roi_hi = snrs >= snr_floor

    wide = bounds == "wide"
    dstar = fit_dstar(cohort.fit_signals, wide=wide)
    dstar_hi = dstar[roi_hi]
    snr_hi = snrs[roi_hi]

    if wide:
        W = load_wide()
        lower = float(W["WIDE_LOW"][1] + W["_RAIL_TOL"] * (W["WIDE_HIGH"][1] - W["WIDE_LOW"][1]))
        upper = float(W["WIDE_HIGH"][1] - W["_RAIL_TOL"] * (W["WIDE_HIGH"][1] - W["WIDE_LOW"][1]))
    else:
        lower = float(R["DSTAR_LOWER_RAIL"])
        upper = float(R["DSTAR_UPPER_RAIL"])

    m, ml, mu = rail_mask(dstar_hi, lower, upper)
    res = RailingResult(
        cohort=cohort.name, n_voxels=int(cohort.n_voxels), n_high_snr=int(roi_hi.sum()),
        snr_floor=snr_floor, bounds=bounds, lower_rail=lower, upper_rail=upper,
        frac_railed=float(np.mean(m)), frac_lower=float(np.mean(ml)),
        frac_upper=float(np.mean(mu)), n_railed=int(np.sum(m)),
        by_snr=_strata(snr_hi, m, lower, upper, dstar_hi),
    )
    # Per-voxel arrays for downstream bootstrap (not serialised by to_dict()).
    res.railed_hi = m
    res.dstar_hi = dstar_hi
    res.snr_hi = snr_hi
    return res
