"""Divergence diagnostic for the CP3 T3a/T3b miss (reproducible).

CP3 found that Gnomon reproduces Fashion's NLLS railing rate (T1), the quantile
fix (T3c), and the flow-vs-railed-NLLS behavior (T4), but NOT the *severe* marginal
Gaussian under-coverage (T3a 0.30, T3b 0.67) -- Gnomon's clean rebuild gives milder
0.80 / 0.90. This script isolates the two causes, so the divergence report's "likely
cause" is run, not asserted:

  1. REGIME: the under-coverage is concentrated in the high-D* tercile (the Gauge
     wall). Pooling across a prior-spanning cohort dilutes it; concentrating in the
     hard corner sharpens it.
  2. RAILED-VOXEL UNCERTAINTY CONVENTION: Fashion's baseline is "overconfident by
     design" -- it must assign near-zero SD to railed/unidentified D*. Gnomon's
     *honest* CRLB assigns those voxels WIDE intervals (the statistically correct
     behavior). Swapping to the overconfident convention drops coverage toward
     Fashion's numbers.

Run: ``KMP_DUPLICATE_LIB_OK=TRUE PYTHONPATH=. python scripts/divergence_diagnostic.py``
"""
from __future__ import annotations

import numpy as np
from scipy.stats import norm

from gnomon import cohort, nlls, metrics as M


def laplace_dstar_coverage(co, fit, sd_dstar):
    """D* central-0.95 coverage (pooled + per-true-D*-tercile) for a given D* SD."""
    z = norm.ppf([0.025, 0.975])
    mean = fit.params[:, 1]
    lo, hi = mean + sd_dstar * z[0], mean + sd_dstar * z[1]
    ind = (co.params_true[:, 1] >= lo) & (co.params_true[:, 1] <= hi)
    g = M.tercile_groups(co.params_true[:, 1])
    terc = {int(k): round(float(ind[g == k].mean()), 3) for k in (0, 1, 2)}
    return round(float(ind.mean()), 3), terc


def main():
    co = cohort.make_headline_cohort(n_noise=200)
    est = nlls.NLLSEstimator(co.bvalues)
    fit = est.fit(co.signals, sigma=co.sigma)         # per-voxel known sigma
    sd_honest = fit.sigma[:, 1]                         # honest capped CRLB (default)
    # Overconfident-by-design: railed/unidentified D* gets a small floored SD.
    sd_oc = np.where(fit.dstar_railed, 0.003, sd_honest)

    print(f"synthetic D* railing rate: {fit.dstar_railed.mean():.3f}")
    p, t = laplace_dstar_coverage(co, fit, sd_honest)
    print(f"HONEST capped CRLB   -> pooled {p}, terciles(low/mid/high D*) {t}")
    p, t = laplace_dstar_coverage(co, fit, sd_oc)
    print(f"OVERCONFIDENT floor  -> pooled {p}, terciles(low/mid/high D*) {t}")
    print("\nInterpretation: Fashion's 0.30/0.67 require BOTH a high-D* cohort AND the "
          "overconfident floored-SD convention; the honest CRLB + prior-spanning "
          "cohort gives milder marginal coverage while the mechanism still holds.")


if __name__ == "__main__":
    main()
