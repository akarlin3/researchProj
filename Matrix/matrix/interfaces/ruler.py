"""Ruler interface (consumed component: **Fashion**, the calibrated error bar).

Contract
--------
``calibrate(mu, raw_sigma, truth=None) -> dict`` turns a *raw* per-voxel posterior
spread into a *calibrated* one, plus a coverage/ECE readout::

    {"sigma": {param: calibrated_sigma}, "interval": {param: (lo, hi)},
     "coverage": {param: {level: empirical_cov}}, "ece": {param: float}}

The real Fashion ruler (``uq.calib.coverage`` / ``uq.calib.ece`` / the skew-aware
``uq.bayesian.mcmc_uncertainty`` posterior, in review @ *NMR in Biomedicine*) drops in
behind this same signature **without touching loop.py**. See ``ASSUMPTIONS.md`` §Fashion.

Implementations here:
  * ``PassThroughRuler``  — CP1 stub: identity (calibrated == raw). Not a ruler.
  * ``PlaceholderRuler``  — CP2 placeholder: a transparent, clearly-labelled NOT-Fashion
                            calibrating ruler. Provisional.
"""
from __future__ import annotations

import numpy as np

LEVELS = (0.50, 0.80, 0.95)
_Z = {0.50: 0.674, 0.80: 1.282, 0.95: 1.960}


class PassThroughRuler:
    """CP1 stub — returns the raw spread unchanged. NOT Fashion."""

    name = "passthrough-ruler"
    label = "NOT-Fashion (CP1 identity stub)"
    provisional = True

    def calibrate(self, mu, raw_sigma, truth=None):
        sigma = {k: np.asarray(v, float).copy() for k, v in raw_sigma.items()}
        interval = {k: (np.asarray(mu[k]) - 1.96 * sigma[k],
                        np.asarray(mu[k]) + 1.96 * sigma[k]) for k in sigma}
        return dict(sigma=sigma, interval=interval, coverage={}, ece={})


class PlaceholderRuler:
    """CP2 placeholder calibrating ruler — clearly labelled, **NOT** Fashion.

    A transparent stand-in for Fashion's calibration ruler: it rescales the raw spread
    by a per-parameter factor so nominal intervals achieve nominal empirical coverage on
    the twin, and reports coverage + ECE the way Fashion's ``uq.calib`` does. The real
    Fashion ruler replaces this object; the loop is unchanged.

    ``cal_factor`` is the only "knowledge" injected; with ``truth`` available (synthetic
    twin) it is estimated empirically, else it falls back to a documented constant.
    """

    name = "placeholder-ruler"
    label = "NOT-Fashion (CP2 transparent placeholder; provisional pending Fashion @ NMRB)"
    provisional = True

    # documented fallback rescale per param when no truth is available
    FALLBACK = dict(D=1.0, Dstar=1.6, f=1.4)

    def calibrate(self, mu, raw_sigma, truth=None):
        sigma, interval, coverage, ece = {}, {}, {}, {}
        for k in raw_sigma:
            m = np.asarray(mu[k], float)
            s_raw = np.maximum(np.asarray(raw_sigma[k], float), 1e-9)
            if truth is not None and k in truth:
                t = np.asarray(truth[k], float)
                # factor that makes the standardised residual unit-variance
                z = (m - t) / s_raw
                factor = float(np.sqrt(np.mean(z * z)))
                factor = min(max(factor, 0.25), 8.0)
            else:
                factor = self.FALLBACK.get(k, 1.0)
            s_cal = s_raw * factor
            sigma[k] = s_cal
            interval[k] = (m - 1.96 * s_cal, m + 1.96 * s_cal)
            if truth is not None and k in truth:
                t = np.asarray(truth[k], float)
                cov = {lv: float(np.mean(np.abs(m - t) <= _Z[lv] * s_cal)) for lv in LEVELS}
                coverage[k] = cov
                ece[k] = float(np.mean([abs(cov[lv] - lv) for lv in LEVELS]))
        return dict(sigma=sigma, interval=interval, coverage=coverage, ece=ece)
