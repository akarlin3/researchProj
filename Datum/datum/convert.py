r"""Convention adapter between the Gauge substrate and the Caliper estimators.

The two reused dependencies speak different IVIM conventions:

* **Gauge** (the substrate): parameter columns ``(D, D*, f)`` in *physical* units
  (D ~ 5e-4..3e-3 mm^2/s, D* ~ 1e-2..1e-1 mm^2/s), signals from
  ``f*exp(-b*D*) + (1-f)*exp(-b*D)``.
* **Caliper** (the estimators + ruler): parameter columns ``(D, f, D*)`` in units
  of *1e-3 mm^2/s* (D ~ 0.5..3, D* ~ 10..100), signals from
  ``f*exp(-b*D*·1e-3) + (1-f)*exp(-b*D·1e-3)``.

So a Gauge cohort is mapped into Caliper's convention by **reordering** the columns
and **scaling D and D\* by 1e3**; the (S0=1) magnitude signals are identical under
both forward models (verified to 0.0 max abs error -- see ``tests/test_convert.py``).
The whole benchmark then runs in one convention, and the load-bearing ruler numbers
(coverage, coverage-gap, ECE) are scale-invariant regardless.
"""
from __future__ import annotations

import numpy as np

# Column indices.
GAUGE_D, GAUGE_DSTAR, GAUGE_F = 0, 1, 2          # gauge.cohort params order (D, D*, f)
CAL_D, CAL_F, CAL_DSTAR = 0, 1, 2                 # caliper PARAM_NAMES order (D, f, D*)

_SCALE = 1.0e3                                    # physical mm^2/s -> 1e-3 mm^2/s


def gauge_to_caliper(params_gauge: np.ndarray) -> np.ndarray:
    """(N, 3) Gauge ``(D, D*, f)`` physical -> (N, 3) Caliper ``(D, f, D*)`` 1e-3."""
    p = np.asarray(params_gauge, dtype=float)
    out = np.empty_like(p)
    out[:, CAL_D] = p[:, GAUGE_D] * _SCALE
    out[:, CAL_F] = p[:, GAUGE_F]
    out[:, CAL_DSTAR] = p[:, GAUGE_DSTAR] * _SCALE
    return out


def caliper_to_gauge(params_caliper: np.ndarray) -> np.ndarray:
    """Inverse of :func:`gauge_to_caliper`."""
    p = np.asarray(params_caliper, dtype=float)
    out = np.empty_like(p)
    out[:, GAUGE_D] = p[:, CAL_D] / _SCALE
    out[:, GAUGE_DSTAR] = p[:, CAL_DSTAR] / _SCALE
    out[:, GAUGE_F] = p[:, CAL_F]
    return out
