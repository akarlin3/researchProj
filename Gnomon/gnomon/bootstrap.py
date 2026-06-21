"""Bootstrap confidence intervals for load-bearing numbers.

Guardrail: bootstrap CIs on every load-bearing number (railing rate, coverages,
ECE/sharpness gaps). Clean-room, numpy-only resampling, seeded from
``manifest.BOOTSTRAP``. The unit of resampling is the per-observation contribution:
for a coverage, pass the in/out indicator; for a paired gap, pass the per-unit
difference. CP3 uses these to decide PASS/FAIL (claimed value inside the CI, or point
within frozen tolerance; directional gaps require the CI to exclude 0).
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.stats import norm


@dataclass
class CI:
    point: float
    lo: float
    hi: float
    method: str

    def excludes_zero(self):
        return self.lo > 0.0 or self.hi < 0.0

    def contains(self, value):
        return self.lo <= value <= self.hi


def bootstrap_ci(values, n_boot=2000, ci=0.95, seed=0, method="percentile"):
    """Bootstrap CI for the **mean** of ``values`` (per-observation contributions)."""
    x = np.asarray(values, dtype=float)
    n = len(x)
    point = float(np.mean(x))
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(n_boot, n))
    boot = x[idx].mean(axis=1)
    lo_q, hi_q = (1 - ci) / 2, 1 - (1 - ci) / 2
    if method == "percentile":
        lo, hi = np.quantile(boot, [lo_q, hi_q])
    elif method == "bca":
        lo, hi = _bca(x, boot, point, lo_q, hi_q)
    else:
        raise ValueError(f"unknown method {method!r}")
    return CI(point=point, lo=float(lo), hi=float(hi), method=method)


def _bca(x, boot, point, lo_q, hi_q):
    # Bias-correction and acceleration (DiCiccio & Efron 1996) for the mean.
    z0 = norm.ppf(np.mean(boot < point)) if 0 < np.mean(boot < point) < 1 else 0.0
    n = len(x)
    jack = (np.sum(x) - x) / (n - 1)               # leave-one-out means
    jbar = jack.mean()
    num = np.sum((jbar - jack) ** 3)
    den = 6.0 * (np.sum((jbar - jack) ** 2) ** 1.5)
    a = num / den if den != 0 else 0.0
    zl, zh = norm.ppf(lo_q), norm.ppf(hi_q)
    adj = lambda z: norm.cdf(z0 + (z0 + z) / (1 - a * (z0 + z)))
    return np.quantile(boot, [adj(zl), adj(zh)])
