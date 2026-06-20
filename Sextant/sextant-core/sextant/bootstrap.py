"""Voxel-level bootstrap confidence intervals for the railing fraction.

The railing fraction is a mean of a per-voxel boolean; its sampling uncertainty
over the (single-subject) voxel population is quantified by resampling voxels
with replacement. Seeded (``GLOBAL_SEED``) and reproducible.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .seeding import GLOBAL_SEED, make_rng

DEFAULT_B = 5000


@dataclass
class BootCI:
    point: float
    lo: float
    hi: float
    level: float
    n_boot: int
    seed: int

    def to_dict(self) -> dict:
        return {k: getattr(self, k) for k in self.__dataclass_fields__}


def bootstrap_fraction(railed_bool, *, n_boot: int = DEFAULT_B, level: float = 0.95,
                       seed: int = GLOBAL_SEED) -> BootCI:
    """Percentile bootstrap CI for the mean of a boolean per-voxel array."""
    x = np.asarray(railed_bool, bool)
    n = x.size
    if n == 0:
        raise ValueError("empty railed array")
    rng = make_rng(seed)
    point = float(x.mean())
    # vectorised resample: draw indices (n_boot, n)
    idx = rng.integers(0, n, size=(n_boot, n))
    means = x[idx].mean(axis=1)
    a = (1.0 - level) / 2.0
    lo, hi = np.percentile(means, [100 * a, 100 * (1 - a)])
    return BootCI(point=point, lo=float(lo), hi=float(hi),
                  level=level, n_boot=n_boot, seed=seed)
