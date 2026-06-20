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
                       seed: int = GLOBAL_SEED, chunk: int = 200) -> BootCI:
    """Percentile bootstrap CI for the mean of a boolean per-voxel array.

    Resamples are drawn in chunks so peak memory is ``O(chunk * n)`` rather than
    ``O(n_boot * n)`` (the latter is multiple GB for 50k voxels x 5000 resamples).
    """
    x = np.asarray(railed_bool, bool).astype(np.float64)
    n = x.size
    if n == 0:
        raise ValueError("empty railed array")
    rng = make_rng(seed)
    point = float(x.mean())
    means = np.empty(n_boot, float)
    done = 0
    while done < n_boot:
        m = min(chunk, n_boot - done)
        idx = rng.integers(0, n, size=(m, n))
        means[done:done + m] = x[idx].mean(axis=1)
        done += m
    a = (1.0 - level) / 2.0
    lo, hi = np.percentile(means, [100 * a, 100 * (1 - a)])
    return BootCI(point=point, lo=float(lo), hi=float(hi),
                  level=level, n_boot=n_boot, seed=seed)
