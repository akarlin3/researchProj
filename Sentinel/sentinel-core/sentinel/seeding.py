"""Deterministic seeding for projSentinel (mirrors Minos/Procrustes convention).

Every stochastic draw in projSentinel goes through ``make_rng`` so the whole
separation experiment is bit-reproducible from ``GLOBAL_SEED``. No bare
``np.random`` anywhere.
"""
from __future__ import annotations

import numpy as np

GLOBAL_SEED = 20260623


def make_rng(seed: int | None = None) -> np.random.Generator:
    """Return a fresh ``numpy`` Generator seeded from ``GLOBAL_SEED`` (or override)."""
    return np.random.default_rng(GLOBAL_SEED if seed is None else seed)
