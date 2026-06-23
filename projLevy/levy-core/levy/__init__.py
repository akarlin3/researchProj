"""Levy-Core: identifiability of the fractional order under the diffusion-MRI
signal-decay forward model.

The contribution is the Fisher-information / Cramer-Rao layer for the fractional
order (stretched-exponential alpha; joint CTRW/fractional Bloch-Torrey (alpha, beta)
in Phase 3), estimated *jointly* with D and S0 under finite b-value sampling and
Rician noise -- an object the upstream Ouroboros tooling does not contain (see
POSITIONING.md), and whose likelihood differs from the fBm/Hurst CRB of
Coeurjolly-Istas 2001 (trajectory/MSD increments) by being b-indexed signal
attenuation.
"""
from __future__ import annotations

from .seeding import GLOBAL_SEED, make_rng

__all__ = ["GLOBAL_SEED", "make_rng"]
