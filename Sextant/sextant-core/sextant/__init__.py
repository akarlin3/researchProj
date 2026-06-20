"""Sextant — boundary-railing as the primary, assumption-free IVIM diagnostic.

A *sextant* measures an angle to the horizon; here we measure how often a
conventional non-linear least-squares (NLLS) bi-exponential IVIM fit runs aground
on its own parameter *horizon* — the optimisation bounds. That a large fraction
of fits rail to a bound is a fact about the optimiser, not about any calibration
model, so it cannot be "overextended": it needs no ground truth and no
trust-in-our-noise-model argument.

This package re-centres Fashion's analysis: the boundary-railing computation is
reused **read-only** from Fashion (see :mod:`sextant.fashion_reuse`), promoted to
the primary claim, and replicated on open human-abdominal diffusion MRI. The
calibration "ruler" (coverage / ECE / sharpness) is demoted to a scoped
secondary section (:mod:`sextant.ruler`).
"""
from .seeding import GLOBAL_SEED, make_rng

__all__ = ["GLOBAL_SEED", "make_rng"]
__version__ = "0.1.0"
