"""Gnomon -- a clean-room ruler for reproduce-or-refute of Fashion's IVIM calibration.

Gnomon independently reimplements Fashion's calibration ruler (IVIM forward model,
boundary-railing NLLS diagnostic, Laplace + MCMC Bayesian posteriors, MAF amortized
posterior, and a coverage/ECE/sharpness scorer) **from spec, not from Fashion's or
Caliper's source**, in order to answer one question: does a from-scratch rebuild
reproduce Fashion's load-bearing numbers? The pinned targets and frozen tolerances
live in :mod:`gnomon.manifest`; the CP3 driver is :mod:`gnomon.reproduce`.

The heavy numerical modules (``forward``, ``cohort``, ``nlls``, ``bayes``, ``flow``,
``metrics``, ``bootstrap``, ``osipi``, ``reproduce``) are imported lazily so that
``import gnomon`` and the manifest work in a numpy-free environment.
"""
from __future__ import annotations

from . import manifest, _paths

__all__ = ["manifest", "_paths", "__version__"]
__version__ = "0.0.1"  # CP1 scaffold (PROVISIONAL; pre-reproduction)
