"""Vernier -- calibration-aware acquisition design for IVIM diffusion MRI.

Vernier asks an acquisition-design question its CRLB-optimal predecessors do not:
at *matched scan-time* and *matched Cramer-Rao precision*, do different diffusion
b-value schemes yield differently-**calibrated** uncertainty *after* conformal
correction? Variance-optimal (CRLB) design and information-gain (BED/EIG) design
both optimise the point/precision; Vernier optimises the trustworthiness of the
error bar and its downstream decision value.

The package is numpy-only and reuses Caliper read-only (the calibration ruler,
conformal wrappers, synthetic IVIM forward model, and reference estimator). It
imports neither the in-review Fashion/Gauge/Minos code nor any clinical data --
the feasibility gate is publication-independent and PHI-free by construction.

Modules
-------
* ``vernier._paths``   -- read-only Caliper path wiring (the single dependency chokepoint).
* ``vernier.schemes``  -- b-value scheme registry, scan-time model, segmented-fit validation.
* ``vernier.crlb``     -- IVIM Fisher information matrix and Cramer-Rao lower bounds.
"""
from __future__ import annotations

__version__ = "0.0.1"

__all__ = ["__version__"]
