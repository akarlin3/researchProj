"""Caliper -- an IVIM uncertainty-quantification calibration toolkit.

The public surface is intentionally small:

* ``caliper.metrics``  -- the model-agnostic calibration ruler (numpy only).
* ``caliper.forward``  -- bi-exponential IVIM forward model + synthetic cohorts.
* ``caliper.conformal``-- split-conformal / CQR / Mondrian coverage correction.
* ``caliper.estimator_reference`` -- an over-confident segmented-fit IVIM
  estimator (numpy only); the torch-free device-under-test for the conformal
  layer.
* ``caliper.estimator_maf`` -- a masked-autoregressive-flow posterior (needs
  the optional ``[estimator]`` extra: torch).

Everything except the torch MAF is imported eagerly; the torch-dependent
estimator is imported lazily so the core stays numpy-only.
"""
from __future__ import annotations

from . import conformal, estimator_reference, forward, metrics

__all__ = [
    "metrics",
    "forward",
    "conformal",
    "estimator_reference",
    "__version__",
]
__version__ = "0.1.0"
