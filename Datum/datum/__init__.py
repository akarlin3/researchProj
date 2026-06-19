"""Datum -- the IVIM uncertainty-quantification calibration benchmark.

Datum defines a *fixed task* (methods x data x metrics), a *curated baseline
panel*, *reference numbers*, and a *submission interface* for scoring a new IVIM
uncertainty method's calibration. It is built **on Fashion's calibration ruler**
(packaged read-only by Caliper) and runs on a synthetic data substrate (Gauge's
cohort now; Lattice when it is built).

Because Fashion's ruler is *in review*, every reference number Datum produces is
**PROVISIONAL** until the ruler locks. See ``datum.manifest`` (the pins),
``datum.provisional`` (the stamping), and ``ASSUMPTIONS.md`` (the SOLID vs
assumption-dependent split).

Module map
----------
* ``datum.task``       -- the frozen, versioned task spec (``TASK_V1``).
* ``datum.substrate``  -- read-only data adapters (Gauge cohort, OSIPI DRO, Lattice stub).
* ``datum.ruler``      -- read-only scoring adapter over ``caliper.metrics``.
* ``datum.baselines``  -- the curated baseline registry.
* ``datum.manifest``   -- machine-readable assumption pins + provisional policy.
* ``datum.provisional``-- PROVISIONAL stamping for ruler-derived numbers.
* ``datum._paths``     -- read-only sibling bootstrap for Caliper/Gauge.
"""
from __future__ import annotations

__version__ = "0.1.0"

from datum import manifest  # noqa: F401  (cheap, no heavy deps)

__all__ = ["__version__", "manifest"]
