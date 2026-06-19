"""Echo -- ground-truth-free *scale calibration* of conformal intervals via repeatability.

Echo validates that a deployed conformal IVIM interval is the right SIZE to capture
measurement irreproducibility -- a precision claim, provably blind to accuracy/bias and
distinct from Gauge's published width-RANK-tracks-repeatability check. See README.md and
ASSUMPTIONS.md. Every result that consumes Fashion/Minos/Gauge is PROVISIONAL.
"""
from __future__ import annotations

__version__ = "0.1.0"

from . import statistic, harness, provenance  # noqa: F401  (invivo imported on demand)

__all__ = ["statistic", "harness", "provenance", "__version__"]
