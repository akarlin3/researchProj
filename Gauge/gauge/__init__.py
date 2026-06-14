"""Gauge: distribution-free conformal coverage for IVIM parameter maps.

Gauge 01 (foundation) ships a self-contained bi-exponential IVIM forward model,
a labeled synthetic cohort, NLLS / quantile base estimators, and split-conformal
+ CQR prediction intervals with an empirical coverage-validation harness.
"""

__all__ = ["forward", "cohort", "estimators", "conformal", "monitor",
           "robustness", "invivo"]
