"""Lattice -- a UQ-calibration digital reference object (DRO) for IVIM.

Lattice provides synthetic ground-truth IVIM cohorts, five clean-room signal
generators (bi-exponential plus four alternative models that reduce to it at a
continuity limit), and a standardised calibration-evaluation interface so any
uncertainty-quantification method can be scored on a common reference.

It is a *resource*, not a research result: the cohorts and generators are solid
now and depend on no publication. See ``docs/DRO_SPEC.md`` and
``docs/POSITIONING.md``.
"""

from __future__ import annotations

from . import generators, cohort, evaluate, publication, osipi
from .cohort import (
    Cohort,
    make_cohort,
    sample_params,
    continuity_residual,
    DEFAULT_BVALUES,
    DEFAULT_SEED,
    PARAM_NAMES,
    PARAM_RANGES,
)
from .generators import FAMILIES
from .evaluate import (
    IVIMQuantileEstimator,
    DEFAULT_QUANTILE_LEVELS,
    to_scorer_inputs,
)

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "generators",
    "cohort",
    "evaluate",
    "publication",
    "osipi",
    "Cohort",
    "make_cohort",
    "sample_params",
    "continuity_residual",
    "DEFAULT_BVALUES",
    "DEFAULT_SEED",
    "PARAM_NAMES",
    "PARAM_RANGES",
    "FAMILIES",
    "IVIMQuantileEstimator",
    "DEFAULT_QUANTILE_LEVELS",
    "to_scorer_inputs",
]
