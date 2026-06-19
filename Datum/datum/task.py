"""The fixed Datum benchmark task -- methods x data x metrics, frozen and versioned.

This is the definition that makes Datum a *benchmark* rather than a script: a
single frozen, versioned specification of the data substrate, the prediction
contract, the metrics (Fashion's ruler-as-standard), the conditioning axis, and
the curated baseline panel. A submission is scored against exactly this spec, so
reference numbers are comparable across methods and over time.

Bumping the substrate, levels, or metric set is a new task VERSION (and triggers
re-validation of every reference number).
"""
from __future__ import annotations

from dataclasses import dataclass, field

from datum.baselines import panel_keys
from datum.manifest import RULER

# Quantile levels a submission must predict; they span Fashion's nominal levels so
# central intervals can be read off at any of them. Ascending in (0, 1).
QUANTILE_LEVELS = (0.005, 0.025, 0.05, 0.10, 0.16, 0.25, 0.50,
                   0.75, 0.84, 0.90, 0.95, 0.975, 0.995)


@dataclass(frozen=True)
class TaskSpec:
    name: str
    version: str
    substrate: str                  # key into datum.substrate.SUBSTRATES
    n_train: int
    n_cal: int
    n_test: int
    seed: int
    quantile_levels: tuple
    nominal_levels: tuple           # Fashion ruler LEVELS (the standard)
    alpha: float                    # central-interval miss rate for headline coverage
    param_names: tuple
    metrics: tuple
    conditioning: str               # the identifiability-wall axis
    baselines: tuple
    n_bootstrap: int                # bootstrap resamples for reference-number CIs
    ruler: str = field(default_factory=lambda: f"{RULER['name']} v{RULER['version']}")

    def __post_init__(self):
        assert tuple(sorted(self.quantile_levels)) == tuple(self.quantile_levels), \
            "quantile_levels must be ascending"
        assert 0.0 < self.alpha < 1.0
        assert len(self.param_names) == 3, "IVIM params are (D, D*, f)"


TASK_V1 = TaskSpec(
    name="datum-ivim-calibration",
    version="v1",
    substrate="gauge_cohort",
    n_train=3000,
    n_cal=2000,
    n_test=3000,
    seed=20260613,                  # matches the Gauge cohort seed
    quantile_levels=QUANTILE_LEVELS,
    nominal_levels=tuple(RULER["nominal_levels"]),
    alpha=0.10,                     # headline central interval = 0.90 (a Fashion level)
    param_names=("D", "Dstar", "f"),
    metrics=("coverage", "coverage_gap", "ece", "sharpness",
             "mean_pinball", "mean_interval_score"),
    conditioning="Dstar_tercile",   # the high-D* identifiability wall
    baselines=panel_keys(),
    n_bootstrap=1000,
)
