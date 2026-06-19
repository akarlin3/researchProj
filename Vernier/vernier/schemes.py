"""vernier.schemes -- IVIM b-value acquisition schemes and their scan-time cost.

Numpy-only, no Caliper dependency (so it imports without wiring). A
:class:`BScheme` bundles a b-value schedule with per-b averages and exposes the
two quantities Vernier's feasibility gate controls for:

* **scan-time** -- proportional to the total number of acquired volumes
  ``sum(averages)``; with a fixed per-volume time this is the scanner-minute cost.
* **segmented-fit support** -- the Caliper reference estimator is a *segmented*
  least-squares fit, which needs a high-b tissue window (>= 2 b-values at/above
  ``high_b_min``) and a low-b perfusion window (>= 2 positive b-values at/below
  ``low_b_max``). Schemes that cannot be segment-fit are rejected up front, so
  the gate never compares a scheme the device-under-test cannot consume.

The registry below provides three **matched-scan-time** schemes (all 11 b-values,
so identical scanner-minutes) that differ only in how they spend that budget --
perfusion-weighted (dense low-b), balanced/clinical, tissue-weighted (dense
high-b). They are the design-rationale candidates for the feasibility sweep;
matched-*CRLB* triples are selected from a wider candidate pool at gate time by
:mod:`vernier.crlb`. A separate set spans scan-time (7/11/15/22 b-values) for the
efficiency-frontier (decision-value-per-scan-minute) analysis.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "DEFAULT_SECONDS_PER_VOLUME",
    "BScheme",
    "MATCHED_SCANTIME",
    "EFFICIENCY_FRONTIER",
    "CANDIDATE_POOL",
]

# Nominal per-volume acquisition time (a single b-value, single average) for a
# multi-slice abdominal DWI volume, in seconds. A constant placeholder used only
# to turn acquisition *counts* into scanner-*minutes*; it cancels in any
# matched-scan-time comparison and is stated explicitly as an assumption when it
# does not (the efficiency frontier).
DEFAULT_SECONDS_PER_VOLUME = 12.0


@dataclass(frozen=True)
class BScheme:
    """An IVIM b-value acquisition scheme.

    Parameters
    ----------
    name : human-readable label.
    bvalues : tuple of b-values in s/mm^2 (ascending, includes b=0).
    averages : optional per-b signal averages; defaults to one average per
        b-value. Length must match ``bvalues``.
    """

    name: str
    bvalues: tuple[float, ...]
    averages: tuple[int, ...] | None = None

    def __post_init__(self) -> None:
        b = np.asarray(self.bvalues, dtype=float)
        if b.ndim != 1 or b.size < 2:
            raise ValueError(f"{self.name}: need >=2 b-values")
        if np.any(np.diff(b) <= 0):
            raise ValueError(f"{self.name}: b-values must be strictly ascending")
        if np.any(b < 0):
            raise ValueError(f"{self.name}: b-values must be non-negative")
        if self.averages is not None:
            a = np.asarray(self.averages)
            if a.shape != b.shape:
                raise ValueError(f"{self.name}: averages must match bvalues length")
            if np.any(a < 1):
                raise ValueError(f"{self.name}: averages must be >= 1")

    # --- arrays ---------------------------------------------------------------
    @property
    def b(self) -> np.ndarray:
        """b-value schedule as a float array (s/mm^2)."""
        return np.asarray(self.bvalues, dtype=float)

    @property
    def averages_arr(self) -> np.ndarray:
        """Per-b averages as an int array (ones if unspecified)."""
        if self.averages is None:
            return np.ones(len(self.bvalues), dtype=int)
        return np.asarray(self.averages, dtype=int)

    # --- cost -----------------------------------------------------------------
    @property
    def n_b(self) -> int:
        """Number of distinct b-values."""
        return len(self.bvalues)

    @property
    def n_acquisitions(self) -> int:
        """Total acquired volumes = sum of per-b averages (the scan-time unit)."""
        return int(self.averages_arr.sum())

    def scan_minutes(self, seconds_per_volume: float = DEFAULT_SECONDS_PER_VOLUME) -> float:
        """Scanner-minutes = n_acquisitions * seconds_per_volume / 60."""
        return self.n_acquisitions * float(seconds_per_volume) / 60.0

    # --- validity -------------------------------------------------------------
    def supports_segmented_fit(self, high_b_min: float = 200.0, low_b_max: float = 50.0) -> bool:
        """True if the scheme has the windows the segmented reference fit needs."""
        b = self.b
        n_high = int(np.sum(b >= high_b_min))
        n_low = int(np.sum((b <= low_b_max) & (b > 0.0)))
        return n_high >= 2 and n_low >= 2

    def require_segmented_fit(self, high_b_min: float = 200.0, low_b_max: float = 50.0) -> "BScheme":
        """Return self if segment-fittable, else raise (fail fast at gate setup)."""
        if not self.supports_segmented_fit(high_b_min, low_b_max):
            raise ValueError(
                f"{self.name}: needs >=2 b>={high_b_min} and >=2 0<b<={low_b_max} "
                f"for a segmented fit; got b={list(self.bvalues)}"
            )
        return self


# --------------------------------------------------------------------------- #
# Matched-scan-time schemes (all 11 b-values, 1 average -> identical minutes).
# Same budget, spent differently: this is the acquisition-design axis the gate
# probes. Each supports the segmented fit.
# --------------------------------------------------------------------------- #
MATCHED_SCANTIME: tuple[BScheme, ...] = (
    BScheme("perfusion-weighted",
            (0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 200.0, 500.0, 800.0)),
    BScheme("balanced-clinical",
            (0.0, 10.0, 20.0, 30.0, 50.0, 80.0, 120.0, 200.0, 400.0, 600.0, 800.0)),
    BScheme("tissue-weighted",
            (0.0, 20.0, 50.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0)),
)

# Efficiency-frontier schemes spanning scan-time (7/11/15/22 b-values). Used for
# decision-value-per-scan-minute, where the per-volume time does NOT cancel.
EFFICIENCY_FRONTIER: tuple[BScheme, ...] = (
    BScheme("sparse-7",
            (0.0, 25.0, 50.0, 100.0, 200.0, 500.0, 800.0)),
    BScheme("clinical-11",
            (0.0, 10.0, 20.0, 30.0, 50.0, 80.0, 120.0, 200.0, 400.0, 600.0, 800.0)),
    BScheme("rich-15",
            (0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0, 150.0, 200.0,
             300.0, 400.0, 550.0, 700.0, 800.0)),
    BScheme("dense-22",
            (0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 100.0, 120.0,
             150.0, 200.0, 250.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0,
             900.0, 1000.0)),
)

# A wider candidate pool (all 11 b-values / matched scan-time) from which
# matched-CRLB triples are selected at gate time. Includes the design-rationale
# schemes plus interpolating variants so a CRLB-matched-but-geometry-different
# triple can be found without tuning toward a divergence outcome.
CANDIDATE_POOL: tuple[BScheme, ...] = MATCHED_SCANTIME + (
    BScheme("perfusion-lean",
            (0.0, 10.0, 25.0, 40.0, 55.0, 70.0, 90.0, 200.0, 350.0, 550.0, 800.0)),
    BScheme("tissue-lean",
            (0.0, 15.0, 40.0, 90.0, 200.0, 300.0, 420.0, 540.0, 660.0, 780.0, 900.0)),
    BScheme("wide-ends",
            (0.0, 15.0, 30.0, 45.0, 60.0, 90.0, 250.0, 400.0, 550.0, 700.0, 850.0)),
    BScheme("mid-heavy",
            (0.0, 20.0, 40.0, 70.0, 110.0, 160.0, 220.0, 300.0, 420.0, 600.0, 800.0)),
)


def by_name(schemes: tuple[BScheme, ...], name: str) -> BScheme:
    """Look up a scheme by name in a tuple (raises KeyError if absent)."""
    for s in schemes:
        if s.name == name:
            return s
    raise KeyError(name)
