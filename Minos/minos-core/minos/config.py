"""Frozen configuration for the Minos-Core synthetic decision model.

All numbers that define the toy live here. The dataclass is frozen so a run cannot
silently mutate its own parameters mid-experiment.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Tuple

from .seeding import GLOBAL_SEED


@dataclass(frozen=True)
class MinosConfig:
    # --- utility (Section 1) ---
    t1: float = 0.0          # spare | treat threshold
    t2: float = 2.0          # treat | escalate threshold
    k_under: float = 2.0     # under-treatment slope (penalised more)
    k_over: float = 1.0      # over-treatment slope

    # --- latent prior, 3-component Gaussian mixture (Section 2) ---
    # Symmetric about the decision midpoint (t1+t2)/2 = 1.0 so the calibrated
    # reported posterior (tau=1) is the decision-optimal error bar (DESIGN 6.2).
    mix_weights: Tuple[float, ...] = (0.35, 0.30, 0.35)
    mix_means: Tuple[float, ...] = (-1.0, 1.0, 3.0)
    mix_stds: Tuple[float, ...] = (1.0, 1.0, 1.0)

    # --- measurement (Section 2) ---
    s: float = 0.5           # intrinsic measurement spread
    tau: float = 1.0         # calibration knob (1 = calibrated)

    # --- shift (Section 3) ---
    delta: float = 0.0       # distribution-shift magnitude
    alpha: float = 0.5       # noise-inflation gain (overconfidence under shift)
    beta: float = 5.0        # downward-bias gain (systematic under-treatment)

    # --- trust-gate (Section 5) ---
    w_train_mean: float = 0.0
    w_train_std: float = 1.0
    q_gate: float = 0.995    # in-distribution quantile -> threshold (0.5% nominal FPR)

    # --- estimation ---
    n_voxels: int = 200_000
    seed: int = GLOBAL_SEED

    def __post_init__(self) -> None:
        assert self.t1 < self.t2, "need t1 < t2"
        assert self.k_under > self.k_over > 0, "asymmetry: under-treatment must cost more"
        assert self.s > 0 and self.tau > 0
        assert len(self.mix_weights) == len(self.mix_means) == len(self.mix_stds)
        assert abs(sum(self.mix_weights) - 1.0) < 1e-9, "mixture weights must sum to 1"
        assert 0.0 < self.q_gate < 1.0
        assert self.n_voxels > 0

    def replace(self, **kw) -> "MinosConfig":
        """Return a copy with overrides (frozen dataclasses need this helper)."""
        from dataclasses import replace as _replace

        return _replace(self, **kw)


DEFAULT = MinosConfig()
