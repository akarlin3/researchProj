"""LoopState — the single shared object threaded through the four loop stages.

Each stage *reads* the fields written by the stage before it and *writes* its own.
The CP1 "loop closes / end-to-end" check is, mechanically, that after a full pass every
field below is populated and no stage raised. ``snapshot`` records a compact per-iteration
summary for the closed-loop evaluation (CP4).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np


@dataclass
class LoopState:
    iteration: int
    n_voxels: int

    # stage 1 — scan
    truth: Optional[dict] = None          # read-only ground-truth snapshot (eval only)
    scan: Optional[np.ndarray] = None     # (n_b, V, n_noise)

    # stage 2 — posterior + ruler (Fashion)
    mu: Optional[dict] = None             # raw posterior centre per param
    raw_sigma: Optional[dict] = None      # raw posterior spread per param
    calib_sigma: Optional[dict] = None    # ruler-calibrated spread per param
    interval: Optional[dict] = None       # {param: (lo, hi)} calibrated interval
    coverage: Optional[dict] = None       # ruler coverage/ECE readout

    # stage 3 — trust gate + action gate (Minos)
    trustworthy: Optional[np.ndarray] = None   # (V,) bool
    action: Optional[np.ndarray] = None        # (V,) int in {SPARE,TREAT,ESCALATE}
    action_ungated: Optional[np.ndarray] = None  # action the gate *would* give w/o trust gate

    # stage 4 — dose replan (Forge-shaped)
    replan: Optional[np.ndarray] = None   # (V,) new prescription (Gy)
    delta_dose: Optional[np.ndarray] = None  # (V,) replan - previous dose

    # provenance: which interface implementations ran this iteration
    components: dict = field(default_factory=dict)

    def is_complete(self) -> bool:
        """True iff every stage wrote its output (the end-to-end check)."""
        return all(v is not None for v in (
            self.scan, self.mu, self.raw_sigma, self.calib_sigma,
            self.trustworthy, self.action, self.replan, self.delta_dose))

    def snapshot(self) -> dict:
        """Compact per-iteration record for the closed-loop summary."""
        from .config import TREAT, SPARE, ESCALATE
        return dict(
            iteration=self.iteration,
            n_treat=int(np.sum(self.action == TREAT)),
            n_spare=int(np.sum(self.action == SPARE)),
            n_escalate=int(np.sum(self.action == ESCALATE)),
            n_untrust=int(np.sum(~self.trustworthy)),
            mean_dose=float(np.mean(self.replan)),
            mean_f_truth=float(np.mean(self.truth["f"])) if self.truth else float("nan"),
            components=dict(self.components),
        )
