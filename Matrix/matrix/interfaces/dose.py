"""Dose-replan interface (consumed component: **Forge**, the dose engine).

Contract
--------
``replan(current_dose, action, state, cfg) -> dict`` consumes the gated per-voxel
decision and returns a new prescription::

    {"dose": new_dose[V], "delta": new_dose - current_dose}

**Forge's real dose engine does not exist yet** — Forge today is only a Monte-Carlo
dose-simulation *feasibility* benchmark (timing + Electron Return Effect); the dose
surrogate is *deferred to 2027, not built*. So this interface ships with a transparent,
analytic placeholder, clearly labelled NOT-Forge. When Forge's engine lands it drops in
behind this signature **without touching loop.py**. See ``ASSUMPTIONS.md`` §Forge.

Implementations here:
  * ``NoOpDoseEngine``      — CP1 stub: returns the current plan unchanged.
  * ``PlaceholderDoseEngine`` — CP3 placeholder: analytic prescription update, NOT-Forge.
"""
from __future__ import annotations

import numpy as np

from ..config import TREAT, SPARE


class NoOpDoseEngine:
    """CP1 stub — never changes the plan. NOT Forge."""

    name = "noop-dose"
    label = "NOT-Forge (CP1 stub)"
    provisional = True

    def replan(self, current_dose, action, state, cfg):
        dose = np.asarray(current_dose, float).copy()
        return dict(dose=dose, delta=np.zeros_like(dose))


class PlaceholderDoseEngine:
    """CP3 placeholder dose engine — transparent analytic surrogate, **NOT** Forge.

    A deliberately simple, isolated stand-in for Forge's Monte-Carlo dose engine:

      * TREAT    -> boost toward ``dose_baseline + dose_boost`` (capped at ``dose_max``);
      * SPARE    -> de-escalate by ``dose_spare_cut`` (floored at ``dose_min``);
      * ESCALATE -> hold the current dose (no change pending review).

    It is purely analytic — no geometry, no transport, no Electron Return Effect. That
    physics is exactly what Forge's real engine will provide; this placeholder makes the
    loop *close* without it and isolates the swap to one object.
    """

    name = "placeholder-dose"
    label = "NOT-Forge (CP3 analytic placeholder; Forge dose engine deferred to 2027)"
    provisional = True

    def replan(self, current_dose, action, state, cfg):
        cur = np.asarray(current_dose, float)
        action = np.asarray(action, int)
        new = cur.copy()
        treat = action == TREAT
        spare = action == SPARE
        new[treat] = np.minimum(cfg.dose_baseline + cfg.dose_boost, cfg.dose_max)
        new[spare] = np.maximum(cur[spare] - cfg.dose_spare_cut, cfg.dose_min)
        # ESCALATE voxels: untouched (hold).
        return dict(dose=new, delta=new - cur)
