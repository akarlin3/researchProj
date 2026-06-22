"""The four-stage synthetic-twin closed loop, wired around a shared ``LoopState``.

    scan  ->  posterior(+ruler)  ->  trust gate + action gate  ->  dose replan  -> (re-scan)

Each stage is a small pure function ``stage_*(...)`` that reads the prior stage's fields
off the ``LoopState`` and writes its own. ``Interfaces`` bundles the three consumed
components (ruler=Fashion, trust+action=Minos, dose=Forge); swapping a real component for
its placeholder is a one-line change *here in the bundle*, never in a stage.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import MatrixConfig
from .fit import fit_scan
from .state import LoopState
from .twin import Twin
from .interfaces.ruler import PassThroughRuler, PlaceholderRuler
from .interfaces.gates import (TrustAllGate, PassThroughActionGate,
                               PlaceholderTrustGate, PlaceholderActionGate)
from .interfaces.dose import NoOpDoseEngine, PlaceholderDoseEngine


@dataclass
class Interfaces:
    """The three consumed-component implementations used for a run."""

    ruler: object
    trust_gate: object
    action_gate: object
    dose_engine: object

    @classmethod
    def passthrough(cls) -> "Interfaces":
        """CP1 wiring: every consumed component is a no-op pass-through stub."""
        return cls(PassThroughRuler(), TrustAllGate(),
                   PassThroughActionGate(), NoOpDoseEngine())

    @classmethod
    def placeholders(cls) -> "Interfaces":
        """CP2/CP3 wiring: clearly-labelled placeholders (NOT-Fashion/Minos/Forge)."""
        return cls(PlaceholderRuler(), PlaceholderTrustGate(),
                   PlaceholderActionGate(), PlaceholderDoseEngine())

    def provenance(self) -> dict:
        return {role: getattr(obj, "label", getattr(obj, "name", repr(obj)))
                for role, obj in (("ruler", self.ruler), ("trust_gate", self.trust_gate),
                                  ("action_gate", self.action_gate),
                                  ("dose_engine", self.dose_engine))}


# ---------------------------------------------------------------- stages -----
def stage_scan(twin: Twin, cfg: MatrixConfig, rng, state: LoopState) -> None:
    state.truth = twin.truth_snapshot()
    state.scan = twin.scan(rng)


def stage_posterior(cfg: MatrixConfig, state: LoopState, ruler) -> None:
    mu, raw_sigma = fit_scan(state.scan, cfg)
    state.mu, state.raw_sigma = mu, raw_sigma
    cal = ruler.calibrate(mu, raw_sigma, truth=state.truth)
    state.calib_sigma = cal["sigma"]
    state.interval = cal.get("interval")
    state.coverage = dict(coverage=cal.get("coverage", {}), ece=cal.get("ece", {}))


def stage_gates(cfg: MatrixConfig, state: LoopState, trust_gate, action_gate) -> None:
    state.trustworthy = trust_gate.trustworthy(state, cfg)
    state.action = action_gate.act(state, cfg)


def stage_dose(cfg: MatrixConfig, state: LoopState, dose_engine, twin: Twin) -> None:
    out = dose_engine.replan(twin.dose, state.action, state, cfg)
    state.replan = out["dose"]
    state.delta_dose = out["delta"]
    twin.apply_plan(state.replan)        # deliver -> twin evolves -> next scan sees it


# ----------------------------------------------------------------- driver ----
def run_iteration(twin: Twin, cfg: MatrixConfig, ifaces: Interfaces, rng,
                  iteration: int) -> LoopState:
    """One full pass through the four stages; returns the populated state."""
    state = LoopState(iteration=iteration, n_voxels=cfg.n_voxels)
    state.components = ifaces.provenance()
    stage_scan(twin, cfg, rng, state)
    stage_posterior(cfg, state, ifaces.ruler)
    stage_gates(cfg, state, ifaces.trust_gate, ifaces.action_gate)
    stage_dose(cfg, state, ifaces.dose_engine, twin)
    return state


def run_loop(cfg: MatrixConfig, ifaces: Interfaces | None = None,
             n_iter: int | None = None):
    """Build the twin and run the closed loop. Returns ``(twin, [LoopState, ...])``.

    The scan RNG is seeded off ``cfg.seed`` so a whole run is reproducible; the twin's
    ground truth is reproducible from ``cfg.seed`` independently (see ``Twin.build``).
    """
    if ifaces is None:
        ifaces = Interfaces.passthrough()
    n_iter = cfg.n_iter if n_iter is None else n_iter
    twin = Twin.build(cfg)
    rng = np.random.default_rng(cfg.seed + 12345)   # scan noise stream
    states = []
    for it in range(n_iter):
        states.append(run_iteration(twin, cfg, ifaces, rng, it))
    return twin, states
