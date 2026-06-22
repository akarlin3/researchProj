"""Trust gate + action gate interfaces (consumed component: **Minos**).

Minos-core is literally *"the decision value of a calibrated per-voxel error bar (VoC)
plus a trust-gate (VoTG), on a synthetic treat/spare/escalate model."* Matrix consumes
exactly that pair, behind two interfaces:

  * **TrustGate** (Minos VoTG): ``trustworthy(state, cfg) -> bool[V]`` — is the
    calibrated uncertainty trustworthy at this voxel?
  * **ActionGate** (Minos decision): ``act(state, cfg) -> int[V]`` in {SPARE,TREAT,
    ESCALATE} — does the calibrated error bar change treat/spare, and where untrusted,
    is action *suppressed* (forced to ESCALATE)?

The real Minos gates (``minos.gate.{gate_fires,votg,gated_actions}`` + the
``minos.decision.bayes_action`` treat/spare/escalate map, applied half provisional @
PR #49) drop in behind these signatures **without touching loop.py**. See
``ASSUMPTIONS.md`` §Minos.

Implementations here:
  * ``TrustAllGate`` / ``PassThroughActionGate`` — CP1 stubs (no-op).
  * ``PlaceholderTrustGate`` / ``PlaceholderActionGate`` — CP2 placeholders, NOT-Minos.
"""
from __future__ import annotations

import numpy as np

from ..config import SPARE, TREAT, ESCALATE


# ---------------------------------------------------------------- CP1 stubs --
class TrustAllGate:
    """CP1 stub — trusts every voxel. NOT Minos."""

    name = "trust-all"
    label = "NOT-Minos (CP1 stub)"
    provisional = True

    def trustworthy(self, state, cfg):
        return np.ones(state.n_voxels, dtype=bool)


class PassThroughActionGate:
    """CP1 stub — holds every voxel at SPARE (no action). NOT Minos."""

    name = "hold-action"
    label = "NOT-Minos (CP1 stub)"
    provisional = True

    def act(self, state, cfg):
        a = np.full(state.n_voxels, SPARE, dtype=int)
        state.action_ungated = a.copy()
        return a


# -------------------------------------------------------- CP2 placeholders --
class PlaceholderTrustGate:
    """CP2 placeholder trust gate (Minos VoTG) — clearly labelled, NOT Minos.

    Flags a voxel UNtrustworthy when either the calibrated f error bar is too wide
    (``sigma_f > cfg.sigma_f_max``) or the fit lands in the high-D* identifiability
    regime (``Dstar_hat > cfg.dstar_untrust``) — the regime where the f-D* coupling
    makes the perfusion estimate unreliable. The real Minos VoTG replaces this.
    """

    name = "placeholder-trustgate"
    label = "NOT-Minos (CP2 placeholder; provisional pending Minos applied half PR#49)"
    provisional = True

    def trustworthy(self, state, cfg):
        sig_f = np.asarray(state.calib_sigma["f"], float)
        dstar_hat = np.asarray(state.mu["Dstar"], float)
        untrust = (sig_f > cfg.sigma_f_max) | (dstar_hat > cfg.dstar_untrust)
        return ~untrust


class PlaceholderActionGate:
    """CP2 placeholder action gate (Minos treat/spare/escalate) — NOT Minos.

    Bayes-action-shaped decision on the QoI = perfusion fraction f, using the
    *calibrated* interval ``mu_f +/- z*sigma_f``:

      * TREAT    if the interval lies confidently above ``cfg.f_treat``;
      * SPARE    if it lies confidently below ``cfg.f_spare``;
      * ESCALATE otherwise (borderline).

    Then the **trust gate suppresses action**: any UNtrustworthy voxel is forced to
    ESCALATE regardless of the raw decision. ``state.action_ungated`` records the
    pre-suppression decision so the CP4 check can measure the gate's effect.
    """

    name = "placeholder-actiongate"
    label = "NOT-Minos (CP2 placeholder; provisional pending Minos applied half PR#49)"
    provisional = True

    def act(self, state, cfg):
        mu_f = np.asarray(state.mu["f"], float)
        sig_f = np.asarray(state.calib_sigma["f"], float)
        lo = mu_f - cfg.z * sig_f
        hi = mu_f + cfg.z * sig_f

        ungated = np.full(state.n_voxels, ESCALATE, dtype=int)
        ungated[lo > cfg.f_treat] = TREAT
        ungated[hi < cfg.f_spare] = SPARE
        state.action_ungated = ungated.copy()

        gated = ungated.copy()
        gated[~state.trustworthy] = ESCALATE        # trust gate suppresses action
        return gated
