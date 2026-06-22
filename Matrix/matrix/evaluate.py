"""Evaluation metrics + bootstrap CIs for the closed-loop run (load-bearing numbers).

Per guardrail 3, every load-bearing summary number is reported with a seeded bootstrap
confidence interval. These metrics operate only on the synthetic twin's known ground
truth and the per-iteration ``LoopState`` records — they make no clinical claim.
"""
from __future__ import annotations

import numpy as np

from .config import TREAT, SPARE, ESCALATE


def bootstrap_ci(values, stat=np.mean, n_boot=2000, seed=0, alpha=0.05):
    """Seeded bootstrap CI for a 1-D sample. Returns ``(point, lo, hi)``."""
    values = np.asarray(values, float)
    point = float(stat(values))
    if values.size == 0:
        return point, float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    n = values.size
    boots = np.array([stat(values[rng.integers(0, n, n)]) for _ in range(n_boot)])
    lo, hi = np.quantile(boots, [alpha / 2, 1 - alpha / 2])
    return point, float(lo), float(hi)


def auroc(scores, labels):
    """AUROC of continuous ``scores`` against boolean ``labels`` (rank statistic)."""
    scores = np.asarray(scores, float)
    labels = np.asarray(labels, bool)
    npos, nneg = int(labels.sum()), int((~labels).sum())
    if npos == 0 or nneg == 0:
        return float("nan")
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(scores.size, float)
    ranks[order] = np.arange(1, scores.size + 1)
    return (ranks[labels].sum() - npos * (npos + 1) / 2.0) / (npos * nneg)


def trust_gate_metrics(twin, state, cfg, seed=0):
    """Does the trust gate flag the (known) untrustworthy low-SNR zone?"""
    untrust = ~state.trustworthy
    lowsnr = twin.lowsnr
    # continuous score = calibrated sigma_f (what the gate thresholds)
    auc = auroc(state.calib_sigma["f"], lowsnr)
    fire_low = bootstrap_ci(untrust[lowsnr], seed=seed)
    fire_good = bootstrap_ci(untrust[~lowsnr], seed=seed + 1)
    return dict(auroc_untrust_vs_lowsnr=auc,
                fire_rate_lowsnr=fire_low, fire_rate_goodsnr=fire_good)


def suppression_metrics(state, seed=0):
    """The trust gate's effect: action-change rate on untrustworthy voxels, gated vs not.

    'Action-change' = the gate's decision is TREAT or SPARE (i.e. it acts), as opposed
    to ESCALATE (hold for review). The trust gate should drive this to ~0 on
    untrustworthy voxels while leaving trustworthy voxels free to act.
    """
    untrust = ~state.trustworthy
    acts_gated = np.isin(state.action, (TREAT, SPARE))
    acts_ungated = np.isin(state.action_ungated, (TREAT, SPARE))
    return dict(
        # without the gate, this fraction of untrustworthy voxels would be acted on:
        act_rate_untrust_ungated=bootstrap_ci(acts_ungated[untrust], seed=seed),
        # with the gate, it should be ~0:
        act_rate_untrust_gated=bootstrap_ci(acts_gated[untrust], seed=seed + 1),
        # trustworthy voxels remain free to act:
        act_rate_trust_gated=bootstrap_ci(acts_gated[~untrust], seed=seed + 2),
        n_suppressed=int(np.sum(acts_ungated & ~acts_gated)),
    )


def dose_warrant_metrics(state):
    """Dose changes only where warranted (TREAT boosts, SPARE de-escalates, ESCALATE holds)."""
    d = state.delta_dose
    return dict(
        treat_all_boost=bool(np.all(d[state.action == TREAT] > 0)) if np.any(state.action == TREAT) else True,
        spare_all_nonpos=bool(np.all(d[state.action == SPARE] <= 0)) if np.any(state.action == SPARE) else True,
        escalate_all_hold=bool(np.all(d[state.action == ESCALATE] == 0)) if np.any(state.action == ESCALATE) else True,
    )


def convergence_series(states):
    """Per-iteration closed-loop trajectory."""
    snaps = [s.snapshot() for s in states]
    return dict(
        mean_f_truth=[s["mean_f_truth"] for s in snaps],
        n_treat=[s["n_treat"] for s in snaps],
        n_escalate=[s["n_escalate"] for s in snaps],
        mean_dose=[s["mean_dose"] for s in snaps],
    )
