"""Trust-gate: an OOD signal on the acquisition feature, a conservative override
policy, the Value of the Trust-Gate (VoTG), and a shift-detection AUC.

Math: DESIGN.md Sections 5 and 6.4.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import norm

from .config import MinosConfig
from .decision import bayes_action
from .generative import BaseDraws, realise
from .utility import Action
from .voi import realised_utility


def gate_signal(w: np.ndarray, cfg: MinosConfig) -> np.ndarray:
    """Standardized one-sided OOD score ``g(w) = (w - m_w)/s_w`` (training stats).

    For unit-variance Gaussian features this is monotone in the train→deployment log
    density-ratio, so thresholding ``g`` is the density-ratio test.
    """
    return (np.asarray(w, dtype=float) - cfg.w_train_mean) / cfg.w_train_std


def gate_threshold(cfg: MinosConfig) -> float:
    """Threshold ``g*`` = the ``q_gate`` quantile of the in-distribution signal."""
    return float(norm.ppf(cfg.q_gate))


def gate_fires(w: np.ndarray, cfg: MinosConfig) -> np.ndarray:
    """Boolean per-voxel: does the gate fire (signal exceeds threshold)?"""
    return gate_signal(w, cfg) > gate_threshold(cfg)


def gated_actions(base: BaseDraws, cfg: MinosConfig, *, tau: float = 1.0,
                  delta: float = 0.0, shift=False) -> np.ndarray:
    """Posterior actions, overridden to ESCALATE wherever the gate fires."""
    mu, w = realise(base, cfg, delta=delta, shift=shift)
    actions = np.asarray(bayes_action(mu, tau * cfg.s, cfg)).copy()
    actions[gate_fires(w, cfg)] = int(Action.ESCALATE)
    return actions


def expected_utility_gated(base: BaseDraws, cfg: MinosConfig, *, tau: float = 1.0,
                           delta: float = 0.0, shift=False) -> float:
    actions = gated_actions(base, cfg, tau=tau, delta=delta, shift=shift)
    return float(np.mean(realised_utility(actions, base.theta, cfg)))


def votg(base: BaseDraws, cfg: MinosConfig, *, delta: float, tau: float = 1.0) -> float:
    """Value of the Trust-Gate ``VoTG(delta) = EU(gated) - EU(posterior)`` under a
    homogeneous shift ``delta`` applied to every voxel."""
    from .voi import expected_utility

    eu_gated = expected_utility_gated(base, cfg, tau=tau, delta=delta, shift=True)
    eu_post = expected_utility("posterior", base, cfg, tau=tau, delta=delta, shift=True)
    return eu_gated - eu_post


def _auc(scores: np.ndarray, labels: np.ndarray) -> float:
    """AUC via the Mann-Whitney U statistic (ties handled by average ranks)."""
    labels = np.asarray(labels, dtype=bool)
    n_pos = int(labels.sum())
    n_neg = labels.size - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, scores.size + 1)
    # average ranks within tied groups
    s_sorted = scores[order]
    i = 0
    while i < s_sorted.size:
        j = i + 1
        while j < s_sorted.size and s_sorted[j] == s_sorted[i]:
            j += 1
        if j - i > 1:
            ranks[order[i:j]] = 0.5 * (i + 1 + j)
        i = j
    sum_pos = ranks[labels].sum()
    return float((sum_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


def detection_auc(base: BaseDraws, cfg: MinosConfig, delta: float) -> float:
    """AUC of the gate signal against the latent shift mask, on a mixed population
    (first half shifted at ``delta``, second half in-distribution)."""
    n = base.theta.shape[0]
    mask = np.arange(n) < (n // 2)
    _, w = realise(base, cfg, delta=delta, shift=mask)
    return _auc(gate_signal(w, cfg), mask)
