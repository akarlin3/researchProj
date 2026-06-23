"""The regret-targeted deployment monitor ``M`` — ported from Minos.

``M`` is the utility-stakes-weighted L1 divergence of the reported-coordinate
distribution ``z = (mu - f_treat) / s`` at session ``k`` against the session-0
reference, weighted by ``omega(z) = k_under * phi(z)`` — a kernel localised to ~1
reported-sd of the decision threshold and scaled by the under-treatment cost.

This is the **regret-targeting**: ``M`` reacts to reported-mass changes *near the
decision threshold* (where the action can flip and regret is incurred) and is, by
construction, near-blind to far-from-threshold shifts that cannot move the optimal
action. The two baselines (ACI marginal coverage; WATCH exchangeability of the
whole score stream) target different functionals — that gap is the wedge.
"""
from __future__ import annotations

import numpy as np

from .config import SentinelConfig
from .seeding import make_rng


def z_edges(cfg: SentinelConfig) -> np.ndarray:
    return np.linspace(-cfg.z_range, cfg.z_range, cfg.z_bins + 1)


def _z_density(mu: np.ndarray, cfg: SentinelConfig, edges: np.ndarray) -> np.ndarray:
    z = (np.asarray(mu, float) - cfg.f_treat) / cfg.s_f
    h, _ = np.histogram(z, bins=edges, density=True)
    return h


def stakes_weights(cfg: SentinelConfig, edges: np.ndarray) -> np.ndarray:
    """``omega(z) = k_under * phi(z)`` at bin centres (Minos DESIGN_C 3.3)."""
    centres = 0.5 * (edges[:-1] + edges[1:])
    return cfg.k_under * np.exp(-0.5 * centres ** 2) / np.sqrt(2.0 * np.pi)


def monitor_M(mu_dep: np.ndarray, p_ref: np.ndarray, cfg: SentinelConfig) -> float:
    """``M = sum_b omega(z_b) |p_dep(z_b) - p_ref(z_b)| dz`` — non-negative staleness."""
    edges = z_edges(cfg)
    p_dep = _z_density(mu_dep, cfg, edges)
    omega = stakes_weights(cfg, edges)
    dz = edges[1] - edges[0]
    return float(np.sum(omega * np.abs(p_dep - p_ref)) * dz)


def reference_density(mu_ref: np.ndarray, cfg: SentinelConfig) -> np.ndarray:
    """Frozen session-0 reference histogram of the reported coordinate."""
    return _z_density(mu_ref, cfg, z_edges(cfg))


def calibrate_m_star(f_true: np.ndarray, p_ref: np.ndarray, cfg: SentinelConfig) -> float:
    """``m* = (1 - mon_alpha)`` quantile of ``M`` under the **no-drift** null.

    Draws ``mon_null_seeds`` fresh no-drift sessions (same patient, fresh measurement
    noise, drift switched off) and returns the upper-alpha quantile — controlling the
    no-drift false-alarm rate at ``mon_alpha``. This is the Minos m* construction.
    """
    rng = make_rng(cfg.seed + 5000)
    scores = np.empty(cfg.mon_null_seeds)
    for j in range(cfg.mon_null_seeds):
        mu = f_true + rng.normal(0.0, cfg.s_f, size=f_true.shape)  # no drift
        scores[j] = monitor_M(mu, p_ref, cfg)
    return float(np.quantile(scores, 1.0 - cfg.mon_alpha))
