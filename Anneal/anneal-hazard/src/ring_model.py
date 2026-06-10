"""Nonlocally-coupled ring of N Kuramoto-Sakaguchi oscillators (Wolfrum-Omel'chenko 2011).

    dtheta_k/dt = omega - (1/(2P+1)) * sum_{|j-k|<=P, ring} sin(theta_k - theta_j + alpha)

with alpha = pi/2 - beta, omega = 0 (rotating frame). Coupling radius P (index space);
r = P/N is the coupling fraction. This is THE canonical system in which chimera states
are finite-N chaotic transients whose mean lifetime grows ~exponentially with N.

EXACT O(N) reduction of the top-hat sum:
    sum_{window} sin(theta_k - theta_j + alpha) = (2P+1) * Im[ e^{i(theta_k+alpha)} conj(Wbar_k) ]
    Wbar_k = (1/(2P+1)) sum_{|j-k|<=P} e^{i theta_j}            (circular moving average)
Wbar is computed with a circular uniform filter -> O(N) per evaluation, not O(N*P).
Wbar_k is ALSO the local order parameter; rho_k = |Wbar_k| measures local coherence.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import uniform_filter1d


@dataclass
class RingResult:
    t: np.ndarray            # decimated sample times
    # spatial summaries of the local order parameter rho_k = |Wbar_k| at each sample:
    rho_mean: np.ndarray
    rho_std: np.ndarray      # PRIMARY chimera indicator (large in chimera, ~0 at sync)
    rho_min: np.ndarray
    rho_max: np.ndarray
    rho_p05: np.ndarray
    rho_p95: np.ndarray
    R_global: np.ndarray     # |(1/N) sum e^{i theta}|
    theta_final: np.ndarray
    dt: float
    decimate: int
    N: int
    P: int
    beta: float
    seed: int
    stopped_early: bool


def local_meanfield(theta: np.ndarray, P: int) -> np.ndarray:
    """Circular moving average of e^{i theta} over window radius P -> Wbar_k (complex)."""
    c = np.exp(1j * theta)
    size = 2 * P + 1
    re = uniform_filter1d(c.real, size=size, mode="wrap")
    im = uniform_filter1d(c.imag, size=size, mode="wrap")
    return re + 1j * im


def _deriv(theta: np.ndarray, P: int, alpha: float, omega: float) -> np.ndarray:
    Wbar = local_meanfield(theta, P)
    return omega - np.imag(np.exp(1j * (theta + alpha)) * np.conj(Wbar))


def make_ring_ic(N: int, ic_cfg: dict, rng: np.random.Generator) -> np.ndarray:
    """Initial condition that seeds a chimera: a coherent base (small noise) with a
    contiguous incoherent arc (uniform random) of width frac*N.
    """
    frac = ic_cfg.get("incoherent_frac", 0.5)
    scale = ic_cfg.get("coherent_scale", 0.05)
    theta = rng.normal(0.0, scale, size=N)
    w = int(round(frac * N))
    start = ic_cfg.get("arc_start", None)
    if start is None:
        start = (N - w) // 2  # centered incoherent arc (deterministic given N, frac)
    idx = (np.arange(start, start + w)) % N
    theta[idx] = rng.uniform(-np.pi, np.pi, size=w)
    return theta


def _summarize(theta: np.ndarray, P: int) -> tuple:
    Wbar = local_meanfield(theta, P)
    rho = np.abs(Wbar)
    Rg = np.abs(np.exp(1j * theta).mean())
    p05, p95 = np.percentile(rho, [5, 95])
    return (float(rho.mean()), float(rho.std()), float(rho.min()), float(rho.max()),
            float(p05), float(p95), float(Rg))


def integrate_ring(
    theta0: np.ndarray,
    N: int,
    P: int,
    beta: float,
    dt: float,
    T_max: float,
    decimate: int,
    seed: int,
    omega: float = 0.0,
    stop_eps: float | None = None,
    dt_hold: float = 0.0,
) -> RingResult:
    """Fixed-step RK4. Records spatial summaries of rho_k every `decimate` steps.
    Optional early stop once rho_std < stop_eps continuously for dt_hold (death realized).
    """
    alpha = np.pi / 2.0 - beta
    n_steps = int(round(T_max / dt))
    x = theta0.astype(np.float64).copy()

    cols = {k: [] for k in
            ["t", "mean", "std", "min", "max", "p05", "p95", "R"]}

    below_start = None
    stopped_early = False

    for step in range(n_steps):
        if step % decimate == 0:
            m, s, mn, mx, p05, p95, Rg = _summarize(x, P)
            cols["t"].append(step * dt); cols["mean"].append(m); cols["std"].append(s)
            cols["min"].append(mn); cols["max"].append(mx)
            cols["p05"].append(p05); cols["p95"].append(p95); cols["R"].append(Rg)
            cur_std = s
        else:
            cur_std = None

        if stop_eps is not None:
            if cur_std is None:
                cur_std = float(local_meanfield(x, P).__abs__().std())
            t_now = step * dt
            if cur_std < stop_eps:
                if below_start is None:
                    below_start = t_now
                elif t_now - below_start >= dt_hold:
                    if step % decimate != 0:
                        m, s, mn, mx, p05, p95, Rg = _summarize(x, P)
                        cols["t"].append(t_now); cols["mean"].append(m); cols["std"].append(s)
                        cols["min"].append(mn); cols["max"].append(mx)
                        cols["p05"].append(p05); cols["p95"].append(p95); cols["R"].append(Rg)
                    stopped_early = True
                    break
            else:
                below_start = None

        k1 = _deriv(x, P, alpha, omega)
        k2 = _deriv(x + 0.5 * dt * k1, P, alpha, omega)
        k3 = _deriv(x + 0.5 * dt * k2, P, alpha, omega)
        k4 = _deriv(x + dt * k3, P, alpha, omega)
        x = x + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

    if not stopped_early:
        t_final = n_steps * dt
        if not cols["t"] or cols["t"][-1] != t_final:
            m, s, mn, mx, p05, p95, Rg = _summarize(x, P)
            cols["t"].append(t_final); cols["mean"].append(m); cols["std"].append(s)
            cols["min"].append(mn); cols["max"].append(mx)
            cols["p05"].append(p05); cols["p95"].append(p95); cols["R"].append(Rg)

    a = np.asarray
    return RingResult(
        t=a(cols["t"]), rho_mean=a(cols["mean"]), rho_std=a(cols["std"]),
        rho_min=a(cols["min"]), rho_max=a(cols["max"]),
        rho_p05=a(cols["p05"]), rho_p95=a(cols["p95"]), R_global=a(cols["R"]),
        theta_final=x, dt=dt, decimate=decimate, N=N, P=P, beta=beta, seed=seed,
        stopped_early=stopped_early,
    )
