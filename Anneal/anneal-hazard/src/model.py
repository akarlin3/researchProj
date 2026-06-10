"""Two-population Sakaguchi-Kuramoto, exact mean-field reduction, fixed-step RK4.

Finite-N equation (all-to-all), reduced via the identity
    (1/N) sum_j sin(theta_j - theta_i - alpha) = r * sin(phi - theta_i - alpha)
where r e^{i phi} = (1/N) sum_j e^{i theta_j} is the population order parameter.
This is an EXACT rewrite of the all-to-all sum (no approximation): we still evolve
all N oscillators per population, just in O(N) per step instead of O(N^2).

    dtheta_i^sigma/dt = mu r_sigma sin(phi_sigma - theta_i^sigma - alpha)
                      + nu r_sigma' sin(phi_sigma' - theta_i^sigma - alpha)

with sigma' != sigma, mu = K11 = K22, nu = K12 = K21, alpha = pi/2 - beta, omega = 0.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config_io import Params


def order_parameter(theta: np.ndarray) -> tuple[float, float]:
    """Return (r, phi) for one population's phase array."""
    z = np.exp(1j * theta).mean()
    return float(np.abs(z)), float(np.angle(z))


def _deriv(x: np.ndarray, N: int, p: Params) -> tuple[np.ndarray, float, float]:
    """RHS of the reduced ODE. x = concat(theta1, theta2). Also returns (r1, r2)
    computed at this state (free by-product, used for storage at step start)."""
    th1 = x[:N]
    th2 = x[N:]
    z1 = np.exp(1j * th1).mean()
    z2 = np.exp(1j * th2).mean()
    r1 = np.abs(z1); phi1 = np.angle(z1)
    r2 = np.abs(z2); phi2 = np.angle(z2)
    d = np.empty_like(x)
    # population 1: intra uses (r1,phi1), inter uses (r2,phi2)
    d[:N] = (p.mu * r1 * np.sin(phi1 - th1 - p.alpha)
             + p.nu * r2 * np.sin(phi2 - th1 - p.alpha))
    # population 2: intra uses (r2,phi2), inter uses (r1,phi1)
    d[N:] = (p.mu * r2 * np.sin(phi2 - th2 - p.alpha)
             + p.nu * r1 * np.sin(phi1 - th2 - p.alpha))
    return d, float(r1), float(r2)


@dataclass
class RunResult:
    t: np.ndarray        # decimated sample times
    r1: np.ndarray       # |order param| population 1 at sample times
    r2: np.ndarray       # |order param| population 2 at sample times
    theta_final: np.ndarray  # final full state (concat), for replay/inspection
    dt: float
    decimate: int
    N: int
    A: float
    seed: int
    stopped_early: bool  # True if integration halted before T_max (early-stop hit)


def make_initial_conditions(N: int, ic_cfg: dict, rng: np.random.Generator) -> np.ndarray:
    """Population 1 = coherent (Normal, r~1); population 2 = incoherent (Uniform, r~0)."""
    coh = ic_cfg["coherent"]
    inc = ic_cfg["incoherent"]
    th1 = rng.normal(loc=coh["loc"], scale=coh["scale"], size=N)
    th2 = rng.uniform(low=inc["low"], high=inc["high"], size=N)
    return np.concatenate([th1, th2])


def integrate(
    theta0: np.ndarray,
    p: Params,
    N: int,
    dt: float,
    T_max: float,
    decimate: int,
    seed: int,
    stop_eps: float | None = None,
    dt_hold: float = 0.0,
) -> RunResult:
    """Fixed-step RK4. Records (t, r1, r2) every `decimate` steps.

    Optional early stop: if stop_eps is given, integration halts once |r1-r2| has
    stayed < stop_eps for a continuous span of dt_hold (the death is fully realized).
    Used in the full ensemble to avoid wasting steps after the strictest-eps death;
    stop_eps should be the SMALLEST eps to be re-swept so all larger eps remain valid.
    """
    n_steps = int(round(T_max / dt))
    x = theta0.astype(np.float64).copy()

    t_store = []
    r1_store = []
    r2_store = []

    below_run_start = None  # time when current |r1-r2|<stop_eps run began
    stopped_early = False

    for step in range(n_steps):
        k1, r1, r2 = _deriv(x, N, p)
        if step % decimate == 0:
            t_store.append(step * dt)
            r1_store.append(r1)
            r2_store.append(r2)

        if stop_eps is not None:
            t_now = step * dt
            if abs(r1 - r2) < stop_eps:
                if below_run_start is None:
                    below_run_start = t_now
                elif t_now - below_run_start >= dt_hold:
                    # death fully realized; record this point and stop
                    if step % decimate != 0:
                        t_store.append(t_now)
                        r1_store.append(r1)
                        r2_store.append(r2)
                    stopped_early = True
                    break
            else:
                below_run_start = None

        k2, _, _ = _deriv(x + 0.5 * dt * k1, N, p)
        k3, _, _ = _deriv(x + 0.5 * dt * k2, N, p)
        k4, _, _ = _deriv(x + dt * k3, N, p)
        x = x + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    if not stopped_early:
        # record final state at T_max
        _, r1, r2 = _deriv(x, N, p)
        t_final = n_steps * dt
        if not t_store or t_store[-1] != t_final:
            t_store.append(t_final)
            r1_store.append(r1)
            r2_store.append(r2)

    return RunResult(
        t=np.asarray(t_store),
        r1=np.asarray(r1_store),
        r2=np.asarray(r2_store),
        theta_final=x,
        dt=dt,
        decimate=decimate,
        N=N,
        A=p.A,
        seed=seed,
        stopped_early=stopped_early,
    )
