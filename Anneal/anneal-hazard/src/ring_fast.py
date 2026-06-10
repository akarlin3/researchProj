"""Numba-accelerated RK4 integrator for the top-hat nonlocal ring (Wolfrum-Omel'chenko).

Same model as ring_model.py, but the circular windowed mean field is computed with an
O(N) sliding running sum inside njit code, and the whole RK4 loop is compiled. This is the
~30x speedup the lifetime study needs (chimera transients require large T_max). Validated
against the numpy reference in ring_model.py (see __main__).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit


@njit(cache=True)
def _deriv(theta, P, alpha):
    N = theta.shape[0]
    cr = np.cos(theta)
    ci = np.sin(theta)
    W = 2 * P + 1
    inv = 1.0 / W
    # initial centered window sum for k=0: indices -P..P (mod N)
    sr = 0.0
    si = 0.0
    for d in range(-P, P + 1):
        idx = d % N
        sr += cr[idx]
        si += ci[idx]
    out = np.empty(N)
    for k in range(N):
        mr = sr * inv
        mi = si * inv
        a = np.cos(theta[k] + alpha)
        b = np.sin(theta[k] + alpha)
        # dtheta_k = -Im[ e^{i(theta_k+alpha)} * conj(Wbar) ] = -(b*mr - a*mi)
        out[k] = -(b * mr - a * mi)
        rem = (k - P) % N
        add = (k + 1 + P) % N
        sr += cr[add] - cr[rem]
        si += ci[add] - ci[rem]
    return out


@njit(cache=True)
def _rho_summ(theta, P):
    """Return (mean, std, min, max, R_global) of the local coherence field rho_k=|Wbar_k|."""
    N = theta.shape[0]
    cr = np.cos(theta)
    ci = np.sin(theta)
    W = 2 * P + 1
    inv = 1.0 / W
    sr = 0.0
    si = 0.0
    for d in range(-P, P + 1):
        idx = d % N
        sr += cr[idx]
        si += ci[idx]
    s1 = 0.0
    s2 = 0.0
    rmin = 1e9
    rmax = -1e9
    gr = 0.0
    gi = 0.0
    for k in range(N):
        mr = sr * inv
        mi = si * inv
        rho = (mr * mr + mi * mi) ** 0.5
        s1 += rho
        s2 += rho * rho
        if rho < rmin:
            rmin = rho
        if rho > rmax:
            rmax = rho
        gr += cr[k]
        gi += ci[k]
        rem = (k - P) % N
        add = (k + 1 + P) % N
        sr += cr[add] - cr[rem]
        si += ci[add] - ci[rem]
    mean = s1 / N
    var = s2 / N - mean * mean
    if var < 0.0:
        var = 0.0
    std = var ** 0.5
    Rg = ((gr / N) ** 2 + (gi / N) ** 2) ** 0.5
    return mean, std, rmin, rmax, Rg


@njit(cache=True)
def _integrate(theta0, P, alpha, dt, n_steps, decimate, stop_eps, hold_samples):
    """RK4 loop with decimated storage of rho summaries and optional early stop.

    Early stop is evaluated on the DECIMATED samples (same data the detector sees), so it
    is consistent with detect_death_ring. We stop only after the below-eps run has lasted
    hold_samples + BUFFER decimated samples, guaranteeing the stored run spans >= dt_hold
    (the truncation-at-the-boundary bug that previously caused deaths to be missed)."""
    BUFFER = 6  # extra decimated samples stored past hold completion (clean re-sweep)
    x = theta0.copy()
    cap = n_steps // decimate + 2
    t = np.empty(cap)
    a_mean = np.empty(cap)
    a_std = np.empty(cap)
    a_min = np.empty(cap)
    a_max = np.empty(cap)
    a_R = np.empty(cap)
    m = 0
    below = 0  # consecutive DECIMATED samples with std<stop_eps
    stopped = 0
    use_stop = stop_eps > 0.0

    for step in range(n_steps):
        if step % decimate == 0:
            mean, std, rmn, rmx, Rg = _rho_summ(x, P)
            t[m] = step * dt
            a_mean[m] = mean
            a_std[m] = std
            a_min[m] = rmn
            a_max[m] = rmx
            a_R[m] = Rg
            m += 1
            if use_stop:
                if std < stop_eps:
                    below += 1
                    if below >= hold_samples + BUFFER:
                        stopped = 1
                        break
                else:
                    below = 0

        k1 = _deriv(x, P, alpha)
        k2 = _deriv(x + 0.5 * dt * k1, P, alpha)
        k3 = _deriv(x + 0.5 * dt * k2, P, alpha)
        k4 = _deriv(x + dt * k3, P, alpha)
        x = x + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    if stopped == 0:
        mean, std, rmn, rmx, Rg = _rho_summ(x, P)
        tf = n_steps * dt
        if m == 0 or t[m - 1] != tf:
            t[m] = tf
            a_mean[m] = mean
            a_std[m] = std
            a_min[m] = rmn
            a_max[m] = rmx
            a_R[m] = Rg
            m += 1

    return t[:m], a_mean[:m], a_std[:m], a_min[:m], a_max[:m], a_R[:m], x, stopped


@dataclass
class RingResultFast:
    t: np.ndarray
    rho_mean: np.ndarray
    rho_std: np.ndarray
    rho_min: np.ndarray
    rho_max: np.ndarray
    R_global: np.ndarray
    theta_final: np.ndarray
    N: int
    P: int
    beta: float
    seed: int
    stopped_early: bool


def integrate_ring_fast(theta0, N, P, beta, dt, T_max, decimate, seed,
                        stop_eps=None, dt_hold=0.0):
    alpha = np.pi / 2.0 - beta
    n_steps = int(round(T_max / dt))
    # hold in decimated samples; +1 so a full run spans >= dt_hold (matches detector)
    hold_samples = int(round(dt_hold / (decimate * dt))) + 1
    se = stop_eps if stop_eps is not None else -1.0
    t, mean, std, rmn, rmx, R, xf, stopped = _integrate(
        np.ascontiguousarray(theta0, dtype=np.float64), int(P), float(alpha),
        float(dt), int(n_steps), int(decimate), float(se), int(hold_samples))
    return RingResultFast(t, mean, std, rmn, rmx, R, xf, N, P, beta, seed, bool(stopped))


if __name__ == "__main__":
    # validate against the numpy reference and time it
    import time
    from src.ring_model import make_ring_ic, integrate_ring

    N, P, beta = 128, 19, 0.08
    rng = np.random.default_rng(20260609)
    th0 = make_ring_ic(N, {"incoherent_frac": 0.5, "coherent_scale": 0.05}, rng)

    ref = integrate_ring(th0, N, P, beta, 0.05, 300.0, decimate=10, seed=0)
    fast = integrate_ring_fast(th0, N, P, beta, 0.05, 300.0, decimate=10, seed=0)  # compile
    d_std = np.max(np.abs(ref.rho_std - fast.rho_std))
    d_mean = np.max(np.abs(ref.rho_mean - fast.rho_mean))
    print(f"max|Δrho_std|={d_std:.2e}  max|Δrho_mean|={d_mean:.2e}  (numpy vs numba)")

    t0 = time.time(); integrate_ring(th0, N, P, beta, 0.05, 2000.0, decimate=10, seed=0); t_np = time.time() - t0
    t0 = time.time(); integrate_ring_fast(th0, N, P, beta, 0.05, 2000.0, decimate=10, seed=0); t_nb = time.time() - t0
    print(f"T=2000 N=128: numpy {t_np:.2f}s  numba {t_nb:.3f}s  speedup x{t_np/t_nb:.0f}")
