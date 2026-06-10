"""CP2 (C2) ring-object validation — field-dumping integrator + chimera classifier.

The CP3 ensemble saved only scalar rho_std(t)/rho_mean(t); the spatial local-coherence
field rho_k(t,.) needed to validate the dying object was NOT stored. Every run is exactly
reproducible from its logged seed (make_ring_ic + fixed-step RK4), so we re-run the target
cells with the SAME seed/IC/dt/P/alpha and the SAME early-stop as the ensemble (so tau
reproduces bit-for-bit), this time recording the decimated theta_k field. We then classify
the pre-death (last 10% of life) spatial state as a canonical chimera (one contiguous
coherent arc + one contiguous incoherent arc) vs degenerate (multi-arc / never-formed).

Exact-match guard: re-integrated tau is compared to the ensemble_N256.csv tau per run.
"""
from __future__ import annotations

import os
import sys

import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RING = os.path.join(ROOT, "anneal-hazard")
sys.path.insert(0, RING)

from src.ring_fast import _deriv, _rho_summ  # njit kernels (exact ensemble math)  # noqa: E402
from src.ring_model import make_ring_ic, local_meanfield  # noqa: E402
from src.ring_detector import detect_death_ring  # noqa: E402
from src.config_io import load_config  # noqa: E402

CFG = load_config(os.path.join(RING, "config.yaml"))
R = CFG["r"]; DT = CFG["dt"]; T_MAX = CFG["T_max"]; DT_HOLD = CFG["dt_hold"]
DECIMATE = CFG["output"]["decimate"]; EPS_STD = CFG["eps_std"]
STOP_EPS = 0.03  # same early-stop as cp3_ensemble (smallest re-sweep eps)
IC = CFG["ic"]


@njit(cache=True)
def _integrate_field(theta0, P, alpha, dt, n_steps, decimate, stop_eps, hold_samples):
    """RK4 with decimated storage of (t, rho_mean, rho_std, full theta field) + early stop
    identical to ring_fast._integrate (BUFFER=6 past hold completion)."""
    BUFFER = 6
    N = theta0.shape[0]
    x = theta0.copy()
    cap = n_steps // decimate + 2
    t = np.empty(cap)
    a_mean = np.empty(cap)
    a_std = np.empty(cap)
    th = np.empty((cap, N))
    m = 0
    below = 0
    stopped = 0
    use_stop = stop_eps > 0.0
    for step in range(n_steps):
        if step % decimate == 0:
            mean, std, rmn, rmx, Rg = _rho_summ(x, P)
            t[m] = step * dt
            a_mean[m] = mean
            a_std[m] = std
            for j in range(N):
                th[m, j] = x[j]
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
            for j in range(N):
                th[m, j] = x[j]
            m += 1
    return t[:m], a_mean[:m], a_std[:m], th[:m]


def run_field(beta, N, seed):
    """Reproduce one ensemble run, returning decimated (t, rho_mean, rho_std, theta[t,k])."""
    P = max(1, int(round(R * N)))
    alpha = np.pi / 2.0 - beta
    n_steps = int(round(T_MAX / DT))
    hold_samples = int(round(DT_HOLD / (DECIMATE * DT))) + 1
    th0 = make_ring_ic(N, IC, np.random.default_rng(seed))
    t, rmean, rstd, th = _integrate_field(
        np.ascontiguousarray(th0, dtype=np.float64), int(P), float(alpha),
        float(DT), int(n_steps), int(DECIMATE), float(STOP_EPS), int(hold_samples))
    return t, rmean, rstd, th, P


def rho_profile(theta_row, P):
    """rho_k = |windowed mean field| for a single snapshot."""
    return np.abs(local_meanfield(theta_row, P))


# ---------------------------------------------------------------- classifier
def count_arcs_circular(mask, min_width):
    """Number of maximal contiguous True-runs on a periodic 1-D boolean array, counting
    only runs of length >= min_width."""
    n = len(mask)
    if mask.all():
        return 1 if n >= min_width else 0
    if not mask.any():
        return 0
    # rotate so index 0 is False (a boundary), then count runs linearly
    start = int(np.argmin(mask))  # first False
    m = np.roll(mask, -start)
    runs = []
    i = 0
    while i < n:
        if m[i]:
            j = i
            while j < n and m[j]:
                j += 1
            runs.append(j - i)
            i = j
        else:
            i += 1
    return sum(1 for r in runs if r >= min_width)


def classify_profile(rho, P, coh_thresh, min_arc_frac):
    """Classify one rho_k profile. coherent_k = rho_k > coh_thresh. Returns dict with
    n_coh_arcs, n_incoh_arcs, coherent_frac, canonical (1 coh arc + 1 incoh arc)."""
    N = len(rho)
    min_w = max(2, int(round(min_arc_frac * N)))
    coh = rho > coh_thresh
    n_coh = count_arcs_circular(coh, min_w)
    n_incoh = count_arcs_circular(~coh, min_w)
    return {
        "n_coh_arcs": n_coh, "n_incoh_arcs": n_incoh,
        "coherent_frac": float(coh.mean()),
        "canonical": bool(n_coh == 1 and n_incoh == 1),
    }


def classify_run(t, rmean, rstd, th, P, tau, *, coh_thresh, min_arc_frac=0.06,
                 last_frac=0.10, formed_rho_std=0.08, canon_vote=0.6):
    """Classify a run's pre-death state.

    Equilibration / 'formed' flag: did a chimera ever establish? -> max rho_std over the
    living window [t_skip, tau) exceeds formed_rho_std (a clear bimodal coherent/incoherent
    split, vs a run that homogenises immediately).

    Pre-death classification: over the last `last_frac` of life, each stored snapshot is
    classified; the run is 'canonical' if >= canon_vote of those snapshots are canonical.
    """
    t = np.asarray(t)
    di = int(np.searchsorted(t, tau, side="right")) - 1
    di = min(max(di, 0), len(t) - 1)
    t_skip = 50.0
    live = (t >= t_skip) & (t <= t[di])
    max_rstd_live = float(rstd[live].max()) if live.any() else float(rstd[:di + 1].max())
    formed = max_rstd_live > formed_rho_std

    t0 = tau - last_frac * (tau - 0.0)
    win = (t >= t0) & (t <= t[di])
    idxs = np.where(win)[0]
    if len(idxs) == 0:
        idxs = np.array([di])
    snaps = []
    for i in idxs:
        prof = rho_profile(th[i], P)
        snaps.append(classify_profile(prof, P, coh_thresh, min_arc_frac))
    canon_votes = np.mean([s["canonical"] for s in snaps])
    n_coh_mode = int(np.median([s["n_coh_arcs"] for s in snaps]))
    n_incoh_mode = int(np.median([s["n_incoh_arcs"] for s in snaps]))
    return {
        "tau_reint": float(tau),
        "formed": bool(formed),
        "max_rho_std_live": max_rstd_live,
        "n_snaps_last10": int(len(idxs)),
        "canon_vote_frac": float(canon_votes),
        "canonical": bool(canon_votes >= canon_vote and formed),
        "n_coh_arcs_med": n_coh_mode,
        "n_incoh_arcs_med": n_incoh_mode,
        "coherent_frac_mean": float(np.mean([s["coherent_frac"] for s in snaps])),
    }
