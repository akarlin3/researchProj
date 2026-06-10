"""Diagnostic: chimera lifetime vs N in the ring model (the finite-N chaotic-transient
test). Reports n_died, median tau, tau range, and the sync-floor of rho_std per (N).
"""
from __future__ import annotations

import os
import sys

import numpy as np

from src.ring_model import make_ring_ic
from src.ring_fast import integrate_ring_fast as integrate_ring
from src.ring_detector import detect_death_ring, precollapse_stats_ring

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def main():
    beta = float(sys.argv[1]) if len(sys.argv) > 1 else 0.08
    r = float(sys.argv[2]) if len(sys.argv) > 2 else 0.15
    Ns = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else [32, 48, 64, 96, 128, 192]
    nseed = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    eps = float(sys.argv[5]) if len(sys.argv) > 5 else 0.04
    T_max = float(sys.argv[6]) if len(sys.argv) > 6 else 3000.0
    dt_hold = 50.0
    dt = 0.05
    seeds = [20260609 + k for k in range(nseed)]
    ic = {"incoherent_frac": 0.5, "coherent_scale": 0.05}

    print(f"ring lifetime vs N  beta={beta} r={r} eps={eps} dt_hold={dt_hold} "
          f"T_max={T_max} dt={dt} nseed={nseed}")
    print(f"{'N':>4} {'P':>4} | {'n_died':>7} | {'median tau':>10} | {'tau range':>15} | "
          f"{'chimera rho_std':>15} | {'sync floor':>10}")
    print("-" * 86)
    for N in Ns:
        P = max(1, int(round(r * N)))
        taus, evs, plats, floors = [], [], [], []
        for seed in seeds:
            rng = np.random.default_rng(seed)
            th0 = make_ring_ic(N, ic, rng)
            res = integrate_ring(th0, N, P, beta, dt, T_max, decimate=10, seed=seed,
                                 stop_eps=eps, dt_hold=dt_hold)
            tau, ev = detect_death_ring(res.t, res.rho_std, eps, dt_hold, T_max)
            st = precollapse_stats_ring(res.t, res.rho_std, res.rho_mean, tau, ev, 50.0, T_max)
            taus.append(tau); evs.append(ev); plats.append(st["rho_std_mean"])
            # sync floor: rho_std in the last bit of a dead run (post-collapse)
            if ev:
                floors.append(float(res.rho_std[-1]))
        ndied = int(sum(evs))
        med = float(np.median(taus))
        tmin, tmax = float(np.min(taus)), float(np.max(taus))
        plat = float(np.nanmean(plats))
        floor = float(np.mean(floors)) if floors else float("nan")
        print(f"{N:4d} {P:4d} | {ndied:2d}/{nseed:<4d} | {med:10.0f} | "
              f"[{tmin:6.0f},{tmax:6.0f}] | {plat:15.3f} | {floor:10.4f}")


if __name__ == "__main__":
    main()
