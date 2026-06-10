"""Diagnostic: find the metastable regime — (beta, N) where the chimera FORMS cleanly
yet DIES within T_max with a spread of lifetimes, and lifetime grows with N.
Reports n_died/n, median tau, and plateau r_incoh per (beta, N). Not a deliverable.
"""
from __future__ import annotations

import os
import sys

import numpy as np

from .config_io import load_config, make_params
from .model import make_initial_conditions, integrate
from .detector import detect_death, precollapse_stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
T_SKIP = 50.0


def main():
    cfg = load_config(os.path.join(ROOT, "config.yaml"))
    A = float(sys.argv[1]) if len(sys.argv) > 1 else 0.20
    betas = [float(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [0.07, 0.08, 0.09, 0.095]
    Ns = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else [64, 128, 256]
    nseed = int(sys.argv[4]) if len(sys.argv) > 4 else 5
    seeds = [cfg["seed_base"] + k for k in range(nseed)]

    print(f"metastability  A={A}  eps={cfg['eps']}  dt_hold={cfg['dt_hold']}  "
          f"T_max={cfg['T_max']}  dt={cfg['dt']}  nseed={nseed}")
    print(f"{'beta':>5} {'N':>4} | {'n_died':>6} | {'median tau':>10} | {'tau range':>16} | {'r_incoh plateau':>16}")
    print("-" * 78)
    for beta in betas:
        for N in Ns:
            p = make_params(A, beta)
            taus, evs, rincs = [], [], []
            for seed in seeds:
                rng = np.random.default_rng(seed)
                th0 = make_initial_conditions(N, cfg["ic"], rng)
                res = integrate(theta0=th0, p=p, N=N, dt=cfg["dt"], T_max=cfg["T_max"],
                                decimate=cfg["output"]["decimate"], seed=seed,
                                stop_eps=cfg["eps"], dt_hold=cfg["dt_hold"])
                tau, ev = detect_death(res.t, res.r1, res.r2, cfg["eps"], cfg["dt_hold"], cfg["T_max"])
                st = precollapse_stats(res.t, res.r1, res.r2, tau, ev, T_SKIP, cfg["T_max"])
                taus.append(tau); evs.append(ev); rincs.append(st["r_incoh_mean"])
            ndied = int(sum(evs))
            med = float(np.median(taus))
            tmin, tmax = float(np.min(taus)), float(np.max(taus))
            rmean = float(np.nanmean(rincs))
            print(f"{beta:5.3f} {N:4d} | {ndied:2d}/{nseed:<3d} | {med:10.0f} | "
                  f"[{tmin:6.0f},{tmax:6.0f}] | {rmean:.3f}")
        print()


if __name__ == "__main__":
    main()
