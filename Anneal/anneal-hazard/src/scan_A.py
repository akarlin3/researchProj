"""Diagnostic scan (not a checkpoint deliverable): where does a clean, metastable
chimera live at beta=0.1, N=128? Reports lifetime + plateau quality vs A across seeds.
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
    N = 128
    A_grid = [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70]
    seeds = [cfg["seed_base"] + k for k in range(3)]
    if len(sys.argv) > 1:
        A_grid = [float(x) for x in sys.argv[1].split(",")]

    print(f"scan  beta={cfg['beta']}  N={N}  eps={cfg['eps']}  dt_hold={cfg['dt_hold']}  "
          f"T_max={cfg['T_max']}  dt={cfg['dt']}  seeds={seeds}")
    print(f"{'A':>5} | {'tau (per seed)':>34} | {'med tau':>8} | {'r_incoh plateau mean±std (range)':>34}")
    print("-" * 100)

    for A in A_grid:
        p = make_params(A, cfg["beta"])
        taus, events, plateaus = [], [], []
        for seed in seeds:
            rng = np.random.default_rng(seed)
            th0 = make_initial_conditions(N, cfg["ic"], rng)
            res = integrate(theta0=th0, p=p, N=N, dt=cfg["dt"], T_max=cfg["T_max"],
                            decimate=cfg["output"]["decimate"], seed=seed,
                            stop_eps=cfg["eps"], dt_hold=cfg["dt_hold"])
            tau, ev = detect_death(res.t, res.r1, res.r2, cfg["eps"], cfg["dt_hold"], cfg["T_max"])
            st = precollapse_stats(res.t, res.r1, res.r2, tau, ev, T_SKIP, cfg["T_max"])
            taus.append(tau); events.append(ev); plateaus.append(st)
        med = float(np.median(taus))
        tau_strs = ", ".join(f"{t:6.0f}{'' if e else 'c'}" for t, e in zip(taus, events))
        pm = float(np.mean([s["r_incoh_mean"] for s in plateaus]))
        ps = float(np.mean([s["r_incoh_std"] for s in plateaus]))
        plo = float(np.min([s["r_incoh_min"] for s in plateaus]))
        phi = float(np.max([s["r_incoh_max"] for s in plateaus]))
        print(f"{A:5.2f} | {tau_strs:>34} | {med:8.0f} | "
              f"{pm:.3f}±{ps:.3f} [{plo:.2f},{phi:.2f}]")
    print("-" * 100)
    print("note: 'c' suffix = censored (survived to T_max). 'plateau' = min(r1,r2) over [50, death).")


if __name__ == "__main__":
    main()
