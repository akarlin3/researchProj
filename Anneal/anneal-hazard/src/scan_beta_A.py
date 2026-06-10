"""Diagnostic: locate the chimera in the (beta, A) plane at N=256.
Chimera signature = incoherent population settles to an INTERMEDIATE, sustained r
(not 0, not 1). Reports mean/std of min(r1,r2) over a late window. Not a deliverable.
"""
from __future__ import annotations

import os
import sys

import numpy as np

from .config_io import load_config, make_params
from .model import make_initial_conditions, integrate

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def main():
    cfg = load_config(os.path.join(ROOT, "config.yaml"))
    N = 256
    T = 1500.0
    win_lo, win_hi = 400.0, 1500.0
    betas = [0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.25, 0.30]
    As = [0.05, 0.10, 0.20, 0.30, 0.40, 0.50]
    if len(sys.argv) > 1:
        betas = [float(x) for x in sys.argv[1].split(",")]
    if len(sys.argv) > 2:
        As = [float(x) for x in sys.argv[2].split(",")]
    seed = cfg["seed_base"]

    print(f"(beta,A) chimera map  N={N}  T={T}  window=[{win_lo},{win_hi}]  seed={seed}")
    print("cell = mean(min(r1,r2)) over window ; '*' if intermediate & sustained (chimera candidate)")
    header = "beta\\A | " + " ".join(f"{A:>6.2f}" for A in As)
    print(header)
    print("-" * len(header))
    for beta in betas:
        cells = []
        for A in As:
            p = make_params(A, beta)
            rng = np.random.default_rng(seed)
            th0 = make_initial_conditions(N, cfg["ic"], rng)
            res = integrate(theta0=th0, p=p, N=N, dt=cfg["dt"], T_max=T,
                            decimate=cfg["output"]["decimate"], seed=seed)
            m = (res.t >= win_lo) & (res.t <= win_hi)
            rinc = np.minimum(res.r1, res.r2)[m]
            mean = float(rinc.mean()); std = float(rinc.std())
            chim = (0.30 < mean < 0.95)
            cells.append(f"{mean:.2f}{'*' if chim else ' '}")
        print(f"{beta:5.2f} | " + " ".join(f"{c:>6}" for c in cells))
    print("-" * len(header))


if __name__ == "__main__":
    main()
