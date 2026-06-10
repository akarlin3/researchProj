"""Diagnostic: locate the ring chimera in the (beta, r) plane. Chimera signature =
sustained intermediate rho_std (spatial spread of local coherence): not 0 (sync) and
persisting to the end of the window. Reports mean rho_std over a late window + final value.
"""
from __future__ import annotations

import os
import sys

import numpy as np

from src.ring_model import make_ring_ic, integrate_ring

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 128
    T = float(sys.argv[2]) if len(sys.argv) > 2 else 1500.0
    betas = [float(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else [0.02, 0.05, 0.08, 0.11]
    rs = [float(x) for x in sys.argv[4].split(",")] if len(sys.argv) > 4 else [0.10, 0.15, 0.20, 0.25, 0.30, 0.35]
    seed = 20260609
    dt = 0.05
    win = 300.0
    ic = {"incoherent_frac": 0.5, "coherent_scale": 0.05}

    print(f"ring (beta,r) chimera map  N={N}  T={T}  late-window>={win}  seed={seed}")
    print("cell = mean rho_std [final rho_std] ; '*' = sustained chimera candidate")
    header = "beta\\r | " + " ".join(f"{r:>11.2f}" for r in rs)
    print(header); print("-" * len(header))
    for beta in betas:
        cells = []
        for r in rs:
            P = max(1, int(round(r * N)))
            rng = np.random.default_rng(seed)
            th0 = make_ring_ic(N, ic, rng)
            res = integrate_ring(th0, N, P, beta, dt, T, decimate=10, seed=seed)
            m = res.t >= win
            mean_std = float(res.rho_std[m].mean())
            fin = float(res.rho_std[-1])
            chim = (mean_std > 0.03 and fin > 0.02)
            cells.append(f"{mean_std:.3f}[{fin:.3f}]{'*' if chim else ' '}")
        print(f"{beta:5.3f} | " + " ".join(f"{c:>11}" for c in cells))
    print("-" * len(header))


if __name__ == "__main__":
    main()
