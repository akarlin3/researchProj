"""Diagnostic: plot r1(t), r2(t) trajectories across A at larger N to SEE whether a
sustained intermediate-r chimera plateau exists at beta=0.1. Not a checkpoint deliverable.
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .config_io import load_config, make_params
from .model import make_initial_conditions, integrate

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results")


def main():
    cfg = load_config(os.path.join(ROOT, "config.yaml"))
    N = 256
    A_list = [0.10, 0.15, 0.20, 0.25, 0.30, 0.40]
    seed = cfg["seed_base"]
    Tshow = 1500.0

    fig, axes = plt.subplots(len(A_list), 1, figsize=(11, 2.0 * len(A_list)), sharex=True)
    for ax, A in zip(axes, A_list):
        p = make_params(A, cfg["beta"])
        rng = np.random.default_rng(seed)
        th0 = make_initial_conditions(N, cfg["ic"], rng)
        res = integrate(theta0=th0, p=p, N=N, dt=cfg["dt"], T_max=Tshow,
                        decimate=cfg["output"]["decimate"], seed=seed)
        ax.plot(res.t, res.r1, lw=0.6, color="C0", label="r1")
        ax.plot(res.t, res.r2, lw=0.6, color="C3", label="r2")
        m = res.t >= 200
        ax.set_ylabel(f"A={A}\nr", fontsize=8)
        ax.set_ylim(0, 1.03)
        ax.axhline(float(np.minimum(res.r1, res.r2)[m].mean()), color="C3", ls=":", lw=0.7)
        ax.text(0.99, 0.06,
                f"r2 mean[200,{int(Tshow)}]={np.minimum(res.r1,res.r2)[m].mean():.3f}  "
                f"std={np.minimum(res.r1,res.r2)[m].std():.3f}",
                transform=ax.transAxes, ha="right", fontsize=8)
    axes[0].legend(loc="lower right", fontsize=8, ncol=2)
    axes[0].set_title(f"Trajectory scan  beta={cfg['beta']}, N={N}, seed={seed}")
    axes[-1].set_xlabel("t")
    fig.tight_layout()
    out = os.path.join(RESULTS, "diag_traj.png")
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print("saved", out)


if __name__ == "__main__":
    main()
