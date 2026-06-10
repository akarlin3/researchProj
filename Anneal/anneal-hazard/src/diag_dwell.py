"""Pre-flight (CP3 check 2): visualise rho_std(t) approach-to-collapse on several pilot
runs, so the near-collapse 'dwell' band is set from data (not guessed), and confirm the
decimated resolution (0.5 t.u.) is fine vs the time spent descending into collapse.
"""
from __future__ import annotations

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from src.config_io import load_config
from src.ring_model import make_ring_ic
from src.ring_fast import integrate_ring_fast
from src.ring_detector import detect_death_ring, dwell_stat_ring

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results")

BAND = (0.04, 0.10)  # candidate near-collapse band for rho_std


def main():
    cfg = load_config(os.path.join(ROOT, "config.yaml"))
    beta, N = 0.12, 64
    P = max(1, int(round(cfg["r"] * N)))
    dt, T_max, dt_hold = cfg["dt"], cfg["T_max"], cfg["dt_hold"]
    eps = cfg["eps_std"]; dec = cfg["output"]["decimate"]

    fig, axes = plt.subplots(3, 2, figsize=(13, 9), sharey=True)
    descent_times = []
    for ax, k in zip(axes.flat, range(6)):
        seed = cfg["seed_base"] + k
        rng = np.random.default_rng(seed)
        th0 = make_ring_ic(N, cfg["ic"], rng)
        res = integrate_ring_fast(th0, N, P, beta, dt, T_max, dec, seed,
                                  stop_eps=min(0.03, eps), dt_hold=dt_hold)
        tau, ev = detect_death_ring(res.t, res.rho_std, eps, dt_hold, T_max)
        ax.plot(res.t, res.rho_std, lw=0.7, color="C0")
        ax.axhspan(BAND[0], BAND[1], color="orange", alpha=0.25)
        ax.axhline(eps, color="k", ls=":", lw=0.8)
        if ev:
            ax.axvline(tau, color="C3", ls="--", lw=0.8)
        ax.set_xlim(max(0, (tau if ev else T_max) - 250), (tau if ev else T_max) + 30)
        ax.set_title(f"seed {seed}  τ={tau:.0f} ev={ev}", fontsize=9)
        ax.set_ylim(0, 0.25)
        # terminal committed-descent dwell (final continuous rho_std<band_hi run ending at death)
        if ev:
            dwell = dwell_stat_ring(res.t, res.rho_std, tau, ev, band_hi=BAND[1])
            descent_times.append(dwell)
            ax.text(0.02, 0.9, f"terminal dwell≈{dwell:.1f} t.u.\n({dwell/(dec*dt):.0f} samples)",
                    transform=ax.transAxes, fontsize=8,
                    bbox=dict(boxstyle="round", fc="white", alpha=0.8))
    for ax in axes[-1]:
        ax.set_xlabel("t")
    for ax in axes[:, 0]:
        ax.set_ylabel("rho_std")
    fig.suptitle(f"Approach to collapse  β={beta} N={N}  band={BAND}  decimated dt={dec*dt} t.u.")
    fig.tight_layout()
    out = os.path.join(RESULTS, "diag_dwell.png")
    fig.savefig(out, dpi=120); plt.close(fig)
    print(f"band={BAND}  decimated_dt={dec*dt} t.u.")
    print(f"dwell-in-band (t.u.) across 6 runs: {[round(d,1) for d in descent_times]}")
    if descent_times:
        print(f"min dwell={min(descent_times):.1f} t.u. = {min(descent_times)/(dec*dt):.0f} decimated samples")
    print("saved", out)


if __name__ == "__main__":
    main()
