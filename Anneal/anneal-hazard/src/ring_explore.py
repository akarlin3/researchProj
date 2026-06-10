"""Diagnostic: confirm a chimera forms in the ring model, visualise it (space-time local
order parameter + a phase snapshot), and read off the rho_std level for detector tuning.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from src.ring_model import local_meanfield, make_ring_ic

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results")


def integrate_with_field(theta0, N, P, beta, dt, T_max, store_dt):
    """Integrate and store the full local-coherence field rho_k(t) and theta snapshots."""
    alpha = np.pi / 2.0 - beta
    n_steps = int(round(T_max / dt))
    store_every = max(1, int(round(store_dt / dt)))
    x = theta0.astype(float).copy()

    def deriv(th):
        Wbar = local_meanfield(th, P)
        return -np.imag(np.exp(1j * (th + alpha)) * np.conj(Wbar))

    ts, rho_field, theta_field = [], [], []
    for step in range(n_steps + 1):
        if step % store_every == 0:
            Wbar = local_meanfield(x, P)
            ts.append(step * dt)
            rho_field.append(np.abs(Wbar).copy())
            theta_field.append(np.mod(x + np.pi, 2 * np.pi) - np.pi)
        if step == n_steps:
            break
        k1 = deriv(x); k2 = deriv(x + 0.5 * dt * k1)
        k3 = deriv(x + 0.5 * dt * k2); k4 = deriv(x + dt * k3)
        x = x + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    return np.array(ts), np.array(rho_field), np.array(theta_field)


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    beta = float(sys.argv[2]) if len(sys.argv) > 2 else 0.11
    r = float(sys.argv[3]) if len(sys.argv) > 3 else 0.35
    T = float(sys.argv[4]) if len(sys.argv) > 4 else 2000.0
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 20260609
    P = max(1, int(round(r * N)))
    dt = 0.05

    rng = np.random.default_rng(seed)
    th0 = make_ring_ic(N, {"incoherent_frac": 0.5, "coherent_scale": 0.05}, rng)
    ts, rho, theta = integrate_with_field(th0, N, P, beta, dt, T, store_dt=2.0)

    rho_std = rho.std(axis=1)
    rho_min = rho.min(axis=1)
    rho_mean = rho.mean(axis=1)

    # characterise the chimera plateau (skip transient)
    msk = ts >= 100
    print(f"ring chimera probe: N={N} beta={beta} r={r} (P={P}) T={T} seed={seed}")
    print(f"  rho_std  over t>=100: mean={rho_std[msk].mean():.4f} "
          f"min={rho_std[msk].min():.4f} max={rho_std[msk].max():.4f}")
    print(f"  rho_min  over t>=100: mean={rho_min[msk].mean():.4f} "
          f"(incoherent region depth)")
    print(f"  rho_mean over t>=100: mean={rho_mean[msk].mean():.4f}")
    print(f"  final rho_std={rho_std[-1]:.4f}  (near 0 => collapsed to sync)")

    fig, ax = plt.subplots(3, 1, figsize=(11, 9),
                           gridspec_kw={"height_ratios": [2, 1, 1]})
    im = ax[0].imshow(rho.T, aspect="auto", origin="lower",
                      extent=[ts[0], ts[-1], 0, N], cmap="viridis", vmin=0, vmax=1)
    ax[0].set_ylabel("oscillator index k")
    ax[0].set_title(f"Local coherence ρ_k(t)  —  ring  N={N}, β={beta}, r={r} (P={P}), seed={seed}")
    fig.colorbar(im, ax=ax[0], label="ρ_k")

    ax[1].plot(ts, rho_std, lw=0.8, color="C2", label="std_k ρ_k (chimera indicator)")
    ax[1].plot(ts, rho_min, lw=0.8, color="C1", alpha=0.7, label="min_k ρ_k")
    ax[1].set_ylabel("ρ spread"); ax[1].legend(fontsize=8, loc="right"); ax[1].set_ylim(0, None)

    # phase snapshot at a representative chimera time
    j = np.argmin(np.abs(ts - min(ts[-1], 100.0)))
    xk = np.linspace(-np.pi, np.pi, N, endpoint=False)
    ax[2].plot(xk, theta[j], ".", ms=3, color="C3")
    ax[2].set_ylabel(f"θ_k at t={ts[j]:.0f}"); ax[2].set_xlabel("x_k  (and: t for panels above)")
    ax[2].set_ylim(-np.pi, np.pi)

    fig.tight_layout()
    out = os.path.join(RESULTS, f"ring_probe_N{N}_b{beta}_r{r}_s{seed}.png")
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print("saved", out)


if __name__ == "__main__":
    main()
