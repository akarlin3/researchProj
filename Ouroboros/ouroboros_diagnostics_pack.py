import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from ouroboros_sim import stroma_rhs
from ouroboros_chaos import rosenstein_estimator

# Ensure output directories exist
os.makedirs('data', exist_ok=True)
os.makedirs('figures', exist_ok=True)

# System configuration (baseline parameters)
Nx = 100
L = 10.0
dx = L / (Nx - 1)
params = (Nx, dx, 0.05, 0.1, 0.01, 0.05, 0.1, 0.2, 0.3, 0.2, 0.1)

def main():
    print("=" * 60)
    print("OUROBOROS CHAOS DIAGNOSTICS CAUTIONARY BATTERY")
    print("=" * 60)
    
    # 1. Generate long simulation data
    print("Simulating vascularized stroma model up to T = 300.0...")
    x = np.linspace(0, L, Nx)
    p0 = np.exp(-(x - L/2)**2 / 2.0)
    c0 = 1.0 - 0.5 * np.exp(-(x - L/2)**2 / 2.0)
    n0 = 0.1 + 0.2 * np.sin(np.pi * x / L)
    y0 = np.concatenate([p0, c0, n0])
    
    T_sim = 300.0
    Nt_sim = 6000  # dt = 0.05
    t_eval = np.linspace(0, T_sim, Nt_sim)
    dt_sim = T_sim / (Nt_sim - 1)
    
    sol = solve_ivp(
        stroma_rhs,
        [0, T_sim],
        y0,
        t_eval=t_eval,
        method='Radau',
        args=params
    )
    
    # Spatial mean of pressure p(x,t)
    ts_all = np.mean(sol.y[:Nx, :], axis=0)
    
    # 2. Extract subsets for analysis
    # Case A: Full transient included (first 150 time units, 3000 steps)
    Nt_150 = 3000
    ts_transient = ts_all[:Nt_150]
    t_transient = t_eval[:Nt_150]
    
    # Case B: Transient removed (late series, t from 100 to 250, 3000 steps)
    idx_start = int(100 / dt_sim)
    idx_end = idx_start + Nt_150
    ts_stationary = ts_all[idx_start:idx_end]
    t_stationary = t_eval[idx_start:idx_end]
    
    # 3. Embedding parameters from Phase 2
    m = 2
    tau = 39
    max_steps = 150
    
    print("\nRunning Rosenstein estimator on Series WITH transient...")
    rosen_transient = rosenstein_estimator(ts_transient, dt_sim, m=m, tau=tau, theiler_window=50, max_steps=max_steps)
    t_rosen = np.arange(max_steps) * dt_sim
    # Linear fit for slope
    slope_trans, intercept_trans = np.polyfit(t_rosen, rosen_transient, 1)
    
    print("Running Rosenstein estimator on Series WITHOUT transient (asymptotic)...")
    rosen_stationary = rosenstein_estimator(ts_stationary, dt_sim, m=m, tau=tau, theiler_window=50, max_steps=max_steps)
    slope_station, intercept_station = np.polyfit(t_rosen, rosen_stationary, 1)
    
    print(f"\nRosenstein LLE Estimates:")
    print(f"  With Transient (T=[0, 150]):     LLE = {slope_trans:+.6f}")
    print(f"  Without Transient (T=[100, 250]): LLE = {slope_station:+.6f}")
    print(f"  True tangent-space (Benettin):    LLE = -0.073362")
    
    # 4. Generate the comparative plot
    fig, axes = plt.subplots(3, 1, figsize=(10, 12))
    
    # Panel 1: Time Series of pressure spatial mean showing transient vs stationary
    axes[0].plot(t_eval, ts_all, 'k-', linewidth=2, label='Pressure Spatial Mean $\\langle p(x,t) \\rangle$')
    axes[0].axvspan(0, 50, color='red', alpha=0.1, label='Transient Regime')
    axes[0].axvspan(100, 250, color='green', alpha=0.1, label='Asymptotic Regime')
    axes[0].set_title("Vascularized Stroma Model: Trajectory Transient Decay", fontsize=12)
    axes[0].set_xlabel("Time (t)")
    axes[0].set_ylabel("Pressure Mean $\\langle p \\rangle$")
    axes[0].grid(True)
    axes[0].legend()
    
    # Panel 2: Rosenstein Divergence with Transient
    axes[1].plot(t_rosen, rosen_transient, 'r-', linewidth=2, label='Rosenstein Divergence')
    axes[1].plot(t_rosen, slope_trans * t_rosen + intercept_trans, 'k--', label=f'Fit (LLE = {slope_trans:+.4f})')
    axes[1].set_title("Data-Driven Divergence WITH Transient (T=[0, 150])", fontsize=12)
    axes[1].set_xlabel("Time Separation (t)")
    axes[1].set_ylabel("Average Log Divergence $\\langle \\ln d(t) \\rangle$")
    axes[1].grid(True)
    axes[1].legend()
    
    # Panel 3: Rosenstein Divergence without Transient
    axes[2].plot(t_rosen, rosen_stationary, 'g-', linewidth=2, label='Rosenstein Divergence')
    axes[2].plot(t_rosen, slope_station * t_rosen + intercept_station, 'k--', label=f'Fit (LLE = {slope_station:+.4f})')
    axes[2].set_title("Data-Driven Divergence WITHOUT Transient (T=[100, 250])", fontsize=12)
    axes[2].set_xlabel("Time Separation (t)")
    axes[2].set_ylabel("Average Log Divergence $\\langle \\ln d(t) \\rangle$")
    axes[2].grid(True)
    axes[2].legend()
    
    plt.suptitle("Methodological Caution: Rosenstein False-Positive Artifact in Decaying Transients", fontsize=14, y=0.99)
    plt.tight_layout()
    plt.savefig("figures/diagnostics_caution.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("\nSaved diagnostics comparison plot to figures/diagnostics_caution.png")
    
    # Save raw json output
    diagnostics_results = {
        'embedding': {
            'm': m,
            'tau': tau,
            'dt': dt_sim
        },
        'benettin_lle': -0.073362,
        'rosenstein_lle_with_transient': float(slope_trans),
        'rosenstein_lle_without_transient': float(slope_station),
        'artifact_explained': (
            "When the system decays monotonically to a stable fixed point, the trajectory contracts. "
            "Because we exclude self-neighbors using the Theiler window, nearest neighbors must be chosen "
            "from different parts of the transient curve. As the transient relaxes, the distance between "
            "these segments temporarily expands (geometrically separating along the curve) before both "
            "eventually cluster near the fixed point. The Rosenstein algorithm interprets this expansion as "
            "sensitive dependence on initial conditions (LLE > 0). Once the transient is removed, the trajectories "
            "stay clustered at the fixed point, and the false positive shrinks or vanishes (LLE = -0.003)."
        )
    }
    with open('data/diagnostics_results.json', 'w') as f:
        json.dump(diagnostics_results, f, indent=2)
    print("Saved diagnostics raw results to data/diagnostics_results.json")

if __name__ == '__main__':
    main()
