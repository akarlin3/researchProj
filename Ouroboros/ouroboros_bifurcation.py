import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

# Import model RHS and Laplacian from simulator
try:
    from ouroboros_sim import stroma_rhs
except ImportError:
    def laplacian(u, dx):
        lap = np.zeros_like(u)
        lap[1:-1] = (u[2:] - 2*u[1:-1] + u[:-2]) / (dx**2)
        lap[0] = (2*u[1] - 2*u[0]) / (dx**2)
        lap[-1] = (2*u[-2] - 2*u[-1]) / (dx**2)
        return lap

    def stroma_rhs(t, y, Nx, dx, Dp, Dc, Dn, gam_p, Sp, gam_c, Sc, rho, gam_n):
        p = y[:Nx]
        c = y[Nx:2*Nx]
        n = y[2*Nx:]
        p_clipped = np.clip(p, 0.0, None)
        c_clipped = np.clip(c, 0.0, None)
        n_clipped = np.clip(n, 0.0, 1.5)
        lap_p = laplacian(p_clipped, dx)
        lap_c = laplacian(c_clipped, dx)
        lap_n = laplacian(n_clipped, dx)
        dp_dt = Dp * lap_p - gam_p * p_clipped + Sp * n_clipped * (1.0 - p_clipped)
        dc_dt = Dc * lap_c - gam_c * c_clipped * (p_clipped / (1.0 + p_clipped)) + Sc * n_clipped
        dn_dt = Dn * lap_n + rho * n_clipped * (1.0 - n_clipped) * (c_clipped / (1.0 + c_clipped)) - gam_n * n_clipped * p_clipped
        return np.concatenate([dp_dt, dc_dt, dn_dt])

def run_pde_sim(rho, gam_n, Dn, Nx=100, L=10.0, T=80.0, dt=0.05, seed=42):
    """Integrate the PDE and return the time series of the spatial mean of pressure p"""
    np.random.seed(seed)
    dx = L / (Nx - 1)
    Nt = int(T / dt) + 1
    t_eval = np.linspace(0, T, Nt)
    
    x = np.linspace(0, L, Nx)
    # Standard initial conditions
    p0 = np.exp(-(x - L/2)**2 / 2.0)
    c0 = 1.0 - 0.5 * np.exp(-(x - L/2)**2 / 2.0)
    n0 = 0.1 + 0.2 * np.sin(np.pi * x / L)
    y0 = np.concatenate([p0, c0, n0])
    
    # Parameters
    Dp, Dc = 0.05, 0.1
    gam_p, Sp = 0.05, 0.1
    gam_c, Sc = 0.2, 0.3
    
    sol = solve_ivp(
        stroma_rhs,
        [0, T],
        y0,
        t_eval=t_eval,
        method='Radau',
        args=(Nx, dx, Dp, Dc, Dn, gam_p, Sp, gam_c, Sc, rho, gam_n)
    )
    
    # Pressure is the first Nx elements
    p_sol = sol.y[:Nx, :].T  # shape (Nt, Nx)
    p_mean = np.mean(p_sol, axis=1) # spatial mean over time
    
    return sol.t, p_mean

def classify_dynamics(t, ts, dt, t_trans=40.0):
    """Classify the long-time behavior of the time series"""
    # Discard transient
    idx_start = np.searchsorted(t, t_trans)
    ts_clean = ts[idx_start:]
    
    variance = np.var(ts_clean)
    
    if variance < 1e-6:
        return 'fixed point', [], variance
        
    # Analyze frequency content
    n_pts = len(ts_clean)
    fft_vals = np.fft.rfft(ts_clean - np.mean(ts_clean))
    fft_freqs = np.fft.rfftfreq(n_pts, d=dt)
    fft_power = np.abs(fft_vals)**2
    
    # Find peaks in FFT power
    peaks, _ = find_peaks(fft_power, prominence=np.max(fft_power)*0.05)
    
    # Find local maxima in the time series for bifurcation diagram
    ts_peaks, _ = find_peaks(ts_clean, distance=5)
    peak_vals = ts_clean[ts_peaks]
    
    # If no local peaks found, use the final values (or all values if it is chaotic)
    if len(peak_vals) == 0:
        peak_vals = np.array([ts_clean[-1]])
        
    if len(peaks) <= 2:
        return 'periodic', peak_vals, variance
    elif len(peaks) <= 5:
        return 'quasiperiodic', peak_vals, variance
    else:
        return 'candidate-chaotic', peak_vals, variance

def main():
    print("="*60)
    print("OUROBOROS BIFURCATION SCAN & DYNAMICS CLASSIFICATION")
    print("="*60)
    
    # Slice choice: we sweep gam_n along the most unstable region boundary
    # found in Checkpoint 1 (rho=0.1, Dn=0.001)
    rho_fixed = 0.1
    Dn_fixed = 0.001
    gam_n_vals = np.linspace(0.1, 5.0, 10)
    
    bif_x = []
    bif_y = []
    
    print(f"Running 1D sweep over gam_n in [0.1, 5.0] at rho={rho_fixed}, Dn={Dn_fixed}...")
    
    for gn in gam_n_vals:
        t, p_mean = run_pde_sim(rho_fixed, gn, Dn_fixed, T=300.0, dt=0.05)
        regime, peaks, var = classify_dynamics(t, p_mean, dt=0.05, t_trans=250.0)
        
        print(f"  gam_n={gn:.2f} -> regime: {regime:<15} (var: {var:.2e}, peaks detected: {len(peaks)})")
        
        # If it's a fixed point, plot the final value
        if regime == 'fixed point':
            bif_x.append(gn)
            bif_y.append(p_mean[-1])
        else:
            for pk in peaks:
                bif_x.append(gn)
                bif_y.append(pk)
                
    # Plot Bifurcation Diagram
    plt.figure(figsize=(8, 5))
    plt.scatter(bif_x, bif_y, color='red', s=15, alpha=0.8, edgecolors='none', label='Attractor State(s)')
    plt.title(f"OUROBOROS Bifurcation Scan (rho={rho_fixed}, Dn={Dn_fixed})", fontsize=12)
    plt.xlabel("Vessel Regression Rate (gam_n)", fontsize=10)
    plt.ylabel("Pressure Spatial Mean (p_mean)", fontsize=10)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()
    plt.tight_layout()
    os.makedirs('figures', exist_ok=True)
    plt.savefig('figures/bifurcation_diagram.png', dpi=150)
    plt.close()
    print("\nSaved bifurcation diagram to figures/bifurcation_diagram.png")
    
if __name__ == '__main__':
    main()
