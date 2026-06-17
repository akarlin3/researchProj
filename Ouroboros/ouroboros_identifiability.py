import os
import json
import numpy as np
import matplotlib.pyplot as plt
from ouroboros_fractional_sindy import gl_weights, gl_derivative_time, build_fractional_pde_model
from ouroboros_model_select import select_temporal_order

# Set random seed for reproducibility
np.random.seed(42)

# Grid parameters
Nx = 50
Nt = 500
T = 5.0
L = 10.0
dt = T / (Nt - 1)
dx = L / (Nx - 1)

def linear_rhs(t, y, Nx):
    """
    A simple, well-posed linear RHS for testing fractional time derivatives.
    Couples the three variables linearly:
    dp/dt = -0.5 * p + 0.2 * c
    dc/dt = -0.3 * c + 0.1 * n
    dn/dt = -0.2 * n
    """
    p = y[:Nx]
    c = y[Nx:2*Nx]
    n = y[2*Nx:]
    
    dp = -0.5 * p + 0.2 * c
    dc = -0.3 * c + 0.1 * n
    dn = -0.2 * n
    return np.concatenate([dp, dc, dn])

def solve_fractional_system(alpha, Nt, dt, Nx, y0):
    """
    Solves D^alpha_t u = RHS(u) using the explicit Grünwald-Letnikov scheme.
    """
    w = gl_weights(alpha, Nt)
    y_sol = np.zeros((Nt, 3 * Nx))
    y_sol[0] = y0
    for k in range(1, Nt):
        rhs_val = linear_rhs(0.0, y_sol[k-1], Nx)
        history = np.zeros(3 * Nx)
        for j in range(1, k + 1):
            history += w[j] * y_sol[k - j]
        y_sol[k] = rhs_val * (dt**alpha) - history
    
    # Reshape to (Nt, Nx, 3)
    p = y_sol[:, :Nx]
    c = y_sol[:, Nx:2*Nx]
    n = y_sol[:, 2*Nx:]
    return np.stack([p, c, n], axis=-1)

def add_noise(u, snr_db):
    """
    Adds Gaussian measurement noise to the state variables.
    """
    if snr_db is None or np.isinf(snr_db):
        return u.copy()
    
    u_noisy = u.copy()
    for i in range(u.shape[-1]):
        signal = u[..., i]
        signal_var = np.var(signal)
        if signal_var == 0:
            continue
        noise_var = signal_var * (10 ** (-snr_db / 10.0))
        noise = np.random.normal(0, np.sqrt(noise_var), size=signal.shape)
        u_noisy[..., i] = signal + noise
    return u_noisy

def main():
    print("=" * 60)
    print("OUROBOROS FRACTIONAL SINDY IDENTIFIABILITY (SENSITIVITY & NOISE)")
    print("=" * 60)
    
    # Spatial grid and initial conditions
    x = np.linspace(0, L, Nx)
    p0 = np.exp(-(x - L/2)**2 / 2.0)
    c0 = 1.0 - 0.5 * np.exp(-(x - L/2)**2 / 2.0)
    n0 = 0.1 + 0.2 * np.sin(np.pi * x / L)
    y0 = np.concatenate([p0, c0, n0])
    
    true_alphas = [0.5, 0.7, 0.9]
    candidates = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    frac_space_orders = (0.5, 1.5)  # Dummy/unused space candidates since we fit ODE features
    
    results = {}
    
    # 1. Clean Sensitivity Sweep
    print("\n--- Part 1: Clean Data Sensitivity Sweep ---")
    for alpha in true_alphas:
        print(f"\nSimulating system with true alpha_t = {alpha}...")
        u_clean = solve_fractional_system(alpha, Nt, dt, Nx, y0)
        
        # Run sweep
        sweep, best = select_temporal_order(u_clean, dt, dx, alpha_t_candidates=candidates, 
                                             frac_space_orders=frac_space_orders, threshold=0.01, k_start=20, u_clean=u_clean, true_alpha=alpha)
        
        selected_alpha = best['alpha_t']
        error = abs(selected_alpha - alpha)
        
        # Calculate margin over next best candidate
        r2_vals = [res['r2_avg'] for res in sweep]
        sorted_indices = np.argsort(r2_vals)[::-1]
        best_r2 = r2_vals[sorted_indices[0]]
        next_best_r2 = r2_vals[sorted_indices[1]] if len(sorted_indices) > 1 else -np.inf
        margin = best_r2 - next_best_r2
        
        print(f"True alpha: {alpha} | Selected alpha: {selected_alpha} | Error: {error:.4f} | R2 Margin: {margin:.4f}")
        
        results[str(alpha)] = {
            'clean': {
                'selected_alpha': float(selected_alpha),
                'error': float(error),
                'margin': float(margin),
                'r2_avg': float(best['r2_avg']),
                'sweep_r2s': {str(res['alpha_t']): float(res['r2_avg']) for res in sweep}
            },
            'noise_sweep': {}
        }
    
    # 2. Noise Robustness Sweep
    print("\n--- Part 2: Noise Robustness Sweep ---")
    snrs = [100, 80, 60, 40, 30, 20, 10]
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    for idx, alpha in enumerate(true_alphas):
        u_clean = solve_fractional_system(alpha, Nt, dt, Nx, y0)
        ax = axes[idx]
        
        snr_selections = []
        snr_errors = []
        
        print(f"\nNoise sweep for true alpha = {alpha}:")
        for snr in snrs:
            # Seed-specific noise for reproducibility
            np.random.seed(snr + 42)
            u_noisy = add_noise(u_clean, snr)
            
            sweep, best = select_temporal_order(u_noisy, dt, dx, alpha_t_candidates=candidates, 
                                                 frac_space_orders=frac_space_orders, threshold=0.01, k_start=20, u_clean=u_clean, true_alpha=alpha)
            
            selected_alpha = best['alpha_t']
            error = abs(selected_alpha - alpha)
            snr_selections.append(selected_alpha)
            snr_errors.append(error)
            
            results[str(alpha)]['noise_sweep'][str(snr)] = {
                'selected_alpha': float(selected_alpha),
                'error': float(error),
                'r2_avg': float(best['r2_avg']),
                'sweep_r2s': {str(res['alpha_t']): float(res['r2_avg']) for res in sweep}
            }
            print(f"  SNR = {snr} dB -> Selected alpha: {selected_alpha} (error: {error:.2f}, R2: {best['r2_avg']:.4f})")
            
        # Plotting the selection vs SNR
        ax.plot(snrs, snr_selections, 'o-', linewidth=2, label='Recovered $\\hat{\\alpha}_t$')
        ax.axhline(alpha, color='r', linestyle='--', label=f'True $\\alpha_t = {alpha}$')
        ax.set_title(f"True $\\alpha_t = {alpha}$ Recovery under Noise")
        ax.set_xlabel("SNR (dB)")
        ax.set_ylabel("Recovered Order $\\hat{\\alpha}_t$")
        ax.set_xscale('log')
        ax.invert_xaxis()  # 100 dB on left, 10 dB on right
        ax.set_xticks(snrs)
        ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
        ax.grid(True)
        ax.legend()
        
    plt.suptitle("Fractional-SINDy Temporal Order Recovery Sensitivity vs Noise", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig('figures/fractional_sindy_sensitivity.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("\nSaved sensitivity plot to figures/fractional_sindy_sensitivity.png")
    
    # Determine the breakdown threshold (minimum SNR where true alpha is successfully recovered)
    print("\n--- Summary & Conclusion ---")
    for alpha in true_alphas:
        breakdown_snr = None
        # Start from high SNR to low SNR
        for snr in snrs:
            err = results[str(alpha)]['noise_sweep'][str(snr)]['error']
            if err < 1e-5:
                breakdown_snr = snr
            else:
                break
        
        if breakdown_snr is not None:
            print(f"True alpha = {alpha}: successfully recovered down to SNR = {breakdown_snr} dB.")
        else:
            print(f"True alpha = {alpha}: failed recovery even at SNR = 100 dB.")
            
    # Save raw json results
    with open('data/identifiability_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("Saved raw results to data/identifiability_results.json")

if __name__ == '__main__':
    main()
