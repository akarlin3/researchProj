import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Ensure figures directory exists
os.makedirs('figures', exist_ok=True)

# 1. Load system parameters from simulator or define them
Nx = 100
L = 10.0
dx = L / (Nx - 1)
params = (Nx, dx, 0.05, 0.1, 0.01, 0.05, 0.1, 0.2, 0.3, 0.2, 0.1)

# Import stroma RHS from simulator
try:
    from ouroboros_sim import stroma_rhs
except ImportError:
    # Fallback definition of stroma_rhs and laplacian if import fails
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

# RK4 integrator step
def rk4_step(f, t, y, dt, *args):
    k1 = f(t, y, *args)
    k2 = f(t + dt/2, y + dt/2 * k1, *args)
    k3 = f(t + dt/2, y + dt/2 * k2, *args)
    k4 = f(t + dt, y + dt * k3, *args)
    return y + dt/6 * (k1 + 2*k2 + 2*k3 + k4)

# 2. Benettin Lyapunov Exponent Estimation
def run_benettin_lyapunov(T_trans=20.0, T_run=80.0, dt=0.005, dt_renorm=0.1, eps=1e-8, seed=42):
    np.random.seed(seed)
    x = np.linspace(0, L, Nx)
    p0 = np.exp(-(x - L/2)**2 / 2.0)
    c0 = 1.0 - 0.5 * np.exp(-(x - L/2)**2 / 2.0)
    n0 = 0.1 + 0.2 * np.sin(np.pi * x / L)
    y_curr = np.concatenate([p0, c0, n0])

    # Transient integration
    t = 0.0
    print("Benettin: Running transient integration...")
    while t < T_trans:
        y_curr = rk4_step(stroma_rhs, t, y_curr, dt, *params)
        t += dt

    # Lyapunov integration setup
    v = np.random.randn(3 * Nx)
    v = v / np.linalg.norm(v)

    lyap_sum = 0.0
    lambda_history = []
    t_history = []

    steps_per_renorm = int(dt_renorm / dt)
    n_renorms = int(T_run / dt_renorm)

    y_and_v = np.concatenate([y_curr, v])

    # Coupled RHS for state + tangent vector (Jacobian-free product)
    def coupled_rhs(t, state):
        y_vec = state[:3*Nx]
        v_vec = state[3*Nx:]
        fy = stroma_rhs(t, y_vec, *params)
        f_perturbed = stroma_rhs(t, y_vec + eps * v_vec, *params)
        jv = (f_perturbed - fy) / eps
        return np.concatenate([fy, jv])

    print("Benettin: Running Lyapunov exponent estimation...")
    for r in range(1, n_renorms + 1):
        for _ in range(steps_per_renorm):
            y_and_v = rk4_step(coupled_rhs, t, y_and_v, dt)
            t += dt

        # Renormalize tangent vector
        y = y_and_v[:3*Nx]
        v_vec = y_and_v[3*Nx:]
        v_norm = np.linalg.norm(v_vec)

        lyap_sum += np.log(v_norm)
        v_vec = v_vec / v_norm
        y_and_v = np.concatenate([y, v_vec])

        current_lambda = lyap_sum / (r * dt_renorm)
        lambda_history.append(current_lambda)
        t_history.append(t - T_trans)

    return np.array(t_history), np.array(lambda_history), y_and_v[:3*Nx]

# 3. Data-Driven LLE Estimator (Rosenstein)
def rosenstein_estimator(time_series, dt, m, tau, theiler_window=50, max_steps=100):
    N = len(time_series)
    M = N - (m - 1) * tau
    if M <= 0:
        raise ValueError("Time series too short for embedding parameters")

    X = np.zeros((M, m))
    for i in range(M):
        X[i] = time_series[i : i + m * tau : tau]

    # Find nearest neighbors excluding Theiler window
    nn_indices = []
    for i in range(M):
        dists = np.linalg.norm(X - X[i], axis=1)
        start = max(0, i - theiler_window)
        end = min(M, i + theiler_window + 1)
        dists[start:end] = np.inf
        nn_indices.append(np.argmin(dists))

    # Track divergence
    d_avg = np.zeros(max_steps)
    for k in range(max_steps):
        s = 0.0
        c = 0
        for i in range(M):
            j = nn_indices[i]
            if i + k < M and j + k < M:
                dist = np.linalg.norm(X[i + k] - X[j + k])
                if dist > 0:
                    s += np.log(dist)
                    c += 1
        d_avg[k] = s / c if c > 0 else np.nan

    return d_avg

# 4. Average Mutual Information (AMI)
def compute_ami(time_series, max_lag=50, n_bins=15):
    ami = []
    n = len(time_series)
    for lag in range(1, max_lag + 1):
        x = time_series[:-lag]
        y = time_series[lag:]
        
        # Compute joint and marginal probabilities
        hist_2d, _, _ = np.histogram2d(x, y, bins=n_bins)
        p_joint = hist_2d / np.sum(hist_2d)
        
        p_x = np.sum(p_joint, axis=1)
        p_y = np.sum(p_joint, axis=0)
        
        mi = 0.0
        for i in range(n_bins):
            for j in range(n_bins):
                if p_joint[i, j] > 1e-12 and p_x[i] > 1e-12 and p_y[j] > 1e-12:
                    mi += p_joint[i, j] * np.log2(p_joint[i, j] / (p_x[i] * p_y[j]))
        ami.append(mi)
    return np.array(ami)

# 5. False Nearest Neighbors (FNN)
def compute_fnn(time_series, tau, max_m=5, R_tol=15.0, A_tol=2.0):
    fnn_percentages = []
    std_x = np.std(time_series)
    
    for m in range(1, max_m + 1):
        N = len(time_series)
        M = N - m * tau
        if M <= 0:
            fnn_percentages.append(0.0)
            continue
            
        # Build X in dimension m
        X = np.zeros((M, m))
        for i in range(M):
            X[i] = time_series[i : i + m * tau : tau]
            
        # Find nearest neighbor distance for each point in dimension m
        false_count = 0
        for i in range(M):
            # Exclude self
            dists = np.linalg.norm(X - X[i], axis=1)
            dists[i] = np.inf
            nn_idx = np.argmin(dists)
            dist_m = dists[nn_idx]
            
            # Check false nearest neighbor criteria in dimension m+1
            val_i_next = time_series[i + m * tau]
            val_nn_next = time_series[nn_idx + m * tau]
            diff_next = abs(val_i_next - val_nn_next)
            
            if dist_m > 0:
                cond1 = diff_next / dist_m > R_tol
                cond2 = np.sqrt(dist_m**2 + diff_next**2) / std_x > A_tol
                if cond1 or cond2:
                    false_count += 1
                    
        fnn_percentages.append((false_count / M) * 100.0)
        
    return np.array(fnn_percentages)

# 6. Correlation Dimension (Grassberger-Procaccia)
def compute_correlation_dimension(time_series, m, tau, r_candidates):
    N = len(time_series)
    M = N - (m - 1) * tau
    X = np.zeros((M, m))
    for i in range(M):
        X[i] = time_series[i : i + m * tau : tau]
        
    # Calculate all pairwise distances (upper triangle)
    dists = []
    for i in range(M):
        for j in range(i + 1, M):
            dists.append(np.linalg.norm(X[i] - X[j]))
    dists = np.array(dists)
    
    C_r = []
    for r in r_candidates:
        count = np.sum(dists < r)
        C_r.append(count / len(dists))
        
    return np.array(C_r)

def main():
    print("="*60)
    print("OUROBOROS CHAOS DIAGNOSTICS")
    print("="*60)
    
    # Run Benettin Exponent integration
    t_history, lambda_history, y_final = run_benettin_lyapunov(
        T_trans=20.0, T_run=80.0, dt=0.005, dt_renorm=0.1
    )
    lambda_max = lambda_history[-1]
    
    print(f"\nBenettin Lyapunov Exponent Estimation Result:")
    print(f"Largest Lyapunov Exponent (LLE) lambda_max = {lambda_max:.6f}")
    
    # Save Lyapunov exponent convergence plot
    plt.figure(figsize=(8, 5))
    plt.plot(t_history, lambda_history, 'b-', linewidth=2, label=r'LLE $\lambda_{\max}$')
    plt.axhline(0.0, color='r', linestyle='--', alpha=0.7)
    plt.title(f"Lyapunov Exponent Convergence (Final $\\lambda_{{\\max}}$ = {lambda_max:.5f})", fontsize=12)
    plt.xlabel("Integration Time (t)", fontsize=10)
    plt.ylabel(r"$\lambda$ estimate", fontsize=10)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig('figures/lyapunov_convergence.png', dpi=150)
    plt.close()
    print("Saved Benettin convergence plot to figures/lyapunov_convergence.png")
    
    # Generate time series for data-driven estimators (T=150 to get a longer time-series)
    print("\nGenerating long-term simulator time series for data-driven cross-check...")
    x = np.linspace(0, L, Nx)
    p0 = np.exp(-(x - L/2)**2 / 2.0)
    c0 = 1.0 - 0.5 * np.exp(-(x - L/2)**2 / 2.0)
    n0 = 0.1 + 0.2 * np.sin(np.pi * x / L)
    y0 = np.concatenate([p0, c0, n0])
    
    T_sim = 150.0
    Nt_sim = 3000
    sol = solve_ivp(
        stroma_rhs,
        [0, T_sim],
        y0,
        t_eval=np.linspace(0, T_sim, Nt_sim),
        method='Radau',
        args=params
    )
    
    # Representative scalar observable: spatial mean of pressure p
    ts = np.mean(sol.y[:Nx, :], axis=0)
    dt_sim = T_sim / (Nt_sim - 1)
    
    # Run Mutual Information to choose lag tau
    print("Computing Average Mutual Information (AMI)...")
    ami_vals = compute_ami(ts, max_lag=100, n_bins=15)
    
    # Choose tau: first local minimum or drop below 1/e
    tau = 1
    for i in range(1, len(ami_vals) - 1):
        if ami_vals[i] < ami_vals[i-1] and ami_vals[i] < ami_vals[i+1]:
            tau = i + 1
            break
    if tau == 1: # if no local minimum, use drop below 1/e
        for i in range(len(ami_vals)):
            if ami_vals[i] < ami_vals[0] / np.e:
                tau = i + 1
                break
    print(f"Chosen delay lag (tau) = {tau} steps ({tau * dt_sim:.3f} time units)")
    
    # Save AMI plot
    plt.figure(figsize=(8, 4))
    plt.plot(np.arange(1, 101) * dt_sim, ami_vals, 'k-', linewidth=2)
    plt.axvline(tau * dt_sim, color='r', linestyle='--', label=f'Chosen lag = {tau * dt_sim:.2f} s')
    plt.title("Average Mutual Information vs Delay Lag", fontsize=12)
    plt.xlabel("Delay (time units)", fontsize=10)
    plt.ylabel("Mutual Information (bits)", fontsize=10)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig('figures/ami_lag_selection.png', dpi=150)
    plt.close()
    print("Saved AMI plot to figures/ami_lag_selection.png")
    
    # Run False Nearest Neighbors to choose embedding dimension m
    print("Computing False Nearest Neighbors (FNN)...")
    fnn_vals = compute_fnn(ts, tau, max_m=5)
    chosen_m = 1
    for m_idx, fnn_val in enumerate(fnn_vals):
        if fnn_val < 1.0: # threshold of 1%
            chosen_m = m_idx + 1
            break
    print(f"Chosen embedding dimension (m) = {chosen_m}")
    
    # Save FNN plot
    plt.figure(figsize=(8, 4))
    plt.plot(np.arange(1, 6), fnn_vals, 'mo-', linewidth=2, markersize=8)
    plt.axhline(1.0, color='r', linestyle='--', label='1% Threshold')
    plt.title("Percentage of False Nearest Neighbors vs Embedding Dimension", fontsize=12)
    plt.xlabel("Embedding Dimension m", fontsize=10)
    plt.ylabel("FNN %", fontsize=10)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig('figures/fnn_dimension_selection.png', dpi=150)
    plt.close()
    print("Saved FNN plot to figures/fnn_dimension_selection.png")
    
    # Run Rosenstein estimator
    print("Running Rosenstein data-driven Lyapunov exponent estimator...")
    # Use standard m=3, tau as chosen, or default parameters to be robust
    max_steps = 150
    rosen_curve = rosenstein_estimator(ts, dt_sim, m=chosen_m, tau=tau, theiler_window=50, max_steps=max_steps)
    
    t_rosen = np.arange(max_steps) * dt_sim
    
    plt.figure(figsize=(8, 5))
    plt.plot(t_rosen, rosen_curve, 'g-', linewidth=2, label="Rosenstein curve")
    plt.title("Rosenstein Divergence Curve (Data-Driven)", fontsize=12)
    plt.xlabel("Time separation (t)", fontsize=10)
    plt.ylabel(r"$\langle \ln d(t) \rangle$", fontsize=10)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig('figures/rosenstein_divergence.png', dpi=150)
    plt.close()
    print("Saved Rosenstein plot to figures/rosenstein_divergence.png")
    
    # Perform linear regression to estimate slope (LLE)
    # Fit over early/mid time steps to evaluate
    slope, intercept = np.polyfit(t_rosen, rosen_curve, 1)
    print(f"Data-driven Rosenstein LLE slope estimate: {slope:.6f}")
    
    # Attractor Reconstruction (if chaotic)
    if lambda_max > 0:
        print("\nSystem is chaotic (lambda_max > 0). Reconstructing attractor and correlation dimension...")
        
        # 3D Attractor plot
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection='3d')
        
        # Reconstruct phase space
        M_embed = len(ts) - 2 * tau
        x_coord = ts[:-2*tau]
        y_coord = ts[tau:-tau]
        z_coord = ts[2*tau:]
        
        ax.plot(x_coord, y_coord, z_coord, 'c-', alpha=0.7)
        ax.set_title("Reconstructed Attractor (Delay Embedding, 3D)")
        ax.set_xlabel("x(t)")
        ax.set_ylabel("x(t + tau)")
        ax.set_zlabel("x(t + 2*tau)")
        plt.savefig('figures/attractor_embedding.png', dpi=150)
        plt.close()
        print("Saved reconstructed attractor plot to figures/attractor_embedding.png")
        
        # Correlation Dimension via Grassberger-Procaccia
        r_vals = np.logspace(-4, -1, 30)
        C_r = compute_correlation_dimension(ts, chosen_m, tau, r_vals)
        
        # Fit slope in the scaling region (avoiding noise floor and saturation)
        valid = (C_r > 1e-4) & (C_r < 0.5)
        if np.sum(valid) > 3:
            gp_slope, _ = np.polyfit(np.log(r_vals[valid]), np.log(C_r[valid]), 1)
        else:
            gp_slope = 0.0
            
        print(f"Estimated Correlation Dimension (Grassberger-Procaccia): {gp_slope:.4f}")
        
        plt.figure(figsize=(8, 5))
        plt.loglog(r_vals, C_r, 'o-')
        plt.title(f"Correlation Sum C(r) vs scale r (D2 = {gp_slope:.3f})")
        plt.xlabel("Scale r")
        plt.ylabel("C(r)")
        plt.grid(True, which="both", ls="-")
        plt.tight_layout()
        plt.savefig('figures/correlation_dimension.png', dpi=150)
        plt.close()
        print("Saved correlation dimension plot to figures/correlation_dimension.png")
    else:
        print("\nSystem is NON-CHAOTIC (lambda_max <= 0). Skipping attractor reconstruction and Grassberger-Procaccia dimension.")
        
    print("\nChaos diagnostics complete.")
    print("="*60)

if __name__ == '__main__':
    main()
