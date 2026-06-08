import os
import json
import numpy as np
import matplotlib.pyplot as plt
from ouroboros_fractional_sindy import gl_weights, FractionalPDELibrary, build_fractional_pde_model
from ouroboros_identifiability import solve_fractional_system, add_noise
from ouroboros_model_select import select_temporal_order
from ouroboros_mitigation import weak_gl_derivative_time, select_temporal_order_weak
import pysindy as ps
from scipy.linalg import toeplitz

# Suppress all numpy warnings for clean log and speed
np.seterr(all='ignore')
np.random.seed(42)

# Grid parameters
Nx = 50
Nt = 500
T = 5.0
L = 10.0
dt = T / (Nt - 1)
dx = L / (Nx - 1)
k_start = 20

# Fast STLSQ implementation
def fast_stlsq(X, Y, threshold, alpha=1e-4, max_iter=10):
    n_features = X.shape[1]
    n_targets = Y.shape[1]
    
    norms = np.linalg.norm(X, axis=0)
    norms[norms < 1e-10] = 1.0
    X_normalized = X / norms
    
    coefs = np.zeros((n_targets, n_features))
    
    for d in range(n_targets):
        y = Y[:, d]
        XtX = X_normalized.T @ X_normalized
        Xty = X_normalized.T @ y
        
        beta = np.linalg.solve(XtX + alpha * np.eye(n_features), Xty)
        
        for _ in range(max_iter):
            small_inds = np.abs(beta) < threshold
            beta[small_inds] = 0.0
            big_inds = ~small_inds
            if not np.any(big_inds):
                break
            X_sub = XtX[np.ix_(big_inds, big_inds)]
            Xty_sub = Xty[big_inds]
            beta_sub = np.linalg.solve(X_sub + alpha * np.eye(np.sum(big_inds)), Xty_sub)
            beta_new = np.zeros(n_features)
            beta_new[big_inds] = beta_sub
            if np.allclose(beta, beta_new):
                break
            beta = beta_new
            
        coefs[d] = beta / norms
        
    return coefs

def get_features(u, dx, frac_space_orders):
    library = FractionalPDELibrary(dx=dx, beta_candidates=frac_space_orders, include_interaction=False)
    library.fit(u)
    Theta = library.transform(u)
    return np.asarray(Theta), library

# Optimized temporal Grünwald-Letnikov derivative using Toeplitz matrix multiplication
def gl_derivative_time(u, dt, order):
    u = np.asarray(u)
    nt, nx, n_vars = u.shape
    w = gl_weights(order, nt)
    r = np.zeros(nt)
    r[0] = w[0]
    W = toeplitz(w, r)
    u_flat = u.reshape(nt, -1)
    u_diff = (W @ u_flat).reshape(nt, nx, n_vars)
    return u_diff / (dt**order)

# Optimized spatial Grünwald-Letnikov derivative using Toeplitz matrix multiplication
def gl_derivative_space(u, dx, order):
    u = np.asarray(u)
    if u.ndim == 3:
        nt, nx, n_vars = u.shape
        w = gl_weights(order, nx)
        r = np.zeros(nx)
        r[0] = w[0]
        W_space = toeplitz(w, r)
        u_flat = u.transpose(0, 2, 1).reshape(-1, nx)
        u_diff_flat = u_flat @ W_space.T
        u_diff = u_diff_flat.reshape(nt, n_vars, nx).transpose(0, 2, 1)
        return u_diff / (dx**order)
    else:
        nx, n_vars = u.shape
        w = gl_weights(order, nx)
        r = np.zeros(nx)
        r[0] = w[0]
        W_space = toeplitz(w, r)
        u_diff = (u @ W_space.T)
        return u_diff / (dx**order)

# Monkey-patch the spatial derivative into the imported module
import ouroboros_fractional_sindy
ouroboros_fractional_sindy.gl_derivative_space = gl_derivative_space

def select_temporal_order_pointwise_fast(u, u_clean, dt, dx, candidates, frac_space_orders, threshold, k_start, true_alpha, Theta_test_flat, u_dot_test_clean_flat):
    nt, nx, n_vars = u.shape
    nt_train = int(0.8 * nt)
    
    u_train = u[:nt_train]
    Theta_train, _ = get_features(u_train, dx, frac_space_orders)
    Theta_train_flat = Theta_train[k_start:].reshape(-1, Theta_train.shape[-1])
    
    sweep_results = []
    
    for alpha_t in candidates:
        u_dot = gl_derivative_time(u, dt, alpha_t)
        u_dot_train_slice = u_dot[k_start:nt_train].reshape(-1, n_vars)
        
        coefs = fast_stlsq(Theta_train_flat, u_dot_train_slice, threshold=threshold)
        
        u_dot_pred_flat = Theta_test_flat @ coefs.T
        
        r2_scores = []
        for i in range(n_vars):
            y_true = u_dot_test_clean_flat[:, i]
            y_pred = u_dot_pred_flat[:, i]
            
            ss_res = np.sum((y_true - y_pred)**2)
            ss_tot = np.sum((y_true - np.mean(y_true))**2)
            r2 = 1.0 - (ss_res / (ss_tot + 1e-10))
            r2_scores.append(r2)
            
        r2_avg = np.mean(r2_scores)
        sweep_results.append({
            'alpha_t': alpha_t,
            'r2_avg': r2_avg
        })
        
    best_idx = np.argmax([res['r2_avg'] for res in sweep_results])
    return sweep_results, sweep_results[best_idx]

# Helper functions for vectorized projections
def get_psi_matrix(nt, train_windows, alpha, w, phi):
    n_windows = len(train_windows)
    Psi = np.zeros((nt, n_windows))
    for j, (k_a, k_b) in enumerate(train_windows):
        psi = np.zeros(k_b + 1)
        for k in range(k_a, k_b + 1):
            psi[:k+1] += w[k::-1] * phi[k - k_a]
        Psi[:k_b+1, j] = psi
    return Psi

def get_phi_matrix(nt, train_windows, phi):
    n_windows = len(train_windows)
    Phi_mat = np.zeros((nt, n_windows))
    for j, (k_a, k_b) in enumerate(train_windows):
        Phi_mat[k_a:k_b+1, j] = phi
    return Phi_mat

def select_temporal_order_weak_fast(u, u_clean, dt, dx, candidates, frac_space_orders, threshold, k_start, true_alpha, Phi_mat_train, X_test, Y_test):
    nt, nx, n_vars = u.shape
    nt_train = int(0.8 * nt)
    
    H = 50
    S = 5
    train_windows = []
    
    for k_a in range(0, nt - H + 1, S):
        k_b = k_a + H
        if k_b < nt_train:
            if k_a >= k_start:
                train_windows.append((k_a, k_b))
                
    # Train features projection
    Theta_train, _ = get_features(u[:nt_train], dx, frac_space_orders)
    X_train_temp = np.tensordot(Theta_train, Phi_mat_train, axes=([0], [0])) * dt # (nx, n_features, n_windows)
    X_train = np.transpose(X_train_temp, (2, 0, 1)).reshape(-1, Theta_train.shape[-1])
    
    sweep_results = []
    for alpha_t in candidates:
        w_train = gl_weights(alpha_t, nt)
        Psi_train = get_psi_matrix(nt_train, train_windows, alpha_t, w_train, phi=Phi_mat_train[train_windows[0][0]:train_windows[0][1]+1, 0])
        
        # Project noisy u_train
        Y_train_temp = np.tensordot(u[:nt_train], Psi_train, axes=([0], [0])) # (nx, n_vars, n_windows)
        Y_train = np.transpose(Y_train_temp, (2, 0, 1)).reshape(-1, n_vars) * (dt / (dt**alpha_t))
        
        coefs = fast_stlsq(X_train, Y_train, threshold=threshold)
        
        Y_pred = X_test @ coefs.T
        
        r2_scores = []
        for i in range(n_vars):
            y_true = Y_test[:, i]
            y_pred_i = Y_pred[:, i]
            
            ss_res = np.sum((y_true - y_pred_i)**2)
            ss_tot = np.sum((y_true - np.mean(y_true))**2)
            r2 = 1.0 - (ss_res / (ss_tot + 1e-10))
            r2_scores.append(r2)
            
        r2_avg = np.mean(r2_scores)
        sweep_results.append({
            'alpha_t': alpha_t,
            'r2_avg': r2_avg
        })
        
    best_idx = np.argmax([res['r2_avg'] for res in sweep_results])
    return sweep_results, sweep_results[best_idx]

def run_sanity_checks():
    print("Running sanity checks to verify fast implementation equivalence...")
    alpha = 0.7
    x = np.linspace(0, L, Nx)
    p0 = np.exp(-(x - L/2)**2 / 2.0)
    c0 = 1.0 - 0.5 * np.exp(-(x - L/2)**2 / 2.0)
    n0 = 0.1 + 0.2 * np.sin(np.pi * x / L)
    y0 = np.concatenate([p0, c0, n0])
    u_clean = solve_fractional_system(alpha, Nt, dt, Nx, y0)
    u_noisy = add_noise(u_clean, 40)
    candidates = [0.5, 0.6, 0.7, 0.8, 0.9]
    frac_space_orders = (0.5, 1.5)
    
    nt_train = int(0.8 * Nt)
    u_clean_test = u_clean[nt_train:]
    Theta_test, _ = get_features(u_clean_test, dx, frac_space_orders)
    Theta_test_flat = Theta_test.reshape(-1, Theta_test.shape[-1])
    u_dot_clean = gl_derivative_time(u_clean, dt, alpha)
    u_dot_test_clean = u_dot_clean[nt_train:]
    u_dot_test_clean_flat = u_dot_test_clean.reshape(-1, 3)
    
    # Check Pointwise
    _, best_slow_p = select_temporal_order(u_noisy, dt, dx, alpha_t_candidates=candidates, 
                                           frac_space_orders=frac_space_orders, threshold=0.01, k_start=k_start, u_clean=u_clean, true_alpha=alpha)
    _, best_fast_p = select_temporal_order_pointwise_fast(
        u_noisy, u_clean, dt, dx, candidates, frac_space_orders, 0.01, k_start, alpha, Theta_test_flat, u_dot_test_clean_flat
    )
    
    print(f"Pointwise: Slow best = {best_slow_p['alpha_t']:.2f} (R2={best_slow_p['r2_avg']:.4f}), Fast best = {best_fast_p['alpha_t']:.2f} (R2={best_fast_p['r2_avg']:.4f})")
    assert abs(best_slow_p['alpha_t'] - best_fast_p['alpha_t']) < 1e-5, "Pointwise alpha mismatch!"
    
    # Check Weak-form
    H = 50
    S = 5
    train_windows = []
    test_windows = []
    for k_a in range(0, Nt - H + 1, S):
        k_b = k_a + H
        if k_b < nt_train:
            if k_a >= k_start:
                train_windows.append((k_a, k_b))
        elif k_a >= nt_train:
            if k_b < Nt:
                test_windows.append((k_a, k_b))
    phi = (1.0 - ((2.0 * np.arange(H + 1) / H) - 1.0)**2) ** 4
    phi /= np.sum(phi) * dt
    Phi_mat_train = get_phi_matrix(nt_train, train_windows, phi)
    Phi_mat_test = get_phi_matrix(Nt - nt_train, [(k_a - nt_train, k_b - nt_train) for k_a, k_b in test_windows], phi)
    
    Theta_clean, _ = get_features(u_clean, dx, frac_space_orders)
    Theta_clean_test = Theta_clean[nt_train:]
    X_test_temp = np.tensordot(Theta_clean_test, Phi_mat_test, axes=([0], [0])) * dt
    X_test_weak = np.transpose(X_test_temp, (2, 0, 1)).reshape(-1, Theta_clean.shape[-1])
    
    w_test = gl_weights(alpha, Nt)
    Psi_test = get_psi_matrix(Nt, test_windows, alpha, w_test, phi)
    Y_test_temp = np.tensordot(u_clean, Psi_test, axes=([0], [0]))
    Y_test_weak = np.transpose(Y_test_temp, (2, 0, 1)).reshape(-1, 3) * (dt / (dt**alpha))
    
    _, best_slow_w = select_temporal_order_weak(u_noisy, u_clean, dt, dx, candidates=candidates,
                                                frac_space_orders=frac_space_orders, threshold=0.01, k_start=k_start, true_alpha=alpha)
    _, best_fast_w = select_temporal_order_weak_fast(
        u_noisy, u_clean, dt, dx, candidates, frac_space_orders, 0.01, k_start, alpha, Phi_mat_train, X_test_weak, Y_test_weak
    )
    
    print(f"Weak-form: Slow best = {best_slow_w['alpha_t']:.2f} (R2={best_slow_w['r2_avg']:.4f}), Fast best = {best_fast_w['alpha_t']:.2f} (R2={best_fast_w['r2_avg']:.4f})")
    assert abs(best_slow_w['alpha_t'] - best_fast_w['alpha_t']) < 1e-5, "Weak-form alpha mismatch!"
    print("Sanity checks PASSED successfully!\n")

def main():
    run_sanity_checks()
    
    true_alphas = [0.5, 0.7, 0.9]
    candidates = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    frac_space_orders = (0.5, 1.5)
    threshold = 0.01
    
    # Define fine SNR sweeps based on coarse breakpoints
    snr_sweeps = {
        0.5: {
            'pointwise': np.arange(20, 61, 5), # 20, 25, 30, 35, 40, 45, 50, 55, 60
            'weak': np.arange(20, 61, 5)
        },
        0.7: {
            'pointwise': np.arange(20, 61, 5),
            'weak': np.arange(0, 41, 5) # 0, 5, 10, 15, 20, 25, 30, 35, 40
        },
        0.9: {
            'pointwise': np.arange(40, 81, 5), # 40, 45, 50, 55, 60, 65, 70, 75, 80
            'weak': np.arange(0, 31, 5) # 0, 5, 10, 15, 20, 25, 30
        }
    }
    
    # Number of realizations
    N_trials = 500
    
    # We will compute results and save them
    sweep_results = {}
    
    # Spatial initial conditions for solving the systems
    x = np.linspace(0, L, Nx)
    p0 = np.exp(-(x - L/2)**2 / 2.0)
    c0 = 1.0 - 0.5 * np.exp(-(x - L/2)**2 / 2.0)
    n0 = 0.1 + 0.2 * np.sin(np.pi * x / L)
    y0 = np.concatenate([p0, c0, n0])
    
    for alpha in true_alphas:
        print(f"==================================================")
        print(f"RUNNING FINE SNR SWEEP FOR TRUE ALPHA = {alpha}")
        print(f"==================================================")
        
        # 1. Solve the system
        u_clean = solve_fractional_system(alpha, Nt, dt, Nx, y0)
        
        sweep_results[str(alpha)] = {
            'pointwise': {},
            'weak': {}
        }
        
        # Pointwise precomputations
        nt_train = int(0.8 * Nt)
        u_clean_test = u_clean[nt_train:]
        Theta_test, _ = get_features(u_clean_test, dx, frac_space_orders)
        Theta_test_flat = Theta_test.reshape(-1, Theta_test.shape[-1])
        u_dot_clean = gl_derivative_time(u_clean, dt, alpha)
        u_dot_test_clean = u_dot_clean[nt_train:]
        u_dot_test_clean_flat = u_dot_test_clean.reshape(-1, 3)
        
        # 2. Pointwise Sweep
        pw_snrs = snr_sweeps[alpha]['pointwise']
        print(f"\nPointwise Sweep (SNRs: {pw_snrs}):")
        for snr in pw_snrs:
            selections = []
            np.random.seed(int(snr) + 42)
            
            for trial in range(N_trials):
                u_noisy = add_noise(u_clean, snr)
                best_alpha = select_temporal_order_pointwise_fast(
                    u_noisy, u_clean, dt, dx, candidates, frac_space_orders, threshold, k_start, alpha, Theta_test_flat, u_dot_test_clean_flat
                )[1]['alpha_t']
                selections.append(best_alpha)
                
            selections = np.array(selections)
            success_rate = np.mean(np.abs(selections - alpha) < 1e-5)
            
            errors = np.abs(selections - alpha)
            mean_error = np.mean(errors)
            std_error = np.std(errors)
            
            sweep_results[str(alpha)]['pointwise'][str(snr)] = {
                'success_rate': float(success_rate),
                'mean_error': float(mean_error),
                'std_error': float(std_error),
                'selections': selections.tolist()
            }
            print(f"  SNR = {snr} dB -> Success Rate: {success_rate:.3f} | Error: {mean_error:.4f} ± {std_error:.4f}")
            
        # Weak-form precomputations
        H = 50
        S = 5
        train_windows = []
        test_windows = []
        for k_a in range(0, Nt - H + 1, S):
            k_b = k_a + H
            if k_b < nt_train:
                if k_a >= k_start:
                    train_windows.append((k_a, k_b))
            elif k_a >= nt_train:
                if k_b < Nt:
                    test_windows.append((k_a, k_b))
                    
        phi = (1.0 - ((2.0 * np.arange(H + 1) / H) - 1.0)**2) ** 4
        phi /= np.sum(phi) * dt
        Phi_mat_train = get_phi_matrix(nt_train, train_windows, phi)
        Phi_mat_test = get_phi_matrix(Nt - nt_train, [(k_a - nt_train, k_b - nt_train) for k_a, k_b in test_windows], phi)
        
        Theta_clean, _ = get_features(u_clean, dx, frac_space_orders)
        Theta_clean_test = Theta_clean[nt_train:]
        X_test_temp = np.tensordot(Theta_clean_test, Phi_mat_test, axes=([0], [0])) * dt
        X_test_weak = np.transpose(X_test_temp, (2, 0, 1)).reshape(-1, Theta_clean.shape[-1])
        
        w_test = gl_weights(alpha, Nt)
        Psi_test = get_psi_matrix(Nt, test_windows, alpha, w_test, phi)
        Y_test_temp = np.tensordot(u_clean, Psi_test, axes=([0], [0]))
        Y_test_weak = np.transpose(Y_test_temp, (2, 0, 1)).reshape(-1, 3) * (dt / (dt**alpha))
        
        # 3. Weak-form Sweep
        wk_snrs = snr_sweeps[alpha]['weak']
        print(f"\nWeak-form Sweep (SNRs: {wk_snrs}):")
        for snr in wk_snrs:
            selections = []
            np.random.seed(int(snr) + 42)
            
            for trial in range(N_trials):
                u_noisy = add_noise(u_clean, snr)
                best_alpha = select_temporal_order_weak_fast(
                    u_noisy, u_clean, dt, dx, candidates, frac_space_orders, threshold, k_start, alpha, Phi_mat_train, X_test_weak, Y_test_weak
                )[1]['alpha_t']
                selections.append(best_alpha)
                
            selections = np.array(selections)
            success_rate = np.mean(np.abs(selections - alpha) < 1e-5)
            
            errors = np.abs(selections - alpha)
            mean_error = np.mean(errors)
            std_error = np.std(errors)
            
            sweep_results[str(alpha)]['weak'][str(snr)] = {
                'success_rate': float(success_rate),
                'mean_error': float(mean_error),
                'std_error': float(std_error),
                'selections': selections.tolist()
            }
            print(f"  SNR = {snr} dB -> Success Rate: {success_rate:.3f} | Error: {mean_error:.4f} ± {std_error:.4f}")
            
    # Save results to JSON
    with open('data/fine_snr_sweep_results.json', 'w') as f:
        json.dump(sweep_results, f, indent=2)
    print("\nSaved fine SNR sweep results to data/fine_snr_sweep_results.json")
    
    # Calculate breakdown brackets
    print("\n" + "="*50)
    print("RECOVERY BRACKETS PER ORDER")
    print("="*50)
    
    brackets = {}
    
    for alpha in true_alphas:
        brackets[str(alpha)] = {}
        for method in ['pointwise', 'weak']:
            snrs_tested = sorted([float(k) for k in sweep_results[str(alpha)][method].keys()])
            
            succ_snrs = []
            fail_snrs = []
            
            for snr in snrs_tested:
                rate = sweep_results[str(alpha)][method][str(int(snr))]['success_rate']
                if rate >= 0.95:
                    succ_snrs.append(snr)
                else:
                    fail_snrs.append(snr)
                    
            lowest_succ = min(succ_snrs) if len(succ_snrs) > 0 else None
            highest_fail = max(fail_snrs) if len(fail_snrs) > 0 else None
            
            brackets[str(alpha)][method] = {
                'lowest_succ': lowest_succ,
                'highest_fail': highest_fail
            }
            
            print(f"Alpha = {alpha} | Method: {method:<10} | Bracket: [{highest_fail} dB, {lowest_succ} dB]")
            
    # Generate the updated breakdown comparison plot
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    colors = {'pointwise': '#d62728', 'weak': '#1f77b4'}
    markers = {'pointwise': 'o', 'weak': 's'}
    labels = {'pointwise': 'Naive Pointwise GL', 'weak': 'Weak-Form GL'}
    
    for idx, alpha in enumerate(true_alphas):
        ax = axes[idx]
        for method in ['pointwise', 'weak']:
            snrs = sorted([float(k) for k in sweep_results[str(alpha)][method].keys()])
            rates = [sweep_results[str(alpha)][method][str(int(s))]['success_rate'] for s in snrs]
            ax.plot(snrs, rates, marker=markers[method], color=colors[method], label=labels[method], linewidth=2)
            
        ax.set_title(f"True Order $\\alpha_t = {alpha}$ Recovery")
        ax.set_xlabel("SNR (dB)")
        ax.set_ylabel("Success Rate (Exact Recovery)")
        ax.set_ylim(-0.05, 1.05)
        ax.grid(True, linestyle=':')
        ax.legend()
        
    plt.suptitle("Fractional Order Success Rate vs. SNR (5 dB resolution)", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig('figures/mitigation_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved updated breakdown figure to figures/mitigation_comparison.png")
    
    # Write RESULTS_snr_brackets.md
    with open('RESULTS_snr_brackets.md', 'w') as f:
        f.write("# Finer SNR Brackets Report (Checkpoint 2)\n\n")
        f.write("This report presents the refined order-recovery SNR breakdown brackets for temporal fractional order recovery, evaluated on a fine 5 dB grid using 500 noise realizations per cell.\n\n")
        f.write("## 1. Recovery Success Criteria\n")
        f.write("A recovery trial is declared **successful** if the estimated order $\\hat{\\alpha}_t$ matches the true order $\\alpha_t$ exactly ($|\\hat{\\alpha}_t - \\alpha_t| < 1e-5$). We define the **recovery threshold** as a bracket `[highest SNR that fails, lowest SNR that succeeds]`, where success requires a **$\\ge 95\\%$ success rate** across 500 independent noise realizations.\n\n")
        
        f.write("## 2. Table of Recovery Brackets (5 dB steps)\n\n")
        f.write("| True $\\alpha_t$ | Pointwise GL Bracket [Fail, Succeed] | Weak-Form GL Bracket [Fail, Succeed] | Order-Dependent Improvement |\n")
        f.write("| :---: | :---: | :---: | :--- |\n")
        for alpha in true_alphas:
            pw_f = brackets[str(alpha)]['pointwise']['highest_fail']
            pw_s = brackets[str(alpha)]['pointwise']['lowest_succ']
            wk_f = brackets[str(alpha)]['weak']['highest_fail']
            wk_s = brackets[str(alpha)]['weak']['lowest_succ']
            
            pw_str = f"[{pw_f:.0f} dB, {pw_s:.0f} dB]" if pw_s is not None and pw_f is not None else "N/A"
            wk_str = f"[{wk_f:.0f} dB, {wk_s:.0f} dB]" if wk_s is not None and wk_f is not None else "N/A"
            
            if pw_s is not None and wk_s is not None:
                imp = f"Weak-form improves recovery by {pw_s - wk_s:.0f} dB"
            elif pw_s is None and wk_s is not None:
                imp = f"Weak-form enables recovery down to {wk_s:.0f} dB (Pointwise failed completely)"
            else:
                imp = "No successful recovery at studied range"
                
            f.write(f"| {alpha:.1f} | {pw_str} | {wk_str} | {imp} |\n")
            
        f.write("\n## 3. Full Sweep Results Data\n\n")
        for alpha in true_alphas:
            f.write(f"### True $\\alpha_t = {alpha:.1f}$\n\n")
            f.write("#### Pointwise GL:\n")
            for snr in sorted([int(k) for k in sweep_results[str(alpha)]['pointwise'].keys()]):
                rate = sweep_results[str(alpha)]['pointwise'][str(snr)]['success_rate']
                err = sweep_results[str(alpha)]['pointwise'][str(snr)]['mean_error']
                std = sweep_results[str(alpha)]['pointwise'][str(snr)]['std_error']
                f.write(f"- SNR = {snr} dB: Success Rate = {rate:.3f}, Error = {err:.4f} ± {std:.4f}\n")
            f.write("\n#### Weak-Form GL:\n")
            for snr in sorted([int(k) for k in sweep_results[str(alpha)]['weak'].keys()]):
                rate = sweep_results[str(alpha)]['weak'][str(snr)]['success_rate']
                err = sweep_results[str(alpha)]['weak'][str(snr)]['mean_error']
                std = sweep_results[str(alpha)]['weak'][str(snr)]['std_error']
                f.write(f"- SNR = {snr} dB: Success Rate = {rate:.3f}, Error = {err:.4f} ± {std:.4f}\n")
            f.write("\n")
            
        f.write("* **Plot Citation**: ![Mitigation Comparison Plot](file:///Users/averykarlin/projOuroboros/figures/mitigation_comparison.png)\n")
        
    print("Saved RESULTS_snr_brackets.md")

if __name__ == '__main__':
    main()
