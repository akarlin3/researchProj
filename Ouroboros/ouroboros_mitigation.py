import os
import json
import numpy as np
import matplotlib.pyplot as plt
from ouroboros_fractional_sindy import gl_weights, gl_derivative_time, FractionalPDELibrary, build_fractional_pde_model
from ouroboros_identifiability import solve_fractional_system, add_noise
from ouroboros_model_select import select_temporal_order
import pysindy as ps

# Set random seed for reproducibility
np.random.seed(42)

# Grid parameters (same as identifiability)
Nx = 50
Nt = 500
T = 5.0
L = 10.0
dt = T / (Nt - 1)
dx = L / (Nx - 1)

def tikhonov_smooth(u, lambd):
    """
    Applies Tikhonov regularization to smooth u along time (axis 0).
    u shape: (nt, nx, n_vars)
    """
    if lambd <= 0:
        return u.copy()
        
    nt, nx, n_vars = u.shape
    D2 = np.zeros((nt - 2, nt))
    for i in range(nt - 2):
        D2[i, i] = 1.0
        D2[i, i+1] = -2.0
        D2[i, i+2] = 1.0
        
    A = np.eye(nt) + lambd * (D2.T @ D2)
    
    u_smoothed = np.zeros_like(u)
    for j in range(nx):
        for k in range(n_vars):
            u_smoothed[:, j, k] = np.linalg.solve(A, u[:, j, k])
    return u_smoothed

def weak_gl_derivative_time(u, dt, alpha, k_a, k_b, phi, w):
    """
    Computes the weak Grünwald-Letnikov derivative for window [k_a, k_b].
    u shape: (nt, nx, n_vars)
    phi: test function of shape (k_b - k_a + 1,)
    w: GL weights of shape (nt,)
    """
    nt, nx, n_vars = u.shape
    # Calculate psi
    psi = np.zeros(k_b + 1)
    for k in range(k_a, k_b + 1):
        p_val = phi[k - k_a]
        for i in range(k + 1):
            psi[i] += w[k - i] * p_val
    # Compute weak integral: sum_{i=0}^{k_b} u[i] * psi[i]
    y_weak = np.sum(u[:k_b+1] * psi[:, np.newaxis, np.newaxis], axis=0) * (dt / (dt**alpha))
    return y_weak

def get_weak_features_and_target(u, dt, dx, alpha, library, train_windows, test_windows, u_clean=None, true_alpha=None):
    """
    Computes weak features and targets for SINDy on train and test windows.
    """
    Theta = library.transform(u)  # (Nt, Nx, N_features)
    nt, nx, n_vars = u.shape
    n_features = Theta.shape[-1]
    
    H = 50
    p = 4
    phi = (1.0 - ((2.0 * np.arange(H + 1) / H) - 1.0)**2) ** p
    phi /= np.sum(phi) * dt
    
    w = gl_weights(alpha, nt)
    
    # Train
    X_train_list = []
    Y_train_list = []
    for k_a, k_b in train_windows:
        theta_window = Theta[k_a:k_b+1]
        theta_weak = np.sum(theta_window * phi[:, np.newaxis, np.newaxis], axis=0) * dt
        X_train_list.append(theta_weak)
        
        y_weak = weak_gl_derivative_time(u, dt, alpha, k_a, k_b, phi, w)
        Y_train_list.append(y_weak)
        
    X_train = np.concatenate(X_train_list, axis=0)
    Y_train = np.concatenate(Y_train_list, axis=0)
    
    # Test
    u_test_source = u_clean if u_clean is not None else u
    Theta_test = library.transform(u_test_source)
    
    target_alpha = true_alpha if true_alpha is not None else alpha
    w_test = gl_weights(target_alpha, nt)
    
    X_test_list = []
    Y_test_list = []
    for k_a, k_b in test_windows:
        theta_window = Theta_test[k_a:k_b+1]
        theta_weak = np.sum(theta_window * phi[:, np.newaxis, np.newaxis], axis=0) * dt
        X_test_list.append(theta_weak)
        
        y_weak = weak_gl_derivative_time(u_test_source, dt, target_alpha, k_a, k_b, phi, w_test)
        Y_test_list.append(y_weak)
        
    X_test = np.concatenate(X_test_list, axis=0)
    Y_test = np.concatenate(Y_test_list, axis=0)
    
    return X_train, Y_train, X_test, Y_test

def select_temporal_order_weak(u, u_clean, dt, dx, candidates, frac_space_orders, threshold, k_start, true_alpha=None):
    """
    Performs weak-form model selection.
    """
    nt, nx, n_vars = u.shape
    nt_train = int(0.8 * nt)
    
    H = 50
    S = 5
    train_windows = []
    test_windows = []
    
    for k_a in range(0, nt - H + 1, S):
        k_b = k_a + H
        if k_b < nt_train:
            if k_a >= k_start:
                train_windows.append((k_a, k_b))
        elif k_a >= nt_train:
            if k_b < nt:
                test_windows.append((k_a, k_b))
                
    library = FractionalPDELibrary(dx=dx, beta_candidates=frac_space_orders, include_interaction=False)
    library.fit(u)
    
    sweep_results = []
    for alpha_t in candidates:
        X_train, Y_train, X_test, Y_test = get_weak_features_and_target(
            u, dt, dx, alpha_t, library, train_windows, test_windows, u_clean=u_clean, true_alpha=true_alpha
        )
        
        optimizer = ps.STLSQ(threshold=threshold, alpha=1e-4, normalize_columns=True)
        optimizer.fit(X_train, Y_train)
        coefs = optimizer.coef_
        
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
            'r2_avg': r2_avg,
            'r2_scores': r2_scores
        })
        
    best_idx = np.argmax([res['r2_avg'] for res in sweep_results])
    best_result = sweep_results[best_idx]
    return sweep_results, best_result

def select_temporal_order_tikhonov(u, u_clean, dt, dx, candidates, frac_space_orders, threshold, k_start, snr, true_alpha=None):
    """
    Performs Tikhonov-regularized pointwise model selection.
    """
    lambda_map = {100: 0.01, 80: 0.1, 60: 1.0, 40: 10.0, 30: 30.0, 20: 100.0, 10: 500.0}
    lambd = lambda_map.get(snr, 1.0)
    
    u_smooth = tikhonov_smooth(u, lambd)
    
    sweep, best = select_temporal_order(u_smooth, dt, dx, alpha_t_candidates=candidates, 
                                         frac_space_orders=frac_space_orders, threshold=threshold, k_start=k_start, u_clean=u_clean, true_alpha=true_alpha)
    return sweep, best

def select_temporal_order_ensemble(u, u_clean, dt, dx, candidates, frac_space_orders, threshold, k_start, true_alpha=None, B=15):
    """
    Performs Ensemble bagging model selection.
    """
    nt, nx, n_vars = u.shape
    nt_train = int(0.8 * nt)
    train_indices = np.arange(k_start, nt_train)
    n_train_select = int(0.8 * len(train_indices))
    
    selected_alphas = []
    
    for b in range(B):
        # Subsample indices with replacement
        sub_idx = np.random.choice(train_indices, size=n_train_select, replace=True)
        u_train_sub = u[sub_idx]
        
        best_alpha = None
        best_r2 = -np.inf
        
        for alpha_t in candidates:
            u_dot = gl_derivative_time(u, dt, alpha_t)
            u_dot_train_sub = u_dot[sub_idx]
            
            # Fit SINDy
            model, library = build_fractional_pde_model(u_train_sub, dt, dx, frac_space_orders, x_dot=u_dot_train_sub, threshold=threshold)
            
            # Predict on clean test set
            coefs = np.asarray(model.coefficients())
            u_clean_test = u_clean[nt_train:]
            Theta_test = library.transform(u_clean_test)
            Theta_test_flat = np.asarray(Theta_test).reshape(u_clean_test.shape[0] * u_clean_test.shape[1], -1)
            
            u_dot_clean = gl_derivative_time(u_clean, dt, true_alpha if true_alpha is not None else alpha_t)
            u_dot_test_clean = u_dot_clean[nt_train:]
            u_dot_test_clean_flat = u_dot_test_clean.reshape(u_dot_test_clean.shape[0] * u_dot_test_clean.shape[1], -1)
            
            u_dot_pred_flat = Theta_test_flat @ coefs.T
            
            r2_scores = []
            for i in range(n_vars):
                y_true = u_dot_test_clean_flat[:, i]
                y_pred = u_dot_pred_flat[:, i]
                ss_res = np.sum((y_true - y_pred)**2)
                ss_tot = np.sum((y_true - np.mean(y_true))**2)
                r2_scores.append(1.0 - ss_res / (ss_tot + 1e-10))
            r2_avg = np.mean(r2_scores)
            
            if r2_avg > best_r2:
                best_r2 = r2_avg
                best_alpha = alpha_t
                
        selected_alphas.append(best_alpha)
        
    final_alpha = np.median(selected_alphas)
    closest_candidate = candidates[np.argmin([abs(c - final_alpha) for c in candidates])]
    return closest_candidate

def main():
    print("=" * 60)
    print("OUROBOROS NOISE-MITIGATION COMPARISON SWEEP")
    print("=" * 60)
    
    # Spatial grid and initial conditions
    x = np.linspace(0, L, Nx)
    p0 = np.exp(-(x - L/2)**2 / 2.0)
    c0 = 1.0 - 0.5 * np.exp(-(x - L/2)**2 / 2.0)
    n0 = 0.1 + 0.2 * np.sin(np.pi * x / L)
    y0 = np.concatenate([p0, c0, n0])
    
    true_alphas = [0.5, 0.7, 0.9]
    candidates = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    frac_space_orders = (0.5, 1.5)
    threshold = 0.01
    k_start = 20
    
    snrs = [100, 80, 60, 40, 30, 20, 10]
    
    # Store errors for plotting: {alpha: {method: [errors_across_snrs]}}
    plot_data = {
        str(a): {
            'naive': [],
            'weak': [],
            'ensemble': [],
            'tikhonov': []
        } for a in true_alphas
    }
    
    mitigation_results = {}
    
    for alpha in true_alphas:
        print(f"\nRunning sweeps for true alpha = {alpha}...")
        u_clean = solve_fractional_system(alpha, Nt, dt, Nx, y0)
        
        mitigation_results[str(alpha)] = {
            'naive': {},
            'weak': {},
            'ensemble': {},
            'tikhonov': {}
        }
        
        for snr in snrs:
            print(f"  SNR = {snr} dB")
            np.random.seed(snr + 42)
            u_noisy = add_noise(u_clean, snr)
            
            # 1. Naive Pointwise
            _, best_naive = select_temporal_order(u_noisy, dt, dx, alpha_t_candidates=candidates,
                                                   frac_space_orders=frac_space_orders, threshold=threshold, k_start=k_start, u_clean=u_clean, true_alpha=alpha)
            err_naive = abs(best_naive['alpha_t'] - alpha)
            plot_data[str(alpha)]['naive'].append(err_naive)
            mitigation_results[str(alpha)]['naive'][str(snr)] = {
                'selected_alpha': float(best_naive['alpha_t']),
                'error': float(err_naive),
                'r2_avg': float(best_naive['r2_avg'])
            }
            
            # 2. Weak-form
            _, best_weak = select_temporal_order_weak(u_noisy, u_clean, dt, dx, candidates=candidates,
                                                      frac_space_orders=frac_space_orders, threshold=threshold, k_start=k_start, true_alpha=alpha)
            err_weak = abs(best_weak['alpha_t'] - alpha)
            plot_data[str(alpha)]['weak'].append(err_weak)
            mitigation_results[str(alpha)]['weak'][str(snr)] = {
                'selected_alpha': float(best_weak['alpha_t']),
                'error': float(err_weak),
                'r2_avg': float(best_weak['r2_avg'])
            }
            
            # 3. Ensemble
            selected_ensemble = select_temporal_order_ensemble(u_noisy, u_clean, dt, dx, candidates=candidates,
                                                               frac_space_orders=frac_space_orders, threshold=threshold, k_start=k_start, true_alpha=alpha, B=5)
            err_ensemble = abs(selected_ensemble - alpha)
            plot_data[str(alpha)]['ensemble'].append(err_ensemble)
            mitigation_results[str(alpha)]['ensemble'][str(snr)] = {
                'selected_alpha': float(selected_ensemble),
                'error': float(err_ensemble)
            }
            
            # 4. Tikhonov
            _, best_tikhonov = select_temporal_order_tikhonov(u_noisy, u_clean, dt, dx, candidates=candidates,
                                                               frac_space_orders=frac_space_orders, threshold=threshold, k_start=k_start, snr=snr, true_alpha=alpha)
            err_tikhonov = abs(best_tikhonov['alpha_t'] - alpha)
            plot_data[str(alpha)]['tikhonov'].append(err_tikhonov)
            mitigation_results[str(alpha)]['tikhonov'][str(snr)] = {
                'selected_alpha': float(best_tikhonov['alpha_t']),
                'error': float(err_tikhonov),
                'r2_avg': float(best_tikhonov['r2_avg'])
            }
            
            print(f"    Naive   -> Selected: {best_naive['alpha_t']:.2f} (error: {err_naive:.2f})")
            print(f"    Weak    -> Selected: {best_weak['alpha_t']:.2f} (error: {err_weak:.2f})")
            print(f"    Ensemble-> Selected: {selected_ensemble:.2f} (error: {err_ensemble:.2f})")
            print(f"    Tikhonov-> Selected: {best_tikhonov['alpha_t']:.2f} (error: {err_tikhonov:.2f})")

    # Save results to JSON
    with open('data/mitigation_results.json', 'w') as f:
        json.dump(mitigation_results, f, indent=2)
    print("\nSaved raw mitigation results to data/mitigation_results.json")

    # Generate Figure
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    methods = ['naive', 'weak', 'ensemble', 'tikhonov']
    markers = {'naive': 'o', 'weak': 's', 'ensemble': '^', 'tikhonov': 'd'}
    colors = {'naive': '#d62728', 'weak': '#1f77b4', 'ensemble': '#2ca02c', 'tikhonov': '#9467bd'}
    labels = {'naive': 'Naive Pointwise GL', 'weak': 'Weak-Form GL', 'ensemble': 'Ensemble-SINDy', 'tikhonov': 'Tikhonov-Regularized GL'}
    
    for idx, alpha in enumerate(true_alphas):
        ax = axes[idx]
        for m in methods:
            ax.plot(snrs, plot_data[str(alpha)][m], marker=markers[m], color=colors[m], label=labels[m], linewidth=2)
        ax.set_title(f"True Order $\\alpha_t = {alpha}$")
        ax.set_xlabel("SNR (dB)")
        ax.set_ylabel("Order Recovery Error $|\\hat{\\alpha}_t - \\alpha_t|$")
        ax.set_xscale('log')
        ax.invert_xaxis()
        ax.set_xticks(snrs)
        ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
        ax.grid(True)
        ax.legend()
        
    plt.suptitle("Order Recovery Error vs SNR for Noise-Mitigation Methods", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig('figures/mitigation_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved comparison figure to figures/mitigation_comparison.png")
    
    # Print final summary statistics
    print("\n" + "="*60)
    print("BREAKDOWN SNR SUMMARY Table")
    print("="*60)
    print(f"{'True alpha':<10} | {'Naive (dB)':<10} | {'Weak (dB)':<10} | {'Ensemble (dB)':<13} | {'Tikhonov (dB)':<13}")
    print("-"*65)
    
    for alpha in true_alphas:
        breakdowns = {}
        for m in methods:
            breakdown_snr = None
            for snr in snrs:
                err = mitigation_results[str(alpha)][m][str(snr)]['error']
                # Correct recovery means error is 0
                if err < 1e-5:
                    breakdown_snr = snr
                else:
                    break
        print(f"{alpha:<10} | {breakdowns['naive']:<10} | {breakdowns['weak']:<10} | {breakdowns['ensemble']:<13} | {breakdowns['tikhonov']:<13}")
        
    print("\nJSON_START")
    print(json.dumps(mitigation_results, indent=2))
    print("JSON_END")

if __name__ == '__main__':
    main()
