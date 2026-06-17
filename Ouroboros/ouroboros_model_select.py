import warnings
import numpy as np
from ouroboros_fractional_sindy import gl_derivative_time, build_fractional_pde_model

# Suppress numerical warnings from singular matrices/convergences
warnings.filterwarnings('ignore')

def select_temporal_order(u, dt, dx, alpha_t_candidates=(0.6, 0.8, 1.0), frac_space_orders=(0.5, 1.5), threshold=0.05, k_start=10, u_clean=None, true_alpha=None):
    """
    Performs a held-out R2 sweep over temporal fractional orders.
    Splits the time series: first 80% for training, remaining 20% for testing.
    Uses k_start index to discard the initial step-change artifact in GL derivative.
    If u_clean is provided, test evaluation (R2) is scored against the clean reference.
    If true_alpha is also provided, the clean reference target is computed at true_alpha.
    """
    nt, nx, n_vars = u.shape
    nt_train = int(0.8 * nt)
    
    u_train = u[:nt_train]
    u_test = u[nt_train:]
    
    sweep_results = []
    
    for alpha_t in alpha_t_candidates:
        # 1. Compute time derivative for the entire time series (to prevent boundary/truncation errors)
        u_dot = gl_derivative_time(u, dt, alpha_t)
        
        u_dot_train = u_dot[:nt_train]
        u_dot_test = u_dot[nt_train:]
        
        # Slice training data to avoid the initial step change spike at t=0
        u_train_slice = u_train[k_start:]
        u_dot_train_slice = u_dot_train[k_start:]
        
        # 2. Fit the fractional SINDy model on the sliced training set
        model, library = build_fractional_pde_model(u_train_slice, dt, dx, frac_space_orders, x_dot=u_dot_train_slice, threshold=threshold)
        
        # 3. Predict on the test set
        # Get coefficients: shape (n_vars, n_features)
        coefs = np.asarray(model.coefficients())
        
        # If u_clean is provided, compute features and targets from u_clean (test set)
        if u_clean is not None:
            u_clean_test = u_clean[nt_train:]
            Theta_test = library.transform(u_clean_test)
            Theta_test_flat = np.asarray(Theta_test).reshape(u_clean_test.shape[0] * u_clean_test.shape[1], -1)
            
            target_alpha = true_alpha if true_alpha is not None else alpha_t
            u_dot_clean = gl_derivative_time(u_clean, dt, target_alpha)
            u_dot_test_clean = u_dot_clean[nt_train:]
            u_dot_test_clean_flat = u_dot_test_clean.reshape(u_dot_test_clean.shape[0] * u_dot_test_clean.shape[1], -1)
            u_dot_test_flat = u_dot_test_clean_flat
        else:
            # Compute design matrix on test set: shape (nt_test, nx, n_features)
            Theta_test = library.transform(u_test)
            Theta_test_flat = np.asarray(Theta_test).reshape(u_test.shape[0] * u_test.shape[1], -1)
            u_dot_test_flat = u_dot_test.reshape(u_dot_test.shape[0] * u_dot_test.shape[1], -1)
        
        # Compute predictions: shape (nt_test * nx, n_vars)
        u_dot_pred_flat = Theta_test_flat @ coefs.T
        
        # 4. Calculate R2 score for each variable
        r2_scores = []
        for i in range(n_vars):
            y_true = u_dot_test_flat[:, i]
            y_pred = u_dot_pred_flat[:, i]
            
            ss_res = np.sum((y_true - y_pred)**2)
            ss_tot = np.sum((y_true - np.mean(y_true))**2)
            
            r2 = 1.0 - (ss_res / (ss_tot + 1e-10))
            r2_scores.append(r2)
            
        r2_avg = np.mean(r2_scores)
        n_active = np.sum(coefs != 0.0)
        
        sweep_results.append({
            'alpha_t': alpha_t,
            'r2_avg': r2_avg,
            'r2_scores': r2_scores,
            'n_active': n_active,
            'model': model,
            'library': library,
            'u_dot_test_flat': u_dot_test_flat,
            'u_dot_pred_flat': u_dot_pred_flat
        })
        
    # Select best alpha_t based on highest average test R2
    best_idx = np.argmax([res['r2_avg'] for res in sweep_results])
    best_result = sweep_results[best_idx]
    
    return sweep_results, best_result

if __name__ == '__main__':
    # Small test
    import json
    data = np.load('data/ouroboros_synth.npz')
    u = data['u']
    dt = data['t'][1] - data['t'][0]
    dx = data['x'][1] - data['x'][0]
    
    results, best = select_temporal_order(u, dt, dx)
    print("\nSweep Table:")
    print(f"{'alpha_t':<10} | {'test_R2':<10} | {'n_active':<10}")
    print("-" * 38)
    for res in results:
        print(f"{res['alpha_t']:<10.2f} | {res['r2_avg']:<10.4f} | {res['n_active']:<10}")
        
    print("\nBest temporal order:", best['alpha_t'])
