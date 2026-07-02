import os
import json
import numpy as np
import matplotlib.pyplot as plt
from ouroboros_model_select import select_temporal_order
from ouroboros_fractional_sindy import gl_derivative_time

def run_discovery_pipeline():
    # 1. Check if synthetic dataset exists; if not, run simulation
    if not os.path.exists('data/ouroboros_synth.npz') or not os.path.exists('figures/sim_fields.png'):
        print("Dataset or sanity figure missing. Running simulation first...")
        from ouroboros_sim import run_simulation
        run_simulation()
        
    # 2. Load dataset
    data = np.load('data/ouroboros_synth.npz')
    t = data['t']
    x = data['x']
    u = data['u']
    metadata = json.loads(data['metadata'].item())
    
    dt = t[1] - t[0]
    dx = x[1] - x[0]
    
    print(f"Loaded simulated fields of shape {u.shape}.")
    print(f"Grid details: Nx = {len(x)}, Nt = {len(t)}, dx = {dx:.4f}, dt = {dt:.4f}")
    
    # 3. Perform temporal order selection sweep
    alpha_t_candidates = (0.6, 0.8, 1.0)
    frac_space_orders = (0.5, 1.5)
    
    print("\nRunning Fractional-SINDy temporal-order sweep...")
    results, best = select_temporal_order(u, dt, dx, alpha_t_candidates, frac_space_orders, threshold=0.05, k_start=10)
    
    # 4. Print Comparison Table
    print("\n" + "="*50)
    print(f"{'alpha_t':<10} | {'test_R2 (Avg)':<15} | {'n_active':<10}")
    print("-"*50)
    for res in results:
        print(f"{res['alpha_t']:<10.2f} | {res['r2_avg']:<15.5f} | {res['n_active']:<10}")
    print("="*50)
    
    print(f"\nSelected temporal order (highest held-out R2): alpha_t = {best['alpha_t']:.2f}")
    print("Discovered SINDy Equations for selected model:")
    best['model'].print()
    
    # 5. Plot LHS Validation (Predicted vs True LHS derivatives)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    var_names = ['Pressure (p)', 'Oxygen (c)', 'Vessel Density (n)']
    colors = ['red', 'blue', 'green']
    
    u_dot_test_flat = best['u_dot_test_flat']
    u_dot_pred_flat = best['u_dot_pred_flat']
    
    for i in range(3):
        ax = axes[i]
        y_true = u_dot_test_flat[:, i]
        y_pred = u_dot_pred_flat[:, i]
        
        ax.scatter(y_true, y_pred, color=colors[i], alpha=0.3, s=15, label='Test points')
        # Diagonal reference line
        min_val = min(y_true.min(), y_pred.min())
        max_val = max(y_true.max(), y_pred.max())
        ax.plot([min_val, max_val], [min_val, max_val], 'k--', label='True LHS = Pred LHS')
        
        # Calculate individual R^2
        r2 = best['r2_scores'][i]
        
        ax.set_title(f"{var_names[i]}\n$R^2$ = {r2:.4f}")
        ax.set_xlabel('True LHS derivative')
        ax.set_ylabel('Predicted LHS derivative')
        ax.grid(True)
        ax.legend()
        
    plt.suptitle(f"LHS Derivative Validation (Test Set) — Selected temporal order $\\alpha_t$ = {best['alpha_t']:.1f}", fontsize=14)
    plt.tight_layout()
    os.makedirs('figures', exist_ok=True)
    plt.savefig('figures/lhs_validation.png', dpi=150)
    plt.close()
    print("Saved LHS validation plot to figures/lhs_validation.png")
    
    # 6. Generate Refutation Report: RESULTS_phase1.md
    report_content = f"""# Phase 1 SINDy Temporal-Order Selection & Refutation Report

This report summarizes the data-driven discovery and temporal order selection results on simulated stromal-pressure dynamics.

## Model Selection Sweep Results

Below is the sweep of temporal fractional orders $\\alpha_t$ against their generalization performance (held-out $R^2$) and complexity (number of active terms in the model):

| Temporal Order ($\\alpha_t$) | Held-out $R^2$ (Average) | Number of Active Terms |
| :--- | :--- | :--- |
"""
    for res in results:
        report_content += f"| {res['alpha_t']:.1f} | {res['r2_avg']:.5f} | {res['n_active']} |\n"
        
    # Check if alpha_t = 1.0 won
    alpha_1_r2 = [r['r2_avg'] for r in results if r['alpha_t'] == 1.0][0]
    best_r2 = best['r2_avg']
    alpha_1_won = (best['alpha_t'] == 1.0)
    
    report_content += f"""
### Key Findings
- **Selected Temporal Order**: $\\alpha_t$ = {best['alpha_t']:.1f} (Held-out $R^2$ = {best_r2:.5f})
- **Did $\\alpha_t = 1.0$ Win?**: **{"YES" if alpha_1_won else "NO"}**
- **Discovered Active Term Count**: {best['n_active']}

## Refutation Conclusion

Ordinary integer-time dynamics ($\\alpha_t = 1.0$) **{"yielded the highest generalization score" if alpha_1_won else "did NOT yield the highest generalization score"}** on the held-out test split. 

Since the simulator's governing law is strictly of integer-order ($\\alpha_t = 1.0$), **{"the fractional SINDy pipeline successfully refuted fractional dynamics and recovered the true integer-order dynamics" if alpha_1_won else "the pipeline selected a fractional derivative model, indicating a potential overfitting or model selection failure"}**. This provides a strong **{"validation" if alpha_1_won else "refutation"}** of the fractional temporal formulation.

## LHS Prediction Validation Plot

The scatter plot below compares the true derivative value against the SINDy prediction for each of the three coupled variables on the held-out test split:

![LHS Validation Plot](figures/lhs_validation.png)

## Discovered Equations for Selected Model ($\\alpha_t$ = {best['alpha_t']:.1f})

```text
"""
    # Capture model print output by temporarily hijacking stdout
    import io, sys
    old_stdout = sys.stdout
    sys.stdout = mystdout = io.StringIO()
    best['model'].print()
    sys.stdout = old_stdout
    
    report_content += mystdout.getvalue()
    report_content += """```
"""
    
    with open('RESULTS_phase1.md', 'w') as f:
        f.write(report_content)
    print("Saved refutation report to RESULTS_phase1.md")

if __name__ == '__main__':
    run_discovery_pipeline()
