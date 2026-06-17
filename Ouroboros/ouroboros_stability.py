import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar

def find_hss(rho, gam_n, Dp, Dc, Dn, gam_p, Sp, gam_c, Sc):
    """Solve for the homogeneous steady state (n*, p*, c*) in the range (0, 1.5)"""
    def residual(n):
        if n <= 0.0 or n >= 1.5:
            return 999.0
        p = (Sp * n) / (gam_p + Sp * n)
        if p == 0.0:
            return 999.0
        c = (Sc * n * (1.0 + p)) / (gam_c * p)
        # Solve for the vessel growth/regression balance
        return rho * (1.0 - n) * (c / (1.0 + c)) - gam_n * p
    
    try:
        # Bracket the root in (1e-6, 1.5 - 1e-6)
        sol = root_scalar(residual, bracket=[1e-6, 1.5 - 1e-6], method='brentq')
        n_val = sol.root
        p_val = (Sp * n_val) / (gam_p + Sp * n_val)
        c_val = (Sc * n_val * (1.0 + p_val)) / (gam_c * p_val)
        return n_val, p_val, c_val
    except ValueError:
        # If bracket has same signs, check if there's any valid root or return None
        return None

def compute_jacobian(n_val, p_val, c_val, rho, gam_n, gam_p, Sp, gam_c, Sc):
    """Construct the local Jacobian matrix at the homogeneous steady state"""
    J11 = -gam_p - Sp * n_val
    J12 = 0.0
    J13 = Sp * (1.0 - p_val)

    J21 = -gam_c * c_val / (1.0 + p_val)**2
    J22 = -gam_c * p_val / (1.0 + p_val)
    J23 = Sc

    J31 = -gam_n * n_val
    J32 = rho * n_val * (1.0 - n_val) / (1.0 + c_val)**2
    J33 = -rho * n_val * c_val / (1.0 + c_val)

    return np.array([
        [J11, J12, J13],
        [J21, J22, J23],
        [J31, J32, J33]
    ])

def analyze_stability(rho, gam_n, Dn, Dp=0.05, Dc=0.1, gam_p=0.05, Sp=0.1, gam_c=0.2, Sc=0.3):
    """Find HSS, compute eigenvalues over wavenumbers k, and classify stability"""
    hss = find_hss(rho, gam_n, Dp, Dc, Dn, gam_p, Sp, gam_c, Sc)
    if hss is None:
        return {
            'status': 'no_hss',
            'max_re': -999.0,
            'max_k': 0.0,
            'max_im': 0.0,
            'hss': None
        }
    
    n_val, p_val, c_val = hss
    J_local = compute_jacobian(n_val, p_val, c_val, rho, gam_n, gam_p, Sp, gam_c, Sc)
    
    k_vals = np.linspace(0.0, 30.0, 300)
    max_re = -np.inf
    max_k = 0.0
    max_im = 0.0
    
    # Compute eigenvalues for each wavenumber
    for k in k_vals:
        J_k = J_local.copy()
        J_k[0, 0] -= Dp * k**2
        J_k[1, 1] -= Dc * k**2
        J_k[2, 2] -= Dn * k**2
        
        eigenvals = np.linalg.eigvals(J_k)
        # Find eigenvalue with largest real part
        idx = np.argmax(np.real(eigenvals))
        re = np.real(eigenvals[idx])
        im = np.imag(eigenvals[idx])
        
        if re > max_re:
            max_re = re
            max_k = k
            max_im = im
            
    # Classify the instability type
    if max_re <= 0.0:
        status = 'stable'
    else:
        if abs(max_im) > 1e-5:
            status = 'hopf'
        else:
            # Turing if finite k is unstable and k=0 is stable
            J_k0 = J_local.copy()
            eigenvals_k0 = np.linalg.eigvals(J_k0)
            re_k0 = np.max(np.real(eigenvals_k0))
            if re_k0 <= 0.0 and max_k > 0.05:
                status = 'turing'
            else:
                status = 'stationary_homogeneous'
                
    return {
        'status': status,
        'max_re': max_re,
        'max_k': max_k,
        'max_im': max_im,
        'hss': (n_val, p_val, c_val)
    }

def main():
    print("="*60)
    print("OUROBOROS LINEAR STABILITY ANALYSIS")
    print("="*60)
    
    # Define grid parameter ranges (1000 points total)
    rho_vals = np.linspace(0.1, 5.0, 10)
    gam_n_vals = np.linspace(0.1, 5.0, 10)
    Dn_vals = np.linspace(0.001, 0.05, 10)
    
    results = []
    stable_count = 0
    hopf_count = 0
    turing_count = 0
    stationary_homogeneous_count = 0
    no_hss_count = 0
    
    max_growth_rate = -np.inf
    most_unstable_params = None
    
    # Run 3D sweep
    print(f"Sweeping 3D parameter grid: rho x gam_n x Dn ({len(rho_vals)}x{len(gam_n_vals)}x{len(Dn_vals)} = 1000 points)...")
    for r in rho_vals:
        for gn in gam_n_vals:
            for dn in Dn_vals:
                res = analyze_stability(r, gn, dn)
                res_dict = {
                    'rho': float(r),
                    'gam_n': float(gn),
                    'Dn': float(dn),
                    'status': res['status'],
                    'max_re': float(res['max_re']),
                    'max_k': float(res['max_k']),
                    'max_im': float(res['max_im']),
                }
                if res['hss'] is not None:
                    res_dict['hss'] = [float(val) for val in res['hss']]
                else:
                    res_dict['hss'] = None
                    
                results.append(res_dict)
                
                # Count statuses
                if res['status'] == 'stable':
                    stable_count += 1
                elif res['status'] == 'hopf':
                    hopf_count += 1
                elif res['status'] == 'turing':
                    turing_count += 1
                elif res['status'] == 'stationary_homogeneous':
                    stationary_homogeneous_count += 1
                elif res['status'] == 'no_hss':
                    no_hss_count += 1
                    
                if res['max_re'] > max_growth_rate:
                    max_growth_rate = res['max_re']
                    most_unstable_params = res_dict
                    
    print("\nSweep Statistics:")
    print(f"  Stable points:                 {stable_count} / 1000")
    print(f"  Hopf-type points:              {hopf_count} / 1000")
    print(f"  Turing-type points:            {turing_count} / 1000")
    print(f"  Stationary homogeneous points: {stationary_homogeneous_count} / 1000")
    print(f"  No HSS points:                 {no_hss_count} / 1000")
    print(f"Maximum leading Re(lambda):       {max_growth_rate:.6f}")
    if most_unstable_params:
        print(f"Most unstable parameters: rho={most_unstable_params['rho']:.4f}, gam_n={most_unstable_params['gam_n']:.4f}, Dn={most_unstable_params['Dn']:.4f}")
        print(f"  Status: {most_unstable_params['status']}, max_k: {most_unstable_params['max_k']:.4f}, max_im: {most_unstable_params['max_im']:.4f}")
        
    # Save raw results
    os.makedirs('data', exist_ok=True)
    with open('data/stability_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("\nSaved stability results to data/stability_results.json")
    
    # Generate multi-panel 2D stability map plots
    # We will pick 3 slices of Dn: minimum (0.001), mid (0.022), and maximum (0.05)
    selected_Dns = [Dn_vals[0], Dn_vals[len(Dn_vals)//2], Dn_vals[-1]]
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    for idx, dn_slice in enumerate(selected_Dns):
        ax = axes[idx]
        
        # Extract 2D grid of rho vs gam_n for this Dn
        grid_data = np.zeros((len(rho_vals), len(gam_n_vals)))
        status_colors = np.zeros((len(rho_vals), len(gam_n_vals)))
        
        for r_idx, r in enumerate(rho_vals):
            for gn_idx, gn in enumerate(gam_n_vals):
                # Find matching result
                match = next((item for item in results if abs(item['rho'] - r) < 1e-5 and abs(item['gam_n'] - gn) < 1e-5 and abs(item['Dn'] - dn_slice) < 1e-5), None)
                if match:
                    grid_data[r_idx, gn_idx] = match['max_re']
                    
        # Create a contour/heatmap
        im = ax.imshow(
            grid_data, 
            extent=[gam_n_vals[0], gam_n_vals[-1], rho_vals[0], rho_vals[-1]], 
            origin='lower',
            cmap='coolwarm',
            aspect='auto',
            vmin=-0.1, vmax=0.1
        )
        
        ax.set_title(f"Dn = {dn_slice:.4f} (Dn/Dp = {dn_slice/0.05:.2f})")
        ax.set_xlabel("Vessel Regression Rate (gam_n)")
        if idx == 0:
            ax.set_ylabel("Vessel Growth Rate (rho)")
        ax.grid(True, linestyle='--', alpha=0.5)
        fig.colorbar(im, ax=ax, label=r"Leading $\text{Re}(\lambda)$")
        
    plt.suptitle("OUROBOROS Linear Stability Map (Leading Growth Rate)", fontsize=14, y=1.02)
    plt.tight_layout()
    os.makedirs('figures', exist_ok=True)
    plt.savefig('figures/stability_map.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved stability map plot to figures/stability_map.png")

    # Checkpoint 2: Hardened stability edge checks
    print("\n" + "="*60)
    print("CHECKPOINT 2: HARDENING STABILITY CLAIM")
    print("="*60)
    
    # 1. Past-the-corner points
    past_points = [
        {'rho': 0.05, 'gam_n': 5.0, 'Dn': 0.001, 'label': 'Lower growth rate (rho=0.05)'},
        {'rho': 0.1,  'gam_n': 6.0, 'Dn': 0.001, 'label': 'Higher regression rate (gam_n=6.0)'},
        {'rho': 0.1,  'gam_n': 5.0, 'Dn': 0.0005, 'label': 'Lower vessel diffusion (Dn=0.0005)'},
        {'rho': 0.05, 'gam_n': 6.0, 'Dn': 0.0005, 'label': 'Extrapolated corner (rho=0.05, gam_n=6.0, Dn=0.0005)'}
    ]
    
    print("\nEvaluating linear stability just past the swept-box corner:")
    hardened_results = []
    for pt in past_points:
        res = analyze_stability(pt['rho'], pt['gam_n'], pt['Dn'])
        print(f"  {pt['label']}:")
        print(f"    Leading Re(lambda) = {res['max_re']:.6f} | status: {res['status']}")
        if res['hss']:
            print(f"    HSS (n*, p*, c*): ({res['hss'][0]:.4f}, {res['hss'][1]:.4f}, {res['hss'][2]:.4f})")
        hardened_results.append({
            'rho': pt['rho'],
            'gam_n': pt['gam_n'],
            'Dn': pt['Dn'],
            'label': pt['label'],
            'max_re': res['max_re'],
            'status': res['status']
        })

    # 2. Non-swept parameters variations
    # Baseline: gam_c = 0.2, Sc = 0.3
    # We test variations of gam_c in [0.05, 0.5] and Sc in [0.1, 0.6]
    gam_c_vals = [0.05, 0.5]
    Sc_vals = [0.1, 0.6]
    
    # We test these variations at two physical regimes: 
    # (i) baseline swept params (rho=0.2, gam_n=0.1, Dn=0.01)
    # (ii) corner swept params (rho=0.1, gam_n=5.0, Dn=0.001)
    regimes = [
        {'rho': 0.2, 'gam_n': 0.1, 'Dn': 0.01, 'name': 'Baseline Swept Params'},
        {'rho': 0.1, 'gam_n': 5.0, 'Dn': 0.001, 'name': 'Corner Swept Params'}
    ]
    
    print("\nVarying non-swept parameters (gam_c, Sc):")
    nonswept_results = []
    for regime in regimes:
        print(f"  In {regime['name']} (rho={regime['rho']}, gam_n={regime['gam_n']}, Dn={regime['Dn']}):")
        for gc in gam_c_vals:
            for sc in Sc_vals:
                res = analyze_stability(
                    regime['rho'], regime['gam_n'], regime['Dn'],
                    gam_c=gc, Sc=sc
                )
                print(f"    gam_c={gc:.2f}, Sc={sc:.2f} -> Leading Re(lambda) = {res['max_re']:.6f} | status: {res['status']}")
                nonswept_results.append({
                    'regime': regime['name'],
                    'rho': regime['rho'],
                    'gam_n': regime['gam_n'],
                    'Dn': regime['Dn'],
                    'gam_c': gc,
                    'Sc': sc,
                    'max_re': res['max_re'],
                    'status': res['status']
                })
                
    # Save the hardened results to stability_hardened.json
    hardened_data = {
        'past_corner_points': hardened_results,
        'non_swept_variations': nonswept_results
    }
    with open('data/stability_hardened.json', 'w') as f:
        json.dump(hardened_data, f, indent=2)
    print("\nSaved hardened stability results to data/stability_hardened.json")
    
if __name__ == '__main__':
    main()

