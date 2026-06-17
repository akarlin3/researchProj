import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

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
    
    # Clip to prevent negative concentrations and runaway growth
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

def run_simulation(Nx=100, Nt=200, L=10.0, T=5.0, seed=42):
    np.random.seed(seed)
    dx = L / (Nx - 1)
    t_eval = np.linspace(0, T, Nt)
    dt = T / (Nt - 1)
    
    x = np.linspace(0, L, Nx)
    
    # Initial Conditions
    p0 = np.exp(-(x - L/2)**2 / 2.0)
    c0 = 1.0 - 0.5 * np.exp(-(x - L/2)**2 / 2.0)
    n0 = 0.1 + 0.2 * np.sin(np.pi * x / L)
    y0 = np.concatenate([p0, c0, n0])
    
    # Parameters
    params = {
        'Dp': 0.05,
        'Dc': 0.1,
        'Dn': 0.01,
        'gam_p': 0.05,
        'Sp': 0.1,
        'gam_c': 0.2,
        'Sc': 0.3,
        'rho': 0.2,
        'gam_n': 0.1
    }
    
    # Solve stiff ODE system
    sol = solve_ivp(
        stroma_rhs,
        [0, T],
        y0,
        t_eval=t_eval,
        method='Radau',
        args=(Nx, dx, params['Dp'], params['Dc'], params['Dn'], params['gam_p'], params['Sp'], params['gam_c'], params['Sc'], params['rho'], params['gam_n'])
    )
    
    # Extract solutions
    p_sol = sol.y[:Nx, :].T  # shape (Nt, Nx)
    c_sol = sol.y[Nx:2*Nx, :].T
    n_sol = sol.y[2*Nx:, :].T
    
    u = np.stack([p_sol, c_sol, n_sol], axis=-1)  # shape (Nt, Nx, 3)
    
    # Create output directories
    os.makedirs('data', exist_ok=True)
    os.makedirs('figures', exist_ok=True)
    
    # Save npz
    metadata = {
        'equations': {
            'p': 'dp/dt = Dp*p_xx - gam_p*p + Sp*n*(1-p)',
            'c': 'dc/dt = Dc*c_xx - gam_c*c*p/(1+p) + Sc*n',
            'n': 'dn/dt = Dn*n_xx + rho*n*(1-n)*c/(1+c) - gam_n*n*p'
        },
        'parameters': params,
        'grid': {'Nx': Nx, 'Nt': Nt, 'L': L, 'T': T}
    }
    
    np.savez(
        'data/ouroboros_synth.npz',
        t=sol.t,
        x=x,
        u=u,
        metadata=json.dumps(metadata)
    )
    
    # Create sanity plot
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    plot_times = [0, int(Nt*0.25), int(Nt*0.5), int(Nt*0.75), Nt-1]
    
    for idx, (ax, name, color) in enumerate(zip(axes, ['Pressure (p)', 'Oxygen (c)', 'Vessel Density (n)'], ['red', 'blue', 'green'])):
        for t_idx in plot_times:
            ax.plot(x, u[t_idx, :, idx], label=f't={sol.t[t_idx]:.2f}')
        ax.set_title(name)
        ax.set_xlabel('Space (x)')
        ax.set_ylabel('Value')
        ax.grid(True)
        if idx == 0:
            ax.legend()
            
    plt.tight_layout()
    plt.savefig('figures/sim_fields.png', dpi=150)
    plt.close()
    
    print(f"Simulation completed. Saved data to data/ouroboros_synth.npz and plot to figures/sim_fields.png")

if __name__ == '__main__':
    run_simulation()
