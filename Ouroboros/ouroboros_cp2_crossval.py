"""
Checkpoint 2: Cross-validate OLD vs NEW pipeline on a shared disagreement cell.

Cell: alpha_t = 0.9, weak-form, SNR = 10 dB.
  * OLD (mitigation.py) single-realization result: selected 0.9, error 0.0
    -> "weak rescues alpha=0.9 down to 10 dB".
  * NEW (fine sweep, 500 realizations): success rate 0.494 at 10 dB
    -> "10 dB is a coin-flip; true bracket [15,20]".

We run, on the SAME 500 noise draws (seed = 10+42 = 52, identical draw order):
  (A) OLD metric  -> ouroboros_mitigation.select_temporal_order_weak  (ps.STLSQ)
  (B) NEW metric  -> ouroboros_fine_snr_sweep.select_temporal_order_weak_fast (fast_stlsq)
and also report the OLD single-realization (first draw) outcome to reproduce the
original "10 dB success" claim.

This isolates: is the flip driven by the STLSQ implementation (A vs B differ over
the same 500 draws) or purely by realization count (A==B over 500, but N=1 was a
lucky draw)?
"""
import numpy as np

from ouroboros_identifiability import solve_fractional_system, add_noise
from ouroboros_fractional_sindy import gl_weights
from ouroboros_mitigation import select_temporal_order_weak  # OLD (ps.STLSQ)
from ouroboros_fine_snr_sweep import (
    Nx, Nt, L, dt, dx, k_start,
    get_features, get_phi_matrix, get_psi_matrix,
    select_temporal_order_weak_fast,  # NEW (fast_stlsq)
)

np.seterr(all="ignore")

ALPHA = 0.9
SNR = 10
N_trials = 500
candidates = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
frac_space_orders = (0.5, 1.5)
threshold = 0.01


def weak_precompute(u_clean, alpha):
    nt_train = int(0.8 * Nt)
    H, S = 50, 5
    train_windows, test_windows = [], []
    for k_a in range(0, Nt - H + 1, S):
        k_b = k_a + H
        if k_b < nt_train:
            if k_a >= k_start:
                train_windows.append((k_a, k_b))
        elif k_a >= nt_train:
            if k_b < Nt:
                test_windows.append((k_a, k_b))
    phi = (1.0 - ((2.0 * np.arange(H + 1) / H) - 1.0) ** 2) ** 4
    phi /= np.sum(phi) * dt
    Phi_mat_train = get_phi_matrix(nt_train, train_windows, phi)
    Phi_mat_test = get_phi_matrix(
        Nt - nt_train, [(k_a - nt_train, k_b - nt_train) for k_a, k_b in test_windows], phi
    )
    Theta_clean, _ = get_features(u_clean, dx, frac_space_orders)
    Theta_clean_test = Theta_clean[nt_train:]
    X_test_temp = np.tensordot(Theta_clean_test, Phi_mat_test, axes=([0], [0])) * dt
    X_test_weak = np.transpose(X_test_temp, (2, 0, 1)).reshape(-1, Theta_clean.shape[-1])
    w_test = gl_weights(alpha, Nt)
    Psi_test = get_psi_matrix(Nt, test_windows, alpha, w_test, phi)
    Y_test_temp = np.tensordot(u_clean, Psi_test, axes=([0], [0]))
    Y_test_weak = np.transpose(Y_test_temp, (2, 0, 1)).reshape(-1, 3) * (dt / (dt ** ALPHA))
    return Phi_mat_train, X_test_weak, Y_test_weak


def main():
    x = np.linspace(0, L, Nx)
    p0 = np.exp(-(x - L / 2) ** 2 / 2.0)
    c0 = 1.0 - 0.5 * np.exp(-(x - L / 2) ** 2 / 2.0)
    n0 = 0.1 + 0.2 * np.sin(np.pi * x / L)
    y0 = np.concatenate([p0, c0, n0])
    u_clean = solve_fractional_system(ALPHA, Nt, dt, Nx, y0)
    Phi_mat_train, X_test_weak, Y_test_weak = weak_precompute(u_clean, ALPHA)

    # Identical 500 draws, same seed/order as the production sweeps.
    np.random.seed(int(SNR) + 42)
    noisy = [add_noise(u_clean, SNR) for _ in range(N_trials)]

    old_sel = np.empty(N_trials)
    new_sel = np.empty(N_trials)
    for i, u_noisy in enumerate(noisy):
        _, best_old = select_temporal_order_weak(
            u_noisy, u_clean, dt, dx, candidates=candidates,
            frac_space_orders=frac_space_orders, threshold=threshold,
            k_start=k_start, true_alpha=ALPHA,
        )
        old_sel[i] = best_old["alpha_t"]
        new_sel[i] = select_temporal_order_weak_fast(
            u_noisy, u_clean, dt, dx, candidates, frac_space_orders,
            threshold, k_start, ALPHA, Phi_mat_train, X_test_weak, Y_test_weak,
        )[1]["alpha_t"]
        if (i + 1) % 100 == 0:
            print(f"  ...{i+1}/{N_trials}")

    def summarize(name, sel):
        err = np.abs(sel - ALPHA)
        succ = err < 1e-5
        rate = float(np.mean(succ))
        fe = err[~succ]
        print(f"{name}: success={rate:.3f} ({int(succ.sum())}/{N_trials}) | "
              f"mean_err={err.mean():.4f} | failed-trial scatter={fe.mean() if fe.size else 0:.4f}"
              f"±{fe.std() if fe.size else 0:.4f}")
        return rate

    print("\n========== CHECKPOINT 2: alpha=0.9 weak @ 10 dB ==========")
    print(f"\n(1) OLD single realization (first draw, seed {SNR+42}): "
          f"selected = {old_sel[0]:.2f}  (reproduces mitigation.py N=1 -> error {abs(old_sel[0]-ALPHA):.2f})")
    print(f"    NEW fast on that same first draw:        selected = {new_sel[0]:.2f}")
    print("\n(2) Over the SAME 500 draws:")
    r_old = summarize("    OLD ps.STLSQ metric ", old_sel)
    r_new = summarize("    NEW fast_stlsq metric", new_sel)
    agree = float(np.mean(old_sel == new_sel))
    print(f"\n    Per-realization agreement (old vs new selected alpha): {agree:.3f}")
    print(f"    => OLD 500-rate={r_old:.3f}  NEW 500-rate={r_new:.3f}")

    import json
    json.dump(
        {"old_sel": old_sel.tolist(), "new_sel": new_sel.tolist(),
         "old_rate": r_old, "new_rate": r_new, "agreement": agree,
         "old_N1_first": float(old_sel[0]), "new_N1_first": float(new_sel[0])},
        open("data/cp2_crossval.json", "w"), indent=2,
    )
    print("\nSaved data/cp2_crossval.json")


if __name__ == "__main__":
    main()
