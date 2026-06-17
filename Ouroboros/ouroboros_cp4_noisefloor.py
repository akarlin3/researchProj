"""
Checkpoint 4: noise-floor selection bias + Wilson confidence intervals.

(A) For each true alpha, compute the order-SELECTION DISTRIBUTION under near-pure noise
    (-5, 0, 5 dB) for the oracle weak and pointwise selectors. This shows the selector's
    bias toward high candidate orders that explains why alpha=0.9 reaches ~48% apparent
    "recovery" at 0 dB while uniform chance over 9 candidates is only ~11%: at the noise
    floor the selector concentrates mass on high orders, so a high TRUE order coincides
    with the bias (inflated success), a low true order does not.

(B) Wilson 95% score intervals for every bracket-edge success rate across all
    systems/methods, drawn from the CP1/CP2/CP3 result JSONs plus the existing oracle
    fine sweep, especially the alpha=0.9 weak @15 dB edge (0.944). States whether each
    fail/succeed edge is robust or criterion-fragile (CI vs the 95% line).
"""
import json
import os
import numpy as np
from multiprocessing import Pool

from ouroboros_identifiability import solve_fractional_system, add_noise
from ouroboros_fractional_sindy import gl_weights
from ouroboros_fine_snr_sweep import (
    Nx, Nt, T, L, dt, dx, k_start,
    get_features, get_phi_matrix, get_psi_matrix, gl_derivative_time,
    select_temporal_order_pointwise_fast, select_temporal_order_weak_fast,
)

np.seterr(all="ignore")

TRUE_ALPHAS = [0.5, 0.7, 0.9]
CANDIDATES = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
FRAC = (0.5, 1.5)
THR = 0.01
N_TRIALS = 500
FLOOR_SNRS = [-5, 0, 5]
Z = 1.959963984540054  # 95%


def _y0():
    x = np.linspace(0, L, Nx)
    p0 = np.exp(-(x - L / 2) ** 2 / 2.0)
    c0 = 1.0 - 0.5 * np.exp(-(x - L / 2) ** 2 / 2.0)
    n0 = 0.1 + 0.2 * np.sin(np.pi * x / L)
    return np.concatenate([p0, c0, n0])


def wilson(k, n, z=Z):
    if n == 0:
        return (0.0, 0.0, 0.0)
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (p, max(0.0, center - half), min(1.0, center + half))


def _precompute_oracle(u_clean, alpha):
    nt_train = int(0.8 * Nt)
    Theta_test, _ = get_features(u_clean[nt_train:], dx, FRAC)
    Theta_test_flat = Theta_test.reshape(-1, Theta_test.shape[-1])
    u_dot_clean = gl_derivative_time(u_clean, dt, alpha)
    u_dot_test_clean_flat = u_dot_clean[nt_train:].reshape(-1, 3)
    H, S = 50, 5
    tr, te = [], []
    for k_a in range(0, Nt - H + 1, S):
        k_b = k_a + H
        if k_b < nt_train:
            if k_a >= k_start:
                tr.append((k_a, k_b))
        elif k_a >= nt_train:
            if k_b < Nt:
                te.append((k_a, k_b))
    phi = (1.0 - ((2.0 * np.arange(H + 1) / H) - 1.0) ** 2) ** 4
    phi /= np.sum(phi) * dt
    Phi_mat_train = get_phi_matrix(nt_train, tr, phi)
    Phi_mat_test = get_phi_matrix(Nt - nt_train, [(a - nt_train, b - nt_train) for a, b in te], phi)
    Theta_clean, _ = get_features(u_clean, dx, FRAC)
    X_test_temp = np.tensordot(Theta_clean[nt_train:], Phi_mat_test, axes=([0], [0])) * dt
    X_test_weak = np.transpose(X_test_temp, (2, 0, 1)).reshape(-1, Theta_clean.shape[-1])
    w_test = gl_weights(alpha, Nt)
    Psi_test = get_psi_matrix(Nt, te, alpha, w_test, phi)
    Y_test_temp = np.tensordot(u_clean, Psi_test, axes=([0], [0]))
    Y_test_weak = np.transpose(Y_test_temp, (2, 0, 1)).reshape(-1, 3) * (dt / (dt ** alpha))
    return dict(Theta_test_flat=Theta_test_flat, u_dot_test_clean_flat=u_dot_test_clean_flat,
                Phi_mat_train=Phi_mat_train, X_test_weak=X_test_weak, Y_test_weak=Y_test_weak)


def run_floor_cell(args):
    alpha, method, snr = args
    u_clean = solve_fractional_system(alpha, Nt, dt, Nx, _y0())
    pc = _precompute_oracle(u_clean, alpha)
    np.random.seed(int(snr) + 42)
    sel = np.empty(N_TRIALS)
    for t in range(N_TRIALS):
        u_noisy = add_noise(u_clean, snr)
        if method == "pointwise":
            sel[t] = select_temporal_order_pointwise_fast(
                u_noisy, u_clean, dt, dx, CANDIDATES, FRAC, THR, k_start, alpha,
                pc["Theta_test_flat"], pc["u_dot_test_clean_flat"])[1]["alpha_t"]
        else:
            sel[t] = select_temporal_order_weak_fast(
                u_noisy, u_clean, dt, dx, CANDIDATES, FRAC, THR, k_start, alpha,
                pc["Phi_mat_train"], pc["X_test_weak"], pc["Y_test_weak"])[1]["alpha_t"]
    hist = {f"{c:.1f}": int(np.sum(np.abs(sel - c) < 1e-5)) for c in CANDIDATES}
    return (alpha, method, int(snr), {
        "hist": hist,
        "success_rate": float(np.mean(np.abs(sel - alpha) < 1e-5)),
        "frac_high": float(np.mean(sel >= 0.8 - 1e-9)),
        "modal": max(hist, key=hist.get),
    })


def collect_wilson():
    """Wilson CIs for bracket-edge cells from all available result JSONs."""
    rows = []  # (system, method, snr, k, n, note)

    def add(system, method, snr, rate, n=500, note=""):
        rows.append((system, method, int(snr), int(round(rate * n)), n, note))

    # Primary oracle fine sweep + cp1 extended weak
    if os.path.exists("data/fine_snr_sweep_results.json"):
        fs = json.load(open("data/fine_snr_sweep_results.json"))
        for a in ["0.5", "0.7", "0.9"]:
            for m in ["pointwise", "weak"]:
                cells = fs[a][m]
                snrs = sorted(int(s) for s in cells)
                # the succeed edge and the fail edge bounding it
                succ = [s for s in snrs if cells[str(s)]["success_rate"] >= 0.95]
                if succ:
                    ls = min(succ)
                    add(f"primary oracle a={a}", m, ls, cells[str(ls)]["success_rate"], note="succeed edge")
                    below = [s for s in snrs if s < ls]
                    if below:
                        hf = max(below)
                        add(f"primary oracle a={a}", m, hf, cells[str(hf)]["success_rate"], note="fail edge")
    if os.path.exists("data/cp1_extended_weak.json"):
        ext = json.load(open("data/cp1_extended_weak.json"))
        for a in ["0.5", "0.7", "0.9"]:
            for s in ["15"]:
                if s in ext[a]:
                    add(f"primary oracle a={a}", "weak", int(s), ext[a][s]["success_rate"],
                        note="15 dB edge (flagged)")
    # CP1 realistic, CP2 VdP, CP3 fairgrid — bracket edges
    for path, sysname, struct in [
        ("data/cp1_realistic_selection.json", "realistic", "cp1"),
        ("data/cp2_benchmark_vdp.json", "VdP oracle", "cp2"),
        ("data/cp3_fairgrid.json", "fair grid", "cp3"),
    ]:
        if not os.path.exists(path):
            continue
        d = json.load(open(path))
        sweep = d.get("sweep", d)
        for a in ["0.5", "0.7", "0.9"]:
            methods = sweep[a].keys()
            for m in methods:
                cells = sweep[a][m]
                snrs = sorted(int(s) for s in cells)
                succ = [s for s in snrs if cells[str(s)]["success_rate"] >= 0.95]
                if succ:
                    ls = min(succ)
                    add(f"{sysname} a={a}", m, ls, cells[str(ls)]["success_rate"], note="succeed edge")
                    below = [s for s in snrs if s < ls]
                    if below:
                        hf = max(below)
                        add(f"{sysname} a={a}", m, hf, cells[str(hf)]["success_rate"], note="fail edge")
    return rows


def main():
    jobs = [(a, m, s) for a in TRUE_ALPHAS for m in ["pointwise", "weak"] for s in FLOOR_SNRS]
    print(f"CP4 noise-floor: {len(jobs)} cells x {N_TRIALS}")
    with Pool(processes=8) as pool:
        out = pool.map(run_floor_cell, jobs)
    floor = {str(a): {"pointwise": {}, "weak": {}} for a in TRUE_ALPHAS}
    for alpha, method, snr, res in out:
        floor[str(alpha)][method][str(snr)] = res
        print(f"  a={alpha} {method:9s} SNR={snr:>3} -> succ={res['success_rate']:.3f} "
              f"modal={res['modal']} frac>=0.8={res['frac_high']:.3f}")

    wrows = collect_wilson()
    payload = {"floor": floor, "wilson_rows": [
        {"system": s, "method": m, "snr": snr, "k": k, "n": n, "note": note,
         "p": wilson(k, n)[0], "ci_lo": wilson(k, n)[1], "ci_hi": wilson(k, n)[2]}
        for (s, m, snr, k, n, note) in wrows]}
    with open("data/cp4_noisefloor.json", "w") as f:
        json.dump(payload, f, indent=2)
    print("Saved data/cp4_noisefloor.json")
    write_md(floor, payload["wilson_rows"])
    print("Saved RESULTS_noise_floor_selection.md")


def write_md(floor, wrows):
    Ln = []
    A = Ln.append
    A("# Noise-Floor Selection Bias and Wilson Confidence Intervals (Checkpoint 4)\n")
    A("## 1. Order-selection distribution at the noise floor (oracle scoring)\n")
    A("Uniform chance over the 9 candidate orders is $1/9 \\approx 11.1\\%$. The tables "
      "below give, at near-pure-noise SNRs, the fraction of 500 realizations selecting "
      "each candidate order, the modal selection, and the fraction selecting a *high* "
      "order ($\\hat\\alpha \\ge 0.8$). A concentration of mass at high orders means a "
      "high *true* order (e.g. 0.9) inherits an inflated apparent success rate purely "
      "from the selector's high-order bias — not from genuine recovery.\n")
    for a in TRUE_ALPHAS:
        A(f"### True $\\alpha_t = {a:.1f}$\n")
        for m, ml in [("pointwise", "Pointwise GL"), ("weak", "Weak-form GL")]:
            A(f"#### {ml}:")
            A("| SNR (dB) | " + " | ".join(f"{c:.1f}" for c in CANDIDATES) +
              " | modal | $P(\\hat\\alpha\\ge0.8)$ | success |")
            A("| :---: | " + " | ".join(":--:" for _ in CANDIDATES) + " | :--: | :--: | :--: |")
            for snr in FLOOR_SNRS:
                r = floor[str(a)][m][str(snr)]
                cells = " | ".join(f"{r['hist'][f'{c:.1f}']/N_TRIALS:.2f}" for c in CANDIDATES)
                A(f"| {snr} | {cells} | {r['modal']} | {r['frac_high']:.2f} | {r['success_rate']:.3f} |")
            A("")
    A("## 2. Wilson 95% CIs for bracket-edge success rates\n")
    A("Each fail/succeed edge with its Wilson score interval. An edge is "
      "**criterion-fragile** if its 95% CI straddles the 0.95 line.\n")
    A("| System | Method | SNR | success $\\hat p$ | Wilson 95% CI | note | vs 0.95 |")
    A("| :--- | :--- | :---: | :---: | :---: | :--- | :---: |")
    for r in wrows:
        straddle = "fragile (straddles)" if r["ci_lo"] < 0.95 < r["ci_hi"] else \
                   ("robust <0.95" if r["ci_hi"] < 0.95 else "robust ≥0.95")
        A(f"| {r['system']} | {r['method']} | {r['snr']} | {r['p']:.3f} | "
          f"[{r['ci_lo']:.3f}, {r['ci_hi']:.3f}] | {r['note']} | {straddle} |")
    A("")
    with open("RESULTS_noise_floor_selection.md", "w") as f:
        f.write("\n".join(Ln))


if __name__ == "__main__":
    main()
