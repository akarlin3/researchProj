"""C1 claims 2,3,4 + refutation probes. Independent fitter = lifelines WeibullFitter."""
import json
import numpy as np
import pandas as pd
from lifelines import WeibullFitter
from scipy.special import gamma
import warnings
warnings.filterwarnings("ignore")

BASE = "/Users/averykarlin/annealMusic"
camp = [json.loads(l) for l in open(f"{BASE}/absorption_results/absorption_campaign.jsonl")]
camp = pd.DataFrame([d for d in camp if d["A"] == 0.5])
red = pd.DataFrame([json.loads(l) for l in open(f"{BASE}/reduced_results/reduced_runs.jsonl")])
df = camp.merge(red, on=["N", "seed"], suffixes=("", "_red"))
df["tau"] = df["t_abs"].astype(float)
df["event"] = (~df["abs_censored"].astype(bool)).astype(int)
Ns = [4, 8, 16, 24, 32, 48, 64]


def fitL(tau, event):
    wf = WeibullFitter().fit(np.asarray(tau, float), event_observed=np.asarray(event, int))
    s = wf.summary
    return float(wf.rho_), float(s.loc["rho_", "coef lower 95%"]), float(s.loc["rho_", "coef upper 95%"])


def cv_from_k(k):
    return np.sqrt(gamma(1 + 2/k) / gamma(1 + 1/k)**2 - 1)


def stratify_report(feature, nbins, label):
    print(f"\n{'='*70}\nStratify by {label} into {nbins} quantile bins\n{'='*70}")
    ks = []
    n_excl1 = 0
    n_cells = 0
    n_k_gt1 = 0
    failed = []
    for N in Ns:
        sub = df[df["N"] == N].copy()
        # quantile bins by feature within this N
        try:
            sub["bin"] = pd.qcut(sub[feature], nbins, labels=False, duplicates="drop")
        except Exception:
            sub["bin"] = pd.qcut(sub[feature].rank(method="first"), nbins, labels=False)
        for b in sorted(sub["bin"].dropna().unique()):
            cell = sub[sub["bin"] == b]
            if len(cell) < 5 or cell["event"].sum() < 3:
                continue
            n_cells += 1
            k, klo, khi = fitL(cell["tau"].values, cell["event"].values)
            ks.append(k)
            if k > 1:
                n_k_gt1 += 1
            if klo > 1:
                n_excl1 += 1
            else:
                failed.append((N, int(b), len(cell), round(k, 3), round(klo, 3), round(khi, 3)))
    print(f"cells={n_cells}  k>1: {n_k_gt1}/{n_cells}  CI-excludes-1: {n_excl1}/{n_cells}")
    if ks:
        print(f"k range: [{min(ks):.2f}, {max(ks):.2f}]")
    if failed:
        print(f"CELLS WITH CI INCLUDING 1 (k,CI): {failed}")
    return ks, n_cells, n_excl1, failed


# CLAIM 2: stratify by absdphi0, 4 bins
ks2, nc2, ne2, f2 = stratify_report("absdphi0", 4, "|Delta phi_0| (absdphi0)  [CLAIM 2]")
print("CLAIM 2: all 28 cells k>1, CI excl 1, k in [1.88,5.12]")

# CLAIM 3: stratify by t_capture, 4 bins  (DECISIVE)
ks3, nc3, ne3, f3 = stratify_report("t_capture", 4, "t_capture (reduced predicted lifetime)  [CLAIM 3]")
print("CLAIM 3: all 28 cells k>1, CI excl 1, k in [2.00,5.00]")

# CLAIM 4: within-bin CV of t_abs vs Weibull-k CV
print(f"\n{'='*70}\nCLAIM 4: within-bin CV(t_abs) vs sqrt(Gamma...) prediction\n{'='*70}")
rows = []
for N in Ns:
    sub = df[df["N"] == N].copy()
    sub["bin"] = pd.qcut(sub["absdphi0"], 4, labels=False, duplicates="drop")
    for b in sorted(sub["bin"].dropna().unique()):
        cell = sub[sub["bin"] == b]
        cv_emp = cell["t_abs"].std(ddof=1) / cell["t_abs"].mean()
        k, klo, khi = fitL(cell["tau"].values, cell["event"].values)
        rows.append((N, int(b), cv_emp, cv_from_k(k), k))
cvs = [r[2] for r in rows]
print(f"empirical within-bin CV(t_abs) range: [{min(cvs):.3f}, {max(cvs):.3f}]  (claim ~0.21-0.57)")
# correlation between empirical CV and Weibull-predicted CV
emp = np.array([r[2] for r in rows]); pred = np.array([r[3] for r in rows])
print(f"corr(empirical CV, Weibull-k CV) = {np.corrcoef(emp, pred)[0,1]:.3f}")
print(f"mean |emp - pred| = {np.mean(np.abs(emp-pred)):.4f}")
print("sample (N,bin,emp_CV,weibull_CV,k):")
for r in rows[:8]:
    print(f"  N={r[0]:3d} bin={r[1]} emp={r[2]:.3f} pred={r[3]:.3f} k={r[4]:.3f}")
