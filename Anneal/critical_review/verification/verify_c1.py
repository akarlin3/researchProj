"""Independent adversarial verification of CHECKPOINT 1 (A=0.5).

Independent fitter: lifelines.WeibullFitter (autograd MLE, Wald CIs) --
a DIFFERENT implementation/CI method than the author's profile-likelihood scipy fitter.
We also cross-check with a hand-rolled scipy censored-Weibull MLE.

Convention: S(t)=exp(-(t/lambda)^k). event=1 = death observed (NOT abs_censored).
"""
import json
import numpy as np
import pandas as pd
from lifelines import WeibullFitter
from scipy import optimize, stats
import warnings
warnings.filterwarnings("ignore")

BASE = "/Users/averykarlin/annealMusic"

# ----------------------------------------------------------------- load + join
camp = []
with open(f"{BASE}/absorption_results/absorption_campaign.jsonl") as f:
    for line in f:
        d = json.loads(line)
        if d["A"] == 0.5:
            camp.append(d)
camp = pd.DataFrame(camp)

red = []
with open(f"{BASE}/reduced_results/reduced_runs.jsonl") as f:
    for line in f:
        red.append(json.loads(line))
red = pd.DataFrame(red)

print(f"campaign A=0.5 rows: {len(camp)}, reduced rows: {len(red)}")

df = camp.merge(red, on=["N", "seed"], suffixes=("", "_red"))
print(f"joined rows: {len(df)}")
# Join correctness: t_abs (campaign) == t_abs_meas (reduced)
mism = (np.abs(df["t_abs"] - df["t_abs_meas"]) > 1e-9).sum()
print(f"t_abs vs t_abs_meas mismatches: {mism}")
cens_mism = (df["abs_censored"].astype(bool) != df["abs_censored_meas"].astype(bool)).sum()
print(f"abs_censored vs abs_censored_meas mismatches: {cens_mism}")

df["tau"] = df["t_abs"].astype(float)
df["event"] = (~df["abs_censored"].astype(bool)).astype(int)
TMAX = 2000.0
print(f"overall censored frac: {1 - df['event'].mean():.4f}")
print(f"any tau==tmax among events? {((df['tau']>=TMAX)&(df['event']==1)).sum()}")


# --------------------------------------------- independent hand-rolled scipy MLE
def hand_weibull(tau, event):
    """Right-censored Weibull MLE from scratch (log-param), profile-likelihood CI on k.
    Completely separate code path from both author and lifelines."""
    tau = np.asarray(tau, float); event = np.asarray(event, int)
    d = int(event.sum())
    lt = np.log(tau)
    sld = lt[event == 1].sum()

    def negll_k(k):
        # profiled lambda
        tk = tau ** k
        lam_k = tk.sum() / d
        ll = d*np.log(k) + (k-1)*sld - d*np.log(lam_k) - tk.sum()/lam_k
        return -ll

    r = optimize.minimize_scalar(negll_k, bounds=(0.05, 50), method="bounded",
                                 options={"xatol": 1e-8})
    khat = float(r.x)
    lam = float((np.sum(tau**khat)/d)**(1/khat))
    llmax = -negll_k(khat)
    thr = stats.chi2.ppf(0.95, 1)/2
    def gap(k): return (llmax - (-negll_k(k))) - thr
    klo = optimize.brentq(gap, 1e-3, khat) if gap(1e-3) > 0 else np.nan
    khi = optimize.brentq(gap, khat, 200) if gap(200) > 0 else np.nan
    return khat, lam, (klo, khi)


def fit_lifelines(tau, event):
    """lifelines WeibullFitter. NOTE lifelines uses S=exp(-(t/lambda_)^rho_),
    so rho_ == k (shape), lambda_ == scale. Wald CI on rho_."""
    wf = WeibullFitter()
    wf.fit(tau, event_observed=event)
    k = wf.rho_
    lam = wf.lambda_
    ci = wf.confidence_interval_
    # row label for rho_
    klo = wf.summary.loc["rho_", "coef lower 95%"]
    khi = wf.summary.loc["rho_", "coef upper 95%"]
    return float(k), float(lam), (float(klo), float(khi))


def cv_from_k(k):
    from scipy.special import gamma
    return np.sqrt(gamma(1+2/k)/gamma(1+1/k)**2 - 1)


Ns = [4, 8, 16, 24, 32, 48, 64]

print("\n" + "="*70)
print("CLAIM 1: pooled per-N Weibull shape k_abs")
print("="*70)
claim1 = [2.27, 2.22, 2.49, 2.67, 2.44, 2.72, 3.03]
ll_k, hand_k = [], []
for N, ck in zip(Ns, claim1):
    sub = df[df["N"] == N]
    kL, lamL, ciL = fit_lifelines(sub["tau"].values, sub["event"].values)
    kH, lamH, ciH = hand_weibull(sub["tau"].values, sub["event"].values)
    ll_k.append(kL); hand_k.append(kH)
    excl_L = ciL[0] > 1
    excl_H = (ciH[0] > 1)
    print(f"N={N:3d} n={len(sub)} ev={sub['event'].sum():3d} | "
          f"lifelines k={kL:.3f} CI=({ciL[0]:.3f},{ciL[1]:.3f}) excl1={excl_L} | "
          f"hand k={kH:.3f} CI=({ciH[0]:.3f},{ciH[1]:.3f}) excl1={excl_H} | claim={ck}")
print("lifelines k rounded:", [round(x,2) for x in ll_k])
print("hand      k rounded:", [round(x,2) for x in hand_k])
print("claimed:            ", claim1)
