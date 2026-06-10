"""C1 REFUTATION PROBES: try to break k>1 by conditioning harder on collective IC."""
import json
import numpy as np
import pandas as pd
from lifelines import WeibullFitter
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


def probe(feature, nbins, focus_Ns=None):
    if focus_Ns is None:
        focus_Ns = Ns
    print(f"\n--- Probe: bin {feature} into {nbins} bins ---")
    allk = []
    n_incl1 = 0
    n_cells = 0
    for N in focus_Ns:
        sub = df[df["N"] == N].copy()
        try:
            sub["bin"] = pd.qcut(sub[feature], nbins, labels=False, duplicates="drop")
        except Exception:
            continue
        kline = []
        for b in sorted(sub["bin"].dropna().unique()):
            cell = sub[sub["bin"] == b]
            if len(cell) < 5 or cell["event"].sum() < 3:
                continue
            n_cells += 1
            k, klo, khi = fitL(cell["tau"].values, cell["event"].values)
            allk.append(k)
            flag = "  <<CI incl 1" if klo <= 1 else ""
            if klo <= 1:
                n_incl1 += 1
            kline.append(f"b{int(b)}(n={len(cell)}):k={k:.2f}[{klo:.2f},{khi:.2f}]{flag}")
        if N in [4, 8, 64]:  # show detail for a couple N
            print(f"  N={N:3d}: " + " ".join(kline))
    if allk:
        print(f"  => cells={n_cells}, k range [{min(allk):.2f},{max(allk):.2f}], "
              f"min k={min(allk):.2f}, CI-includes-1 in {n_incl1}/{n_cells} cells")
    return min(allk) if allk else None, n_incl1, n_cells


print("="*70)
print("PROBE A: finer t_capture bins (narrower => closer to fixed collective IC)")
print("Does within-stratum k collapse toward 1?")
print("="*70)
for nb in [4, 6, 8, 10]:
    probe("t_capture", nb)

print("\n" + "="*70)
print("PROBE B: alternative collective-IC predictors")
print("="*70)
probe("Rsync0", 4)
probe("Rincoh0", 4)
probe("dphi0", 4)

print("\n" + "="*70)
print("PROBE C: joint 2D bins (absdphi0 x Rincoh0), 3x3 per N -> ~9 cells/N")
print("="*70)
allk = []; n_incl1 = 0; n_cells = 0
for N in Ns:
    sub = df[df["N"] == N].copy()
    sub["b1"] = pd.qcut(sub["absdphi0"], 3, labels=False, duplicates="drop")
    sub["b2"] = pd.qcut(sub["Rincoh0"], 3, labels=False, duplicates="drop")
    for (b1, b2), cell in sub.groupby(["b1", "b2"]):
        if len(cell) < 8 or cell["event"].sum() < 4:
            continue
        n_cells += 1
        k, klo, khi = fitL(cell["tau"].values, cell["event"].values)
        allk.append(k)
        if klo <= 1:
            n_incl1 += 1
            print(f"  N={N} cell=({b1},{b2}) n={len(cell)} k={k:.2f} CI=[{klo:.2f},{khi:.2f}] <<incl 1")
print(f"=> 2D joint: cells={n_cells}, k range [{min(allk):.2f},{max(allk):.2f}], "
      f"min k={min(allk):.2f}, CI-incl-1 in {n_incl1}/{n_cells}")

print("\n" + "="*70)
print("PROBE D: 3D joint (absdphi0 x Rincoh0 x dphi0) 2x2x2 per N, narrowest IC volume")
print("="*70)
allk = []; n_incl1 = 0; n_cells = 0; mins = []
for N in Ns:
    sub = df[df["N"] == N].copy()
    sub["b1"] = pd.qcut(sub["absdphi0"], 2, labels=False, duplicates="drop")
    sub["b2"] = pd.qcut(sub["Rincoh0"], 2, labels=False, duplicates="drop")
    sub["b3"] = pd.qcut(sub["dphi0"], 2, labels=False, duplicates="drop")
    for keys, cell in sub.groupby(["b1", "b2", "b3"]):
        if len(cell) < 8 or cell["event"].sum() < 4:
            continue
        n_cells += 1
        k, klo, khi = fitL(cell["tau"].values, cell["event"].values)
        allk.append(k)
        if klo <= 1:
            n_incl1 += 1
print(f"=> 3D joint: cells={n_cells}, k range [{min(allk):.2f},{max(allk):.2f}], "
      f"min k={min(allk):.2f}, CI-incl-1 in {n_incl1}/{n_cells}")

print("\n" + "="*70)
print("PROBE E: pool all N, finest single-predictor t_capture bins (max stats per narrow IC)")
print("="*70)
for nb in [8, 12, 20]:
    sub = df.copy()
    sub["bin"] = pd.qcut(sub["t_capture"], nb, labels=False, duplicates="drop")
    ks = []; incl = 0; nc = 0
    for b in sorted(sub["bin"].dropna().unique()):
        cell = sub[sub["bin"] == b]
        if len(cell) < 10:
            continue
        nc += 1
        k, klo, khi = fitL(cell["tau"].values, cell["event"].values)
        ks.append(k)
        if klo <= 1:
            incl += 1
    print(f"  {nb} pooled bins: cells={nc}, k range [{min(ks):.2f},{max(ks):.2f}], "
          f"min k={min(ks):.2f}, CI-incl-1 in {incl}/{nc}")
