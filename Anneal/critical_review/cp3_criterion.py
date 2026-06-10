"""CHECKPOINT 3 (C2) — criterion-dependence of the lifetime trend and hazard shape.

C2 worry: the decreasing-median-with-N trend (and hence k(N)) might be an artifact of the
rho_std<0.04 death criterion timing a non-canonical channel. Test: re-detect death on the
SAME saved trajectories (beta=0.130, all five N) under literature-standard alternative
criteria, changing ONLY the label, and compare median tau(N) and Weibull k(N).

Criteria (all hold dt_hold=50, right-censored at T_max=12000), on the stored decimated traces:
  original         : rho_std(t) < 0.04                 (the paper's spatial-homogeneity rule)
  mean_coh_0.78    : rho_mean(t) > 0.78                 literature-standard structural rule —
                     the chimera dies when mean LOCAL coherence reaches the coherent plateau
                     (the incoherent region is lost; Wolfrum-Omel'chenko-style 'coherent region
                     fills the ring'). NB: in THIS system collapse is a RISE of mean coherence
                     (not a drop), so the fixed-level crossing is upward.
  struct_loss_0.08 : rho_std(t) < 0.08   } the chimera's two-region CONTRAST is lost (fires
  struct_loss_0.10 : rho_std(t) < 0.10   } earlier for runs that dissolve into a low-contrast
                     limbo before full homogenisation — the contamination-sensitivity probe.

Traces: N in {32,64,128} from results/traces/ (original campaign); N=192 regenerated
(critical_review/cp3_criterion/); N=256 from the CP2 re-run. All reproduce tau bit-for-bit.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RING = os.path.join(ROOT, "anneal-hazard")
sys.path.insert(0, RING)
from src.survival import fit_all  # noqa: E402
from src.ring_detector import detect_death_ring  # noqa: E402

OUT = os.path.join(HERE, "cp3_criterion")
os.makedirs(OUT, exist_ok=True)
DT_HOLD = 50.0
T_MAX = 12000.0
BETA = 0.130

TRACE_PATHS = {
    32: os.path.join(RING, "results", "traces", "cond_b0.130_N32.npz"),
    64: os.path.join(RING, "results", "traces", "cond_b0.130_N64.npz"),
    128: os.path.join(RING, "results", "traces", "cond_b0.130_N128.npz"),
    192: os.path.join(OUT, "cond_b0.130_N192.npz"),
    256: os.path.join(HERE, "cp2_validation", "cp2_traces_b0.13_N256.npz"),
}


def detect_rise(t, x, level, dt_hold, T_max):
    """First t0 such that x(t) > level for all t in [t0, t0+dt_hold]; else censored."""
    above = x > level
    n = len(t)
    i = 0
    while i < n:
        if not above[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and above[j + 1]:
            j += 1
        if t[j] - t[i] >= dt_hold:
            return float(t[i]), 1
        i = j + 1
    return float(T_max), 0


CRITERIA = {
    "original": lambda t, rs, rm: detect_death_ring(t, rs, 0.04, DT_HOLD, T_MAX),
    "mean_coh_0.78": lambda t, rs, rm: detect_rise(t, rm, 0.78, DT_HOLD, T_MAX),
    "struct_loss_0.08": lambda t, rs, rm: detect_death_ring(t, rs, 0.08, DT_HOLD, T_MAX),
    "struct_loss_0.10": lambda t, rs, rm: detect_death_ring(t, rs, 0.10, DT_HOLD, T_MAX),
}


def load_cell(N):
    d = np.load(TRACE_PATHS[N], allow_pickle=True)
    dec = int(d["decimate"]); dt = float(d["dt"])
    rho_std = d["rho_std"]; rho_mean = d["rho_mean"]
    runs = []
    for i in range(len(rho_std)):
        rs = np.asarray(rho_std[i], float); rm = np.asarray(rho_mean[i], float)
        t = np.arange(len(rs)) * dec * dt
        runs.append((t, rs, rm))
    return runs


def fit(taus, evs):
    fs = fit_all(np.asarray(taus, float), np.asarray(evs, int))
    klo, khi = fs.weibull.k_ci
    return {"k": float(fs.weibull.k), "k_lo": float(klo), "k_hi": float(khi),
            "ci_excl_1": bool(np.isfinite(klo) and klo > 1.0), "lrt_p": float(fs.lrt_p),
            "n": int(fs.n), "n_events": int(fs.n_events)}


def main():
    Ns = sorted(TRACE_PATHS)
    res = {c: {} for c in CRITERIA}
    cells = {N: load_cell(N) for N in Ns}
    print(f"beta={BETA}: loaded N={Ns} (runs each: " +
          ", ".join(f"{N}:{len(cells[N])}" for N in Ns) + ")")
    for cname, fn in CRITERIA.items():
        for N in Ns:
            taus, evs = [], []
            for t, rs, rm in cells[N]:
                tau, ev = fn(t, rs, rm)
                taus.append(tau); evs.append(ev)
            taus = np.array(taus); evs = np.array(evs)
            med = float(np.median(taus[evs == 1])) if evs.sum() else float("nan")
            f = fit(taus, evs)
            res[cname][N] = {"median_tau": med, "censor_frac": float(1 - evs.mean()), **f}

    # tables
    print(f"\n{'='*96}\nMEDIAN LIFETIME vs N  (beta={BETA}) — does the DECREASE persist under each criterion?\n{'='*96}")
    print(f"{'criterion':>18} | " + " ".join(f"N={N:>4}" for N in Ns) + " | trend")
    for c in CRITERIA:
        meds = [res[c][N]["median_tau"] for N in Ns]
        trend = "DECREASING" if all(meds[i] >= meds[i+1] for i in range(len(meds)-1)) else \
                ("decr(32->256)" if meds[0] > meds[-1] else "not-decr")
        print(f"{c:>18} | " + " ".join(f"{m:>6.0f}" for m in meds) + f" | {trend} "
              f"(ratio 256/32={meds[-1]/meds[0]:.2f})")

    print(f"\n{'='*96}\nWEIBULL SHAPE k vs N  (beta={BETA}) — is k>1 and the N-trend criterion-robust?\n{'='*96}")
    print(f"{'criterion':>18} | " + " ".join(f"N={N:>4}" for N in Ns) + " | all CI>1?")
    for c in CRITERIA:
        ks = [res[c][N]["k"] for N in Ns]
        all1 = all(res[c][N]["ci_excl_1"] for N in Ns)
        print(f"{c:>18} | " + " ".join(f"{k:>6.2f}" for k in ks) + f" |   {'YES' if all1 else 'no'}")
    print(f"\n(detail: k [CI], censoring per cell)")
    for c in CRITERIA:
        print(f"  {c}:")
        for N in Ns:
            r = res[c][N]
            print(f"     N={N:>4} k={r['k']:.3f} [{r['k_lo']:.3f},{r['k_hi']:.3f}] "
                  f"CI>1={r['ci_excl_1']} cens={100*r['censor_frac']:.1f}% "
                  f"medTau={r['median_tau']:.0f} LRTp={r['lrt_p']:.1e}")

    with open(os.path.join(OUT, "cp3_criterion.json"), "w") as f:
        json.dump({"beta": BETA, "Ns": Ns, "dt_hold": DT_HOLD, "results": res}, f, indent=2)
    print(f"\n[saved] {os.path.join(OUT, 'cp3_criterion.json')}")
    return res, Ns


if __name__ == "__main__":
    main()
