"""INDEPENDENT adversarial verification of CHECKPOINT 3 (C2 criterion-dependence).

Re-loads raw decimated traces, re-detects death with an INDEPENDENT detector,
and fits a right-censored Weibull with a HAND-ROLLED scipy MLE (direct (logk,loglam)
joint L-BFGS-B optimization of the censored log-likelihood + profile-likelihood CI),
which is a DIFFERENT method than the author's profiled Nelder-Mead fitter.

Survival: S(t)=exp(-(t/lam)^k). event=1 death observed, 0 right-censored at T_max.
"""
import os
import numpy as np
from scipy import optimize, stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)  # worktree root
RING = os.path.join(ROOT, "anneal-hazard")

DT_HOLD = 50.0
T_MAX = 12000.0

TRACE_PATHS = {
    32: os.path.join(RING, "results", "traces", "cond_b0.130_N32.npz"),
    64: os.path.join(RING, "results", "traces", "cond_b0.130_N64.npz"),
    128: os.path.join(RING, "results", "traces", "cond_b0.130_N128.npz"),
    192: os.path.join(HERE, "cp3_criterion", "cond_b0.130_N192.npz"),
    256: os.path.join(HERE, "cp2_validation", "cp2_traces_b0.13_N256.npz"),
}
NS = [32, 64, 128, 192, 256]


# ----------------------------------------------------------- independent detectors
def detect_fall(t, x, level, dt_hold=DT_HOLD, T_max=T_MAX):
    """First t0 such that x(t) < level continuously for >= dt_hold. Independent
    re-implementation using run-length grouping of the boolean mask."""
    below = np.asarray(x) < level
    if not below.any():
        return float(T_max), 0
    # find maximal runs of True
    idx = np.flatnonzero(below)
    # split where index gaps appear
    splits = np.flatnonzero(np.diff(idx) > 1)
    starts = np.concatenate(([idx[0]], idx[splits + 1]))
    ends = np.concatenate((idx[splits], [idx[-1]]))
    for s, e in zip(starts, ends):
        if t[e] - t[s] >= dt_hold:
            return float(t[s]), 1
    return float(T_max), 0


def detect_rise(t, x, level, dt_hold=DT_HOLD, T_max=T_MAX):
    """First t0 such that x(t) > level continuously for >= dt_hold (upward crossing)."""
    above = np.asarray(x) > level
    if not above.any():
        return float(T_max), 0
    idx = np.flatnonzero(above)
    splits = np.flatnonzero(np.diff(idx) > 1)
    starts = np.concatenate(([idx[0]], idx[splits + 1]))
    ends = np.concatenate((idx[splits], [idx[-1]]))
    for s, e in zip(starts, ends):
        if t[e] - t[s] >= dt_hold:
            return float(t[s]), 1
    return float(T_max), 0


# ------------------------------------------------ hand-rolled censored Weibull MLE
def weibull_negll(params, tau, event):
    """negative log-lik for right-censored Weibull, params=(logk, loglam).
    contribution: death -> log f = log k - log lam + (k-1)(log t - log lam) - (t/lam)^k
                  censored -> log S = -(t/lam)^k
    """
    logk, loglam = params
    k = np.exp(logk); lam = np.exp(loglam)
    z = tau / lam
    zk = z ** k
    ll = np.where(event == 1,
                  np.log(k) - np.log(lam) + (k - 1.0) * np.log(z) - zk,
                  -zk)
    return -np.sum(ll)


def fit_weibull_indep(tau, event):
    tau = np.asarray(tau, float); event = np.asarray(event, int)
    d = int(event.sum())
    if d == 0:
        return dict(k=np.nan, lam=np.nan, k_lo=np.nan, k_hi=np.nan, ci_excl_1=False, llmax=np.nan)
    # init: method-of-moments-ish
    mt = np.mean(tau[event == 1])
    x0 = [np.log(1.5), np.log(mt)]
    res = optimize.minimize(weibull_negll, x0, args=(tau, event),
                            method="L-BFGS-B",
                            options={"ftol": 1e-12, "gtol": 1e-9, "maxiter": 5000})
    logk, loglam = res.x
    k_hat = float(np.exp(logk)); lam_hat = float(np.exp(loglam))
    llmax = -res.fun

    # profile-likelihood 95% CI on k: minimize over lam at fixed k
    thr = stats.chi2.ppf(0.95, 1) / 2.0

    def prof_negll_at_k(k):
        # profile out lam analytically: lam^k = sum(tau^k)/d
        tk = tau ** k
        lam_k = tk.sum() / d
        lam = lam_k ** (1.0 / k)
        return weibull_negll([np.log(k), np.log(lam)], tau, event)

    llmax_prof = -prof_negll_at_k(k_hat)

    def gap(k):
        return (llmax_prof - (-prof_negll_at_k(k))) - thr

    klo = khi = np.nan
    try:
        a = k_hat
        while gap(a) < 0 and a > 1e-3:
            a *= 0.5
        if gap(a) > 0:
            klo = optimize.brentq(gap, a, k_hat, xtol=1e-6)
    except Exception:
        pass
    try:
        b = k_hat
        while gap(b) < 0 and b < 1e3:
            b *= 2.0
        if gap(b) > 0:
            khi = optimize.brentq(gap, k_hat, b, xtol=1e-6)
    except Exception:
        pass
    ci_excl_1 = bool(np.isfinite(klo) and klo > 1.0)
    return dict(k=k_hat, lam=lam_hat, k_lo=float(klo), k_hi=float(khi),
                ci_excl_1=ci_excl_1, llmax=float(llmax))


# --------------------------------------------------------- Kaplan-Meier (independent)
def km_lnS_slope(tau, event):
    """KM estimator + log-cumulative-hazard slope (early vs late) as a non-Weibull
    increasing-hazard check. Returns (n_events, ln(-lnS) regression slope, KM lnS curvature)."""
    tau = np.asarray(tau, float); event = np.asarray(event, int)
    order = np.argsort(tau)
    tau, event = tau[order], event[order]
    n = len(tau)
    uniq = np.unique(tau[event == 1])
    S = 1.0
    times, surv = [], []
    for t in uniq:
        at_risk = np.sum(tau >= t)
        d = np.sum((tau == t) & (event == 1))
        if at_risk == 0:
            continue
        S *= (1 - d / at_risk)
        times.append(t); surv.append(S)
    times = np.array(times); surv = np.array(surv)
    # log-cumulative-hazard: ln(-ln S) vs ln t -> slope = Weibull k if Weibull holds
    m = (surv > 0) & (surv < 1) & (times > 0)
    if m.sum() < 3:
        return dict(n_events=int(event.sum()), llh_slope=np.nan)
    x = np.log(times[m]); y = np.log(-np.log(surv[m]))
    slope, intercept = np.polyfit(x, y, 1)
    return dict(n_events=int(event.sum()), llh_slope=float(slope),
                times=times, surv=surv)


def load_cell(N):
    d = np.load(TRACE_PATHS[N], allow_pickle=True)
    dec = int(d["decimate"]); dt = float(d["dt"])
    rs_all = d["rho_std"]; rm_all = d["rho_mean"]
    runs = []
    for i in range(len(rs_all)):
        rs = np.asarray(rs_all[i], float); rm = np.asarray(rm_all[i], float)
        t = np.arange(len(rs)) * dec * dt
        runs.append((t, rs, rm))
    return runs


def run_criterion(cells, detector):
    out = {}
    for N in NS:
        taus, evs, finals_rm, finals_rs = [], [], [], []
        for t, rs, rm in cells[N]:
            tau, ev = detector(t, rs, rm)
            taus.append(tau); evs.append(ev)
            finals_rm.append(rm[-1]); finals_rs.append(rs[-1])
        taus = np.array(taus); evs = np.array(evs)
        med = float(np.median(taus[evs == 1])) if evs.sum() else float("nan")
        fit = fit_weibull_indep(taus, evs)
        km = km_lnS_slope(taus, evs)
        out[N] = dict(median=med, censor_frac=float(1 - evs.mean()),
                      n_events=int(evs.sum()), taus=taus, evs=evs,
                      finals_rm=np.array(finals_rm), finals_rs=np.array(finals_rs),
                      llh_slope=km["llh_slope"], **fit)
    return out


def main():
    cells = {N: load_cell(N) for N in NS}
    print("loaded:", {N: len(cells[N]) for N in NS})

    criteria = {
        "original (rho_std<0.04)": lambda t, rs, rm: detect_fall(t, rs, 0.04),
        "struct_loss_0.08":        lambda t, rs, rm: detect_fall(t, rs, 0.08),
        "struct_loss_0.10":        lambda t, rs, rm: detect_fall(t, rs, 0.10),
        "mean_coh_0.78 (rise)":    lambda t, rs, rm: detect_rise(t, rm, 0.78),
    }

    results = {c: run_criterion(cells, fn) for c, fn in criteria.items()}

    print("\n" + "=" * 90)
    print("MEDIAN tau(N)  [independent detect]")
    print("=" * 90)
    print(f"{'criterion':>26} | " + " ".join(f"N={N:>4}" for N in NS) + " | decreasing?")
    for c, r in results.items():
        meds = [r[N]["median"] for N in NS]
        decr = all(meds[i] >= meds[i + 1] for i in range(len(meds) - 1))
        decr_overall = meds[0] > meds[-1]
        tag = "MONO-DECR" if decr else ("decr(32>256)" if decr_overall else "NOT-DECR")
        print(f"{c:>26} | " + " ".join(f"{m:>6.0f}" for m in meds) + f" | {tag}")

    print("\n" + "=" * 90)
    print("WEIBULL k(N)  [hand-rolled MLE]  (k [CI_lo,CI_hi]  CI>1?  cens%)")
    print("=" * 90)
    for c, r in results.items():
        print(f"  {c}:")
        ks = []
        for N in NS:
            x = r[N]
            ks.append(x["k"])
            print(f"     N={N:>4}  k={x['k']:.3f} [{x['k_lo']:.3f},{x['k_hi']:.3f}] "
                  f"CI>1={str(x['ci_excl_1']):>5}  cens={100*x['censor_frac']:.1f}%  "
                  f"med={x['median']:.0f}  KM_llh_slope={x['llh_slope']:.2f}")
        print(f"     k(N) = [{', '.join(f'{k:.3f}' for k in ks)}]  "
              f"all CI>1: {all(r[N]['ci_excl_1'] for N in NS)}")

    # --------- ARTIFACT PROBE: mean_coh_0.78 censored runs at N=256 ---------
    print("\n" + "=" * 90)
    print("ARTIFACT PROBE: mean_coh>0.78 censored runs (final rho_mean & rho_std<0.04 death?)")
    print("=" * 90)
    mc = results["mean_coh_0.78 (rise)"]
    for N in NS:
        x = mc[N]
        cens_mask = (x["evs"] == 0)
        ncens = int(cens_mask.sum())
        frm = x["finals_rm"][cens_mask]
        # of the censored-by-meancoh runs, did rho_std<0.04 fire? re-detect
        runs = cells[N]
        std_dead = 0
        std_final = []
        rm_final_list = []
        for j, (t, rs, rm) in enumerate(runs):
            if not cens_mask[j]:
                continue
            tau_s, ev_s = detect_fall(t, rs, 0.04)
            std_dead += ev_s
            std_final.append(rs[-1])
            rm_final_list.append(rm[-1])
        rm_final_list = np.array(rm_final_list)
        if ncens:
            print(f"  N={N:>4}: meancoh-censored={ncens} ({100*ncens/300:.1f}%)  "
                  f"final rho_mean: med={np.median(rm_final_list):.3f} "
                  f"[{rm_final_list.min():.3f},{rm_final_list.max():.3f}]  "
                  f"of these rho_std<0.04 DEAD={std_dead}/{ncens}")
        else:
            print(f"  N={N:>4}: meancoh-censored=0")

    # excluding twisted (final rho_mean<0.7) collapses, re-fit mean_coh
    print("\n" + "=" * 90)
    print("mean_coh>0.78 EXCLUDING twisted collapses (drop runs whose final rho_mean<0.70)")
    print("=" * 90)
    ks_excl = []
    for N in NS:
        x = mc[N]
        runs = cells[N]
        keep_tau, keep_ev = [], []
        dropped = 0
        for j, (t, rs, rm) in enumerate(runs):
            if rm[-1] < 0.70:   # twisted / non-sync final state
                dropped += 1
                continue
            keep_tau.append(x["taus"][j]); keep_ev.append(x["evs"][j])
        fit = fit_weibull_indep(np.array(keep_tau), np.array(keep_ev))
        ks_excl.append(fit["k"])
        nev = int(np.sum(np.array(keep_ev) == 1))
        med = float(np.median(np.array(keep_tau)[np.array(keep_ev) == 1])) if nev else np.nan
        print(f"  N={N:>4}: dropped {dropped} twisted; kept {len(keep_tau)}; "
              f"k={fit['k']:.3f} [{fit['k_lo']:.3f},{fit['k_hi']:.3f}] CI>1={fit['ci_excl_1']} "
              f"cens={100*(1-np.mean(keep_ev)):.1f}% med={med:.0f}")
    print(f"  k(N) excl-twisted = [{', '.join(f'{k:.3f}' for k in ks_excl)}]")

    return results


if __name__ == "__main__":
    main()
