"""Survival & hazard analysis (self-implemented, auditable).

All estimators take right-censored data as two arrays:
    tau   : observed time (death time, or T_max if censored)
    event : 1 = death observed, 0 = right-censored

Provides: Kaplan-Meier with Greenwood CI, Nelson-Aalen, Epanechnikov kernel hazard,
censored MLE for exponential and Weibull (profile-likelihood CI on the Weibull shape k),
the Weibull-vs-exponential LRT, and the ln S(t) linearity / runs-test diagnostics.

    Weibull:  S(t) = exp(-(t/lambda)^k),   exponential <=> k = 1.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import optimize, stats


# --------------------------------------------------------------------------- KM
@dataclass
class KMResult:
    times: np.ndarray      # distinct event times
    surv: np.ndarray       # S(t) right after each event time
    var: np.ndarray        # Greenwood variance of S(t)
    ci_lo: np.ndarray
    ci_hi: np.ndarray
    n_at_risk: np.ndarray
    n_events: np.ndarray


def kaplan_meier(tau, event, alpha=0.05) -> KMResult:
    tau = np.asarray(tau, float); event = np.asarray(event, int)
    order = np.argsort(tau, kind="mergesort")
    tau, event = tau[order], event[order]
    n = len(tau)
    uniq = np.unique(tau[event == 1])
    times, surv, var, n_risk, n_ev = [], [], [], [], []
    S = 1.0
    cum = 0.0  # Greenwood running sum
    for t in uniq:
        at_risk = int(np.sum(tau >= t))
        d = int(np.sum((tau == t) & (event == 1)))
        if at_risk == 0:
            continue
        S *= (1.0 - d / at_risk)
        if at_risk > d:
            cum += d / (at_risk * (at_risk - d))
        v = S * S * cum
        times.append(t); surv.append(S); var.append(v)
        n_risk.append(at_risk); n_ev.append(d)
    times = np.array(times); surv = np.array(surv); var = np.array(var)
    z = stats.norm.ppf(1 - alpha / 2)
    se = np.sqrt(var)
    ci_lo = np.clip(surv - z * se, 0, 1)
    ci_hi = np.clip(surv + z * se, 0, 1)
    return KMResult(times, surv, var, ci_lo, ci_hi, np.array(n_risk), np.array(n_ev))


# -------------------------------------------------------------------- Nelson-Aalen
@dataclass
class NAResult:
    times: np.ndarray
    cumhaz: np.ndarray
    var: np.ndarray
    incr: np.ndarray       # d_i / n_i increments (used by the kernel hazard)
    n_at_risk: np.ndarray
    n_events: np.ndarray


def nelson_aalen(tau, event) -> NAResult:
    tau = np.asarray(tau, float); event = np.asarray(event, int)
    uniq = np.unique(tau[event == 1])
    times, H, V, incr, n_risk, n_ev = [], [], [], [], [], []
    h = 0.0; v = 0.0
    for t in uniq:
        at_risk = int(np.sum(tau >= t))
        d = int(np.sum((tau == t) & (event == 1)))
        if at_risk == 0:
            continue
        h += d / at_risk
        v += d / (at_risk * at_risk)
        times.append(t); H.append(h); V.append(v); incr.append(d / at_risk)
        n_risk.append(at_risk); n_ev.append(d)
    return NAResult(np.array(times), np.array(H), np.array(V), np.array(incr),
                    np.array(n_risk), np.array(n_ev))


def epanechnikov_hazard(na: NAResult, grid: np.ndarray, bandwidth: float) -> np.ndarray:
    """Kernel-smoothed hazard h(t) = sum_i K_b(t - t_i) * dH(t_i), Epanechnikov kernel."""
    h = np.zeros_like(grid, dtype=float)
    b = bandwidth
    for ti, di in zip(na.times, na.incr):
        u = (grid - ti) / b
        k = np.where(np.abs(u) <= 1, 0.75 * (1 - u * u) / b, 0.0)
        h += k * di
    return h


# ------------------------------------------------------------- censored MLE fits
def _loglik_exp(tau, event):
    """Return (lambda_hat, loglik_max) for exponential S=exp(-t/lambda)."""
    tau = np.asarray(tau, float); event = np.asarray(event, int)
    d = int(event.sum()); tot = float(tau.sum())
    if d == 0:
        return float("inf"), float("nan")
    lam = tot / d
    ll = -d * np.log(lam) - tot / lam
    return lam, ll


def _loglik_weibull_given_k(k, tau, event, d, sum_log_t_death):
    if not (1e-3 < k < 1e2):   # keep tau**k away from under/overflow
        return -np.inf
    with np.errstate(over="ignore", invalid="ignore"):
        tk = tau ** k
    Stk = tk.sum()
    if d == 0 or not np.isfinite(Stk) or Stk <= 0:
        return -np.inf
    # profiled lambda: lam^k = Stk / d
    lam_k = Stk / d
    # ll = d ln k - d ln(lam^k)/... expand with lam:
    # ll = d ln k + (k-1) sum_death ln t - d k ln lam - sum_all (t/lam)^k
    #    = d ln k + (k-1) S_lnt - d ln(lam^k) - Stk/lam^k
    ll = d * np.log(k) + (k - 1.0) * sum_log_t_death - d * np.log(lam_k) - Stk / lam_k
    return ll


@dataclass
class WeibullFit:
    k: float
    lam: float
    loglik: float
    k_ci: tuple
    lam_ci: tuple = (float("nan"), float("nan"))


def fit_weibull(tau, event) -> WeibullFit:
    tau = np.asarray(tau, float); event = np.asarray(event, int)
    d = int(event.sum())
    sum_log_t_death = float(np.log(tau[event == 1]).sum()) if d else 0.0

    def neg(logk):
        k = np.exp(logk[0])
        return -_loglik_weibull_given_k(k, tau, event, d, sum_log_t_death)

    res = optimize.minimize(neg, x0=[0.0], method="Nelder-Mead",
                            options={"xatol": 1e-6, "fatol": 1e-9})
    k_hat = float(np.exp(res.x[0]))
    lam_hat = float((np.sum(tau ** k_hat) / d) ** (1.0 / k_hat))
    ll_max = _loglik_weibull_given_k(k_hat, tau, event, d, sum_log_t_death)

    # profile-likelihood 95% CI on k: {k : 2(ll_max - ll(k)) <= chi2_{1,0.95}}
    thr = stats.chi2.ppf(0.95, 1) / 2.0

    def gap(k):
        return (ll_max - _loglik_weibull_given_k(k, tau, event, d, sum_log_t_death)) - thr

    klo = khi = float("nan")
    try:
        lo_lo = k_hat
        while gap(lo_lo) < 0 and lo_lo > 1e-3:
            lo_lo *= 0.5
        if gap(lo_lo) > 0:
            klo = optimize.brentq(gap, lo_lo, k_hat)
    except Exception:
        pass
    try:
        hi_hi = k_hat
        while gap(hi_hi) < 0 and hi_hi < 1e3:
            hi_hi *= 2.0
        if gap(hi_hi) > 0:
            khi = optimize.brentq(gap, k_hat, hi_hi)
    except Exception:
        pass

    return WeibullFit(k_hat, lam_hat, float(ll_max), (klo, khi))


@dataclass
class FitSummary:
    n: int
    n_events: int
    exp_lambda: float
    exp_loglik: float
    weibull: WeibullFit
    lrt_stat: float
    lrt_p: float


def fit_all(tau, event) -> FitSummary:
    tau = np.asarray(tau, float); event = np.asarray(event, int)
    lam, ll_e = _loglik_exp(tau, event)
    wf = fit_weibull(tau, event)
    stat = 2.0 * (wf.loglik - ll_e)
    stat = max(stat, 0.0)
    p = float(stats.chi2.sf(stat, 1))
    return FitSummary(len(tau), int(event.sum()), lam, ll_e, wf, stat, p)


# ---------------------------------------------------- ln S(t) linearity + runs test
def lnS_linearity(km: KMResult):
    """Regress ln S(t) on t over event times with S>0; return R^2, slope, runs-test p
    on the sign of residuals (Wald-Wolfowitz). Linear ln S <=> constant hazard."""
    m = km.surv > 0
    t = km.times[m]; lnS = np.log(km.surv[m])
    if len(t) < 3:
        return {"r2": float("nan"), "slope": float("nan"), "runs_p": float("nan"), "n": int(len(t))}
    slope, intercept = np.polyfit(t, lnS, 1)
    pred = slope * t + intercept
    ss_res = float(np.sum((lnS - pred) ** 2))
    ss_tot = float(np.sum((lnS - lnS.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    runs_p = _runs_test(np.sign(lnS - pred))
    return {"r2": r2, "slope": float(slope), "runs_p": runs_p, "n": int(len(t))}


def _runs_test(signs) -> float:
    """Wald-Wolfowitz runs test on a sequence of +/- signs; returns two-sided p."""
    s = signs[signs != 0]
    n = len(s)
    if n < 2:
        return float("nan")
    n_pos = int(np.sum(s > 0)); n_neg = int(np.sum(s < 0))
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    runs = 1 + int(np.sum(s[1:] != s[:-1]))
    mu = 1 + 2 * n_pos * n_neg / n
    var = (2 * n_pos * n_neg * (2 * n_pos * n_neg - n)) / (n * n * (n - 1))
    if var <= 0:
        return float("nan")
    z = (runs - mu) / np.sqrt(var)
    return float(2 * stats.norm.sf(abs(z)))


# --------------------------------------------------------------------- self-test
if __name__ == "__main__":
    rng = np.random.default_rng(0)
    Tcap = 5.0

    # exponential truth, k should be ~1
    x = rng.exponential(1.0, size=3000)
    ev = (x <= Tcap).astype(int); xt = np.minimum(x, Tcap)
    fs = fit_all(xt, ev)
    print(f"[exp truth] events={fs.n_events} exp_lambda={fs.exp_lambda:.3f} "
          f"weibull k={fs.weibull.k:.3f} CI={tuple(round(v,3) for v in fs.weibull.k_ci)} "
          f"LRT p={fs.lrt_p:.3f}")

    # weibull truth k=1.8
    k_true = 1.8
    x = rng.weibull(k_true, size=3000)
    ev = (x <= Tcap).astype(int); xt = np.minimum(x, Tcap)
    fs = fit_all(xt, ev)
    km = kaplan_meier(xt, ev)
    lin = lnS_linearity(km)
    print(f"[weibull k=1.8] events={fs.n_events} weibull k={fs.weibull.k:.3f} "
          f"CI={tuple(round(v,3) for v in fs.weibull.k_ci)} LRT p={fs.lrt_p:.2e} "
          f"| lnS R2={lin['r2']:.4f} runs_p={lin['runs_p']:.3f}")
