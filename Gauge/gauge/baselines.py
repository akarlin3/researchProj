"""Model-based IVIM uncertainty-quantification baselines (Fashion-independent).

These are the *standard* model-based UQ family from the IVIM / probabilistic-ML
literature -- the methods Gauge's distribution-free conformal approach is measured
against. Fashion's normalizing flow is NOT here; it is related work only.

Baselines
---------
* ``SingleGaussianPNN``  -- one probabilistic NN with a (mean, log-var) head
  trained by Gaussian NLL (the simple heteroscedastic baseline).
* ``MDNDeepEnsemble``    -- a deep ensemble of Mixture Density Networks (the
  Casali approach: total predictive uncertainty with an aleatoric/epistemic
  split).
* ``DeepEnsemblePoint``  -- a plain deep ensemble of point regressors (epistemic
  spread only; the deliberately weak baseline).
* ``BayesianIVIM_MCMC``  -- per-voxel Bayesian fit of the bi-exponential model by
  vectorized random-walk Metropolis (posterior intervals; the gold-standard
  Bayesian reference).

Every baseline exposes ``predict_samples(X) -> (N, 3, S)`` predictive draws in
*physical* units for (D, D*, f). Intervals at any level come from
``quantiles_from_samples``; this keeps raw / conformal / conformalized-model-based
comparisons on identical footing in benchmark.py.
"""
import os
import pickle
import time

import numpy as np

from gauge.forward import ivim_signal, DEFAULT_B_VALUES
from gauge.cohort import (generate_cohort, D_RANGE, DSTAR_RANGE, F_RANGE,
                          DEFAULT_SNR_GRID, DEFAULT_SEED)
from gauge.estimators import fit_nlls, IVIMQuantileRegressor

PARAM_NAMES = ("D", "D*", "f")
N_SAMPLES = 500          # predictive draws per point
ALPHAS = (0.05, 0.10, 0.20, 0.30)
_RESULTS_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
_CACHE = os.path.join(_RESULTS_DIR, "predictions.pkl")


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
class _Std:
    """Per-column standardizer (fit on train, applied everywhere)."""

    def fit(self, A):
        A = np.asarray(A, float)
        self.mu = A.mean(0)
        self.sd = A.std(0) + 1e-12
        return self

    def fwd(self, A):
        return (np.asarray(A, float) - self.mu) / self.sd

    def inv(self, A):
        return np.asarray(A, float) * self.sd + self.mu


def quantiles_from_samples(samples, alpha):
    """Lower/upper central (alpha/2, 1-alpha/2) quantiles from (N,3,S) draws."""
    lo = np.quantile(samples, alpha / 2.0, axis=2)
    hi = np.quantile(samples, 1.0 - alpha / 2.0, axis=2)
    return lo, hi


def _phys_from_std(std_samples, ystd):
    """Map standardized-space draws (N,3,S) back to physical units."""
    return std_samples * ystd.sd[None, :, None] + ystd.mu[None, :, None]


# --------------------------------------------------------------------------- #
# torch MLP baselines
# --------------------------------------------------------------------------- #
def _make_mlp(torch, nn, n_in, n_out, hidden=64, seed=None):
    # Seed BEFORE constructing the layers so weight initialisation is a reproducible
    # function of ``seed`` -- NOT of torch's process-global RNG state + call order.
    # (Gauge-CI determinism fix: the legacy code seeded inside ``_train`` AFTER the
    # net was already built, so init drew from the global generator and never tracked
    # the run seed -- the "multi-seed secretly single-init" trap.)
    if seed is not None:
        torch.manual_seed(int(seed))
    return nn.Sequential(
        nn.Linear(n_in, hidden), nn.ReLU(),
        nn.Linear(hidden, hidden), nn.ReLU(),
        nn.Linear(hidden, n_out),
    )


def _train(torch, net, Xtr, Ytr, loss_fn, epochs=600, lr=2e-3, wd=1e-5):
    # Net init is already seeded in ``_make_mlp``; full-batch Adam consumes no RNG.
    opt = torch.optim.Adam(net.parameters(), lr=lr, weight_decay=wd)
    Xt = torch.tensor(Xtr, dtype=torch.float32)
    Yt = torch.tensor(Ytr, dtype=torch.float32)
    net.train()
    for _ in range(epochs):
        opt.zero_grad()
        loss = loss_fn(net(Xt), Yt)
        loss.backward()
        opt.step()
    net.eval()
    return float(loss.item())


class SingleGaussianPNN:
    name = "PNN-Gaussian"

    def __init__(self, seed=0, hidden=64, epochs=600):
        self.seed, self.hidden, self.epochs = seed, hidden, epochs

    def fit(self, X, Y):
        import torch
        import torch.nn as nn
        self.xs, self.ys = _Std().fit(X), _Std().fit(Y)
        Xs, Ys = self.xs.fwd(X), self.ys.fwd(Y)
        self.net = _make_mlp(torch, nn, X.shape[1], 6, self.hidden, seed=self.seed)

        def nll(out, y):
            mean, logvar = out[:, :3], out[:, 3:].clamp(-8, 8)
            return (0.5 * (logvar + (y - mean) ** 2 / torch.exp(logvar))).mean()

        self.final_loss = _train(torch, self.net, Xs, Ys, nll,
                                 epochs=self.epochs)
        return self

    def predict_samples(self, X, n=N_SAMPLES, rng=None):
        import torch
        rng = rng or np.random.default_rng(self.seed)
        with torch.no_grad():
            out = self.net(torch.tensor(self.xs.fwd(X), dtype=torch.float32)).numpy()
        mean, std = out[:, :3], np.exp(0.5 * np.clip(out[:, 3:], -8, 8))
        z = rng.standard_normal((X.shape[0], 3, n))
        draws = mean[:, :, None] + std[:, :, None] * z
        return _phys_from_std(draws, self.ys)


class MDNDeepEnsemble:
    """Deep ensemble of Mixture Density Networks (Casali-style)."""
    name = "MDN-DeepEnsemble"

    def __init__(self, n_members=5, n_comp=2, seed=0, hidden=64, epochs=600):
        self.M, self.K = n_members, n_comp
        self.seed, self.hidden, self.epochs = seed, hidden, epochs

    def fit(self, X, Y):
        import torch
        import torch.nn as nn
        self.xs, self.ys = _Std().fit(X), _Std().fit(Y)
        Xs, Ys = self.xs.fwd(X), self.ys.fwd(Y)
        K = self.K
        out_dim = 3 * K * 3  # per param: K means, K logvars, K logits

        def mdn_nll(out, y):
            o = out.reshape(-1, 3, 3 * K)
            mu = o[:, :, :K]
            logvar = o[:, :, K:2 * K].clamp(-8, 8)
            logit = o[:, :, 2 * K:]
            logw = torch.log_softmax(logit, dim=2)
            yy = y[:, :, None]
            comp = -0.5 * (logvar + (yy - mu) ** 2 / torch.exp(logvar)
                           + np.log(2 * np.pi))
            return -(torch.logsumexp(logw + comp, dim=2)).mean()

        self.nets, self.losses = [], []
        for m in range(self.M):
            net = _make_mlp(torch, nn, X.shape[1], out_dim, self.hidden,
                            seed=self.seed * 100 + m)
            loss = _train(torch, net, Xs, Ys, mdn_nll, epochs=self.epochs)
            self.nets.append(net)
            self.losses.append(loss)
        return self

    def _member_params(self, net, X):
        import torch
        K = self.K
        with torch.no_grad():
            o = net(torch.tensor(self.xs.fwd(X), dtype=torch.float32)).numpy()
        o = o.reshape(-1, 3, 3 * K)
        mu = o[:, :, :K]
        var = np.exp(np.clip(o[:, :, K:2 * K], -8, 8))
        logit = o[:, :, 2 * K:]
        w = np.exp(logit - logit.max(2, keepdims=True))
        w = w / w.sum(2, keepdims=True)
        return mu, var, w  # each (N,3,K)

    def predict_samples(self, X, n=N_SAMPLES, rng=None):
        rng = rng or np.random.default_rng(self.seed + 1)
        N = X.shape[0]
        per = n // self.M
        chunks = []
        for net in self.nets:
            mu, var, w = self._member_params(net, X)   # (N,3,K)
            # choose component per draw via the mixture weights
            cdf = np.cumsum(w, axis=2)
            u = rng.uniform(size=(N, 3, per))
            comp = (u[..., None, :] > cdf[..., :, None]).sum(axis=2)  # (N,3,per)
            comp = np.clip(comp, 0, self.K - 1)
            mu_sel = np.take_along_axis(mu, comp, axis=2)
            sd_sel = np.sqrt(np.take_along_axis(var, comp, axis=2))
            chunks.append(mu_sel + sd_sel * rng.standard_normal((N, 3, per)))
        draws = np.concatenate(chunks, axis=2)
        return _phys_from_std(draws, self.ys)

    def uncertainty_split(self, X):
        """Mean aleatoric vs epistemic predictive variance (physical units)."""
        means, vars_ = [], []
        for net in self.nets:
            mu, var, w = self._member_params(net, X)
            m = (w * mu).sum(2)                       # member predictive mean
            v = (w * (var + mu ** 2)).sum(2) - m ** 2  # member predictive var
            means.append(m)
            vars_.append(v)
        means = np.stack(means)   # (M,N,3) std-space
        vars_ = np.stack(vars_)
        aleatoric = vars_.mean(0)              # E[var]
        epistemic = means.var(0)              # Var[mean]
        scale = self.ys.sd ** 2
        return aleatoric * scale, epistemic * scale


class DeepEnsemblePoint:
    """Plain deep ensemble of point regressors (epistemic spread only)."""
    name = "DeepEnsemble-Point"

    def __init__(self, n_members=5, seed=0, hidden=64, epochs=600):
        self.M, self.seed, self.hidden, self.epochs = n_members, seed, hidden, epochs

    def fit(self, X, Y):
        import torch
        import torch.nn as nn
        self.xs, self.ys = _Std().fit(X), _Std().fit(Y)
        Xs, Ys = self.xs.fwd(X), self.ys.fwd(Y)
        mse = lambda out, y: ((out - y) ** 2).mean()
        self.nets = []
        for m in range(self.M):
            net = _make_mlp(torch, nn, X.shape[1], 3, self.hidden,
                            seed=self.seed * 50 + m)
            _train(torch, net, Xs, Ys, mse, epochs=self.epochs)
            self.nets.append(net)
        return self

    def predict_samples(self, X, n=N_SAMPLES, rng=None):
        import torch
        rng = rng or np.random.default_rng(self.seed + 2)
        preds = []
        with torch.no_grad():
            xt = torch.tensor(self.xs.fwd(X), dtype=torch.float32)
            for net in self.nets:
                preds.append(net(xt).numpy())
        preds = np.stack(preds)                 # (M,N,3) std-space
        mean, std = preds.mean(0), preds.std(0) + 1e-9
        z = rng.standard_normal((X.shape[0], 3, n))
        draws = mean[:, :, None] + std[:, :, None] * z
        return _phys_from_std(draws, self.ys)


# --------------------------------------------------------------------------- #
# Bayesian IVIM via vectorized random-walk Metropolis
# --------------------------------------------------------------------------- #
def _nlls_init_and_noise(signals, b):
    """Per-voxel NLLS point estimate (D,D*,f), S0, and residual-std sigma."""
    N = signals.shape[0]
    theta = np.empty((N, 3))
    s0 = np.empty(N)
    sigma = np.empty(N)
    for i, s in enumerate(signals):
        e = fit_nlls(s, b, return_s0=True)
        theta[i] = (e["D"], e["Dstar"], e["f"])
        s0[i] = e["S0"]
        model = ivim_signal(b, e["D"], e["Dstar"], e["f"], S0=e["S0"])
        sigma[i] = max(float(np.std(s - model)), 1e-3)
    return theta, s0, sigma


class BayesianIVIM_MCMC:
    """Per-voxel Bayesian fit of the bi-exponential model.

    Component-wise (Metropolis-within-Gibbs) random walk, vectorized across
    voxels, with per-dimension step adaptation during burn-in -- this handles
    the wildly different posterior scales of the well-constrained D vs the
    ill-posed D*. The noise level is supplied (we feed the TRUE per-voxel sigma
    so the Bayesian baseline is evaluated at its best, not handicapped by a
    plug-in noise estimate). Uniform priors over the physiological ranges.
    """
    name = "Bayesian-MCMC"

    def __init__(self, seed=0, n_samples=N_SAMPLES, burn=1500, thin=4,
                 step_scale=0.3, target_accept=0.30):
        self.seed, self.n_samples = seed, n_samples
        self.burn, self.thin, self.step_scale = burn, thin, step_scale
        self.target_accept = target_accept
        self.accept_rate = None
        self.accept_per_dim = None
        self.convergence_drift = None  # |mean(1st half)-mean(2nd half)|/post-std

    def fit(self, X, Y):  # per-voxel; nothing to learn from train
        return self

    def predict_samples_for(self, signals, b, theta0, s0_0, sigma, rng):
        Nv = signals.shape[0]
        lb = np.array([D_RANGE[0], DSTAR_RANGE[0], F_RANGE[0], 0.3])
        ub = np.array([D_RANGE[1], DSTAR_RANGE[1], F_RANGE[1], 2.0])
        step = self.step_scale * (ub - lb)              # per-dim, adapted below
        state = np.column_stack([theta0[:, 0], theta0[:, 1], theta0[:, 2], s0_0])
        state = np.clip(state, lb + 1e-9, ub - 1e-9)
        inv2s2 = 1.0 / (2.0 * sigma ** 2)

        def loglik(th):
            model = ivim_signal(b[None, :], th[:, 0:1], th[:, 1:2], th[:, 2:3],
                                S0=th[:, 3:4])
            return -np.sum((model - signals) ** 2, axis=1) * inv2s2

        cur_ll = loglik(state)
        total = self.burn + self.n_samples * self.thin
        kept = []
        adapt_win = 100
        win_acc = np.zeros(4)
        win_n = 0
        glob_acc = np.zeros(4)
        for it in range(total):
            for d in range(4):
                prop = state.copy()
                prop[:, d] = state[:, d] + step[d] * rng.standard_normal(Nv)
                inb = (prop[:, d] >= lb[d]) & (prop[:, d] <= ub[d])
                pll = np.where(inb, loglik(prop), -np.inf)
                a = np.log(rng.uniform(size=Nv)) < (pll - cur_ll)
                state[a] = prop[a]
                cur_ll[a] = pll[a]
                rate = a.mean()
                win_acc[d] += rate
                glob_acc[d] += rate
            win_n += 1
            # per-dim step adaptation during burn-in only
            if it < self.burn and win_n == adapt_win:
                r = win_acc / win_n
                step *= np.exp(0.6 * (r - self.target_accept))
                step = np.clip(step, 1e-6 * (ub - lb), 1.0 * (ub - lb))
                win_acc[:] = 0.0
                win_n = 0
            if it >= self.burn and (it - self.burn) % self.thin == 0:
                kept.append(state[:, :3].copy())
        samples = np.stack(kept, axis=2)                # (Nv,3,n_samples)
        self.accept_per_dim = glob_acc / total
        self.accept_rate = float(self.accept_per_dim.mean())
        # convergence diagnostic: normalized drift between chain halves
        h = samples.shape[2] // 2
        m1, m2 = samples[:, :, :h].mean(2), samples[:, :, h:].mean(2)
        psd = samples.std(2) + 1e-12
        self.convergence_drift = np.mean(np.abs(m1 - m2) / psd, axis=0)  # (3,)
        return samples


# --------------------------------------------------------------------------- #
# build + cache all predictions for the benchmark
# --------------------------------------------------------------------------- #
def build_predictions(alphas=ALPHAS, seed=DEFAULT_SEED, cache_path=None,
                      force=False, verbose=True):
    """Train/run every baseline + conformal bases on the matched cohort.

    Returns a dict of per-(method, split, alpha) interval bounds plus the
    conformal-base predictions, cached to ``cache_path`` (regenerable; not
    committed). A **pure function of ``seed``**: the cohort, NLLS, HGB and MCMC
    streams are seeded as before, and (Gauge-CI) the torch NN bases now take
    per-seed init streams derived from ``seed`` so their weight initialisation
    varies across seeds too. The cache path defaults to a **seed-specific** file
    so a multi-seed sweep never silently reuses another seed's predictions.
    """
    if cache_path is None:
        cache_path = os.path.join(_RESULTS_DIR, f"predictions_seed{int(seed)}.pkl")
    if (not force) and os.path.exists(cache_path):
        with open(cache_path, "rb") as fh:
            return pickle.load(fh)

    os.makedirs(_RESULTS_DIR, exist_ok=True)
    import torch
    torch.use_deterministic_algorithms(True, warn_only=True)

    cohort = generate_cohort(4000, 2000, 3000, seed=seed)
    b = cohort.b
    Xtr, Ytr = cohort.signals["train"], cohort.params["train"]
    splits = {"cal": cohort.signals["cal"], "test": cohort.signals["test"]}

    R = {
        "meta": {"alphas": list(alphas), "snr_grid": list(DEFAULT_SNR_GRID),
                 "sizes": cohort.sizes, "seed": seed, "params": PARAM_NAMES,
                 "b": b},
        "cal_true": cohort.params["cal"], "test_true": cohort.params["test"],
        "cal_snr": cohort.snr["cal"], "test_snr": cohort.snr["test"],
        "methods": [], "diag": {},
    }

    # --- conformal bases (Gauge 01) -------------------------------------
    t0 = time.time()
    nlls = {}
    mcmc_init = {}
    for sname, Xs in splits.items():
        theta, s0, sigma = _nlls_init_and_noise(Xs, b)
        nlls[sname] = theta
        mcmc_init[sname] = (theta, s0, sigma)
        # Store the NLLS S0 and residual-noise estimate so Gauge 03 can build a
        # label-free SNR proxy (SNR_hat = S0_hat / sigma_hat) without re-fitting.
        R[f"nlls_s0_{sname}"] = s0
        R[f"nlls_sigma_{sname}"] = sigma
    R["nlls_cal"], R["nlls_test"] = nlls["cal"], nlls["test"]
    if verbose:
        print(f"[build] NLLS bases done ({time.time()-t0:.0f}s)")

    levels = sorted({a / 2 for a in alphas} | {1 - a / 2 for a in alphas})
    qreg = IVIMQuantileRegressor(levels, random_state=0).fit(Xtr, Ytr)
    R["hgb_levels"] = levels
    for sname, Xs in splits.items():
        for j in range(3):
            for q in levels:
                R[f"hgb_{sname}_{j}_{q:.4f}"] = qreg.predict_quantile(Xs, j, q)
    if verbose:
        print(f"[build] HGB quantile base done ({time.time()-t0:.0f}s)")

    # --- model-based baselines -----------------------------------------
    # Per-seed NN init streams (independent + reproducible from the run seed via
    # SeedSequence). Replaces the legacy hardcoded ``seed=0`` so the deep-ensemble
    # initialisation epistemic spread is actually resampled across seeds, not held
    # fixed at one (favourable) init -- the B.0-probe finding.
    pnn_s, mdn_s, de_s = (int(s.generate_state(1)[0])
                          for s in np.random.SeedSequence(seed).spawn(3))
    nn_models = [
        SingleGaussianPNN(seed=pnn_s),
        MDNDeepEnsemble(n_members=5, n_comp=2, seed=mdn_s),
        DeepEnsemblePoint(n_members=5, seed=de_s),
    ]
    for model in nn_models:
        tt = time.time()
        model.fit(Xtr, Ytr)
        R["methods"].append(model.name)
        for sname, Xs in splits.items():
            samp = model.predict_samples(Xs, rng=np.random.default_rng(seed + 7))
            for a in alphas:
                lo, hi = quantiles_from_samples(samp, a)
                R[f"{model.name}_{sname}_lo_{a}"] = lo
                R[f"{model.name}_{sname}_hi_{a}"] = hi
        if isinstance(model, MDNDeepEnsemble):
            al, ep = model.uncertainty_split(splits["test"])
            R["diag"]["mdn_aleatoric_mean"] = al.mean(0)
            R["diag"]["mdn_epistemic_mean"] = ep.mean(0)
        if verbose:
            print(f"[build] {model.name} done ({time.time()-tt:.0f}s)")

    # --- Bayesian MCMC (fed TRUE per-voxel noise sigma = 1/SNR) ----------
    tt = time.time()
    bayes = BayesianIVIM_MCMC(seed=0)
    R["methods"].append(bayes.name)
    snr_of = {"cal": R["cal_snr"], "test": R["test_snr"]}
    for sname, Xs in splits.items():
        theta, s0, _ = mcmc_init[sname]
        sigma_true = 1.0 / np.asarray(snr_of[sname], float)   # S0 = 1
        samp = bayes.predict_samples_for(
            Xs, b, theta, s0, sigma_true, np.random.default_rng(seed + 11))
        for a in alphas:
            lo, hi = quantiles_from_samples(samp, a)
            R[f"{bayes.name}_{sname}_lo_{a}"] = lo
            R[f"{bayes.name}_{sname}_hi_{a}"] = hi
    R["diag"]["mcmc_accept_rate"] = bayes.accept_rate
    R["diag"]["mcmc_accept_per_dim"] = bayes.accept_per_dim
    R["diag"]["mcmc_convergence_drift"] = bayes.convergence_drift
    if verbose:
        print(f"[build] {bayes.name} done ({time.time()-tt:.0f}s, "
              f"accept={bayes.accept_rate:.2f})")

    with open(cache_path, "wb") as fh:
        pickle.dump(R, fh)
    return R


# --------------------------------------------------------------------------- #
# GATE 0 -- protocol check: each baseline emits intervals; print raw coverage
# --------------------------------------------------------------------------- #
def main():
    from gauge.conformal import empirical_coverage, interval_width
    R = build_predictions(force=os.environ.get("GAUGE_FORCE") == "1")
    alphas = R["meta"]["alphas"]
    test_true = R["test_true"]

    print("=" * 78)
    print("CP0 / GATE 0 -- model-based baselines, matched protocol sanity")
    print("=" * 78)
    m = R["meta"]
    print(f"cohort seed {m['seed']}  sizes {m['sizes']}  "
          f"SNR grid {m['snr_grid']}  b-values {len(m['b'])}")
    print(f"baselines: {R['methods']}")
    print(f"alpha sweep {alphas}  (nominal = 1 - alpha)")
    if "mcmc_accept_rate" in R["diag"]:
        apd = R["diag"].get("mcmc_accept_per_dim")
        drift = R["diag"].get("mcmc_convergence_drift")
        print(f"MCMC mean acceptance: {R['diag']['mcmc_accept_rate']:.2f}  "
              f"per-dim (D,D*,f,S0): "
              f"{np.array2string(apd, precision=2)}")
        print(f"MCMC convergence drift |half1-half2|/post-std (D,D*,f): "
              f"{np.array2string(drift, precision=3)}  (<<1 = mixed)")
    if "mdn_aleatoric_mean" in R["diag"]:
        al, ep = R["diag"]["mdn_aleatoric_mean"], R["diag"]["mdn_epistemic_mean"]
        print("MDN uncertainty split (mean predictive var, physical units):")
        for j, p in enumerate(PARAM_NAMES):
            tot = al[j] + ep[j]
            print(f"   {p:>3}: aleatoric={al[j]:.3e}  epistemic={ep[j]:.3e}  "
                  f"epistemic frac={ep[j]/tot:.2f}")
    print("-" * 78)
    print("RAW marginal coverage on TEST (sanity: each baseline emits "
          "per-parameter intervals)")
    print(f"{'method':>20} | {'param':>4} | " +
          " | ".join(f"a={a}" for a in alphas))
    print("-" * 78)
    for name in R["methods"]:
        for j, p in enumerate(PARAM_NAMES):
            cells = []
            for a in alphas:
                lo = R[f"{name}_test_lo_{a}"][:, j]
                hi = R[f"{name}_test_hi_{a}"][:, j]
                cells.append(f"{empirical_coverage(lo, hi, test_true[:, j]):.3f}")
            print(f"{name:>20} | {p:>4} | " + " | ".join(c.rjust(5) for c in cells))
        print("-" * 78)
    print("GATE 0: PASS -- all baselines emit per-parameter intervals on the "
          "matched cohort.")
    print("(Whether they are WELL-calibrated is the CP1 question, reported "
          "honestly there.)")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    raise SystemExit(main())
