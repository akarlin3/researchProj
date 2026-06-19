"""CP2/CP3 clean-data layer -- a synthetic IVIM cohort with Fashion-calibrated uncertainty.

PROVISIONAL: consumes Fashion's uncertainty generators (in review). See ../ASSUMPTIONS.md.

Strictly synthetic + open data only (IP gate): signals come from Fashion's IVIM simulator
(`uq.ivim_simulator`, pancreatic-anchor priors). No pancData3 / MSK / clinical data.

Produces, for a cohort of voxels with KNOWN truths:
    truth[param]   ground-truth D / Dstar / f         (synthetic; the label)
    est[:, k]      Fashion point estimate             (the reported centre mu)
    sigma[:, k]    Fashion reported uncertainty       (the reported error bar)
    lo/hi[:, k]    skew-aware credible interval        (MCMC generator only)
for param index k in {D:0, Dstar:1, f:2}. This is the (theta, mu, sigma) the
decision-calibration gap (gap_applied.py) and the validity monitor (monitor_applied.py)
operate on.

Results are cached under ../results/cohorts/ as .npz so the (slow) MCMC refit runs once.
"""
from __future__ import annotations

import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # future/
import _paths  # noqa: E402

_paths.add_fashion()

from uq.ivim_fit import make_model  # noqa: E402  (Fashion)
from uq.ivim_simulator import simulate_random, B_SCHEMES, PANCREAS_PRIORS  # noqa: E402
from uq.bayesian import laplace_uncertainty, mcmc_uncertainty  # noqa: E402
from uq.bootstrap import bootstrap_cell  # noqa: E402

PARAM_INDEX = {"D": 0, "Dstar": 1, "f": 2}
DEFAULT_ALGORITHM = "OGC_AmsterdamUMC_biexp"   # standard NLLS biexp (Amsterdam UMC, open OSIPI)
CACHE_DIR = os.path.join(_paths._FUTURE, "results", "cohorts")

# D* and the other params live on very different numeric scales; the decision-calibration
# gap is scale-invariant in the parameter, so downstream code rescales to O(1) units.
# These multipliers turn raw mm^2/s into convenient units (D, D* -> 1e-3 mm^2/s; f stays).
UNIT_SCALE = {"D": 1e3, "Dstar": 1e3, "f": 1.0}


def _cohort_path(generator, n, snr, scheme, seed, algorithm):
    os.makedirs(CACHE_DIR, exist_ok=True)
    tag = f"{generator}_n{n}_snr{int(snr)}_{scheme}_seed{seed}_{algorithm}"
    return os.path.join(CACHE_DIR, tag + ".npz")


def simulate_cohort(generator="laplace", n=3000, snr=40.0, scheme="clinical_sparse",
                    seed=0, algorithm=DEFAULT_ALGORITHM, *, use_cache=True,
                    mcmc_nwalkers=16, mcmc_nsteps=800, mcmc_burn=300, mcmc_thin=5,
                    boot_K=200, verbose=True):
    """Generate (or load) a synthetic IVIM cohort with Fashion uncertainty.

    generator: 'laplace' (fast, Hessian SD), 'mcmc' (skew-aware credible interval --
    Fashion's calibrated recipe), or 'bootstrap' (resampling SD). Returns a dict.
    """
    path = _cohort_path(generator, n, snr, scheme, seed, algorithm)
    if use_cache and os.path.exists(path):
        d = np.load(path, allow_pickle=True)
        out = {k: d[k] for k in d.files}
        out["meta"] = out["meta"].item()
        return out

    b = B_SCHEMES[scheme]
    rng = np.random.default_rng(seed)
    sim = simulate_random(n, b, snr=snr, rng=rng)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = make_model(algorithm, b)
    truth = sim["truth"]  # dict of (n,) arrays D, Dstar, f

    lo = hi = None
    gen_rng = np.random.default_rng(seed + 1)
    if generator == "laplace":
        est, sigma = laplace_uncertainty(model, sim["signals"], b)
    elif generator == "bootstrap":
        est, sigma = bootstrap_cell(model, sim["signals"], b, K=boot_K, rng=gen_rng)
    elif generator == "mcmc":
        est, sigma, lo, hi = mcmc_uncertainty(
            model, sim["signals"], b, nwalkers=mcmc_nwalkers, nsteps=mcmc_nsteps,
            burn=mcmc_burn, thin=mcmc_thin, rng=gen_rng)
    else:
        raise ValueError(f"unknown generator {generator!r}")

    truth_arr = np.stack([truth["D"], truth["Dstar"], truth["f"]], axis=1)  # (n,3)
    out = {
        "truth": truth_arr, "est": np.asarray(est, float), "sigma": np.asarray(sigma, float),
        "lo": np.full_like(est, np.nan) if lo is None else np.asarray(lo, float),
        "hi": np.full_like(est, np.nan) if hi is None else np.asarray(hi, float),
        "bvalues": np.asarray(b, float),
        "meta": dict(generator=generator, n=int(n), snr=float(snr), scheme=scheme,
                     seed=int(seed), algorithm=algorithm,
                     fashion_seed=0, gauge_seed=20260613, provisional=True),
    }
    if use_cache:
        np.savez(path, **{k: (np.array(v, dtype=object) if k == "meta" else v)
                          for k, v in out.items()})
    if verbose:
        fin = np.isfinite(out["sigma"][:, 1]).mean()
        print(f"[data] cohort generator={generator} n={n} snr={snr} {scheme}: "
              f"D* sigma finite frac={fin:.2f}  (cached -> {os.path.basename(path)})")
    return out


def extract_param(cohort, param="Dstar", *, scaled=True, drop_nonfinite=True):
    """Pull (theta, mu, sigma[, lo, hi]) for one parameter, optionally rescaled to O(1).

    Returns dict with theta/mu/sigma (and lo/hi if present), plus the index mask kept.
    """
    k = PARAM_INDEX[param]
    sc = UNIT_SCALE[param] if scaled else 1.0
    theta = cohort["truth"][:, k] * sc
    mu = cohort["est"][:, k] * sc
    sigma = cohort["sigma"][:, k] * sc
    lo = cohort["lo"][:, k] * sc
    hi = cohort["hi"][:, k] * sc
    mask = np.ones(theta.shape, bool)
    if drop_nonfinite:
        mask = np.isfinite(theta) & np.isfinite(mu) & np.isfinite(sigma) & (sigma > 0)
    res = dict(theta=theta[mask], mu=mu[mask], sigma=sigma[mask],
               n_total=int(theta.size), n_kept=int(mask.sum()),
               param=param, unit_scale=sc, mask=mask)
    if np.isfinite(lo[mask]).any():
        res["lo"] = lo[mask]
        res["hi"] = hi[mask]
    return res


if __name__ == "__main__":  # quick smoke
    c = simulate_cohort(generator="laplace", n=200, snr=40, seed=0)
    d = extract_param(c, "Dstar")
    print(f"Dstar: kept {d['n_kept']}/{d['n_total']}; "
          f"theta in [{d['theta'].min():.1f},{d['theta'].max():.1f}] (1e-3 mm^2/s)")
