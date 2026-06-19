# Lattice DRO specification (v0.1)

Lattice is a **digital reference object (DRO)**: a reproducible package of
synthetic IVIM ground truth, the signals generated from it, and the interface to
score uncertainty against it.

## 1. Ground-truth schema

A cohort is a `lattice.cohort.Cohort` dataclass:

| Field | Shape | Meaning |
|-------|-------|---------|
| `family` | str | generator family name |
| `bvalues` | `(n_b,)` | acquisition b-values (s/mm²) |
| `params` | `(n, 3)` | ground truth, columns `(D, Dstar, f)` |
| `extra` | `(n, n_extra)` | family-specific ground-truth params |
| `extra_names` | tuple | names of the `extra` columns |
| `signals_clean` | `(n, n_b)` | noise-free forward signals (S0 = 1) |
| `signals` | `(n, n_b)` | noisy observations at `snr` under `noise` |
| `snr`, `noise`, `seed`, `prior` | scalars | metadata |
| `param_names` | `("D","Dstar","f")` | canonical column order |

Units: `D`, `Dstar` in mm²/s; `f` dimensionless; `b` in s/mm².
Physiological ranges (`lattice.cohort.PARAM_RANGES`): `D ∈ [0.5, 3.0]·10⁻³`,
`Dstar ∈ [10, 100]·10⁻³`, `f ∈ [0.05, 0.40]`.

On disk: `Cohort.save(path)` writes `path.npz` (arrays) + `path.json` (manifest
sufficient to regenerate from the seed).

## 2. Reference parameter distributions (priors)

- `uniform` — independent uniform over each physiological range.
- `realistic` — a Gaussian copula with a documented latent correlation
  (`_REALISTIC_CORR`: mild +D/f, −Dstar/f coupling) and the same uniform
  marginals. Default.

The base `(D, Dstar, f)` draw depends only on `(seed, prior, n)` and is therefore
**identical across families** at a fixed seed — the property the continuity gates
rely on.

## 3. Generator families (the alternative models)

All return normalized `S(b)/S0`. See `lattice.generators` and the README table.
Continuity limits: lognormal `cv=0` (exact), stretched `β=1` (exact), triexp
`g=0` (exact), gamma `k→∞` (asymptotic).

## 4. b-value scheme

`DEFAULT_BVALUES` — 22 points, dense at low b (perfusion) and sparse at high b
(tissue): `[0,5,10,15,20,30,40,50,60,80,100,120,150,200,250,300,400,500,600,700,800,1000]`.
Any custom `bvalues` array may be passed to `make_cohort`.

## 5. Noise model

`rician` (default), `gaussian`, or `none`. `sigma = S0 / snr` defined at b=0.

## 6. Seeds / reproducibility

`make_cohort(seed=...)` drives three offset streams (base, extra, noise) so
results are bit-reproducible and base params are family-invariant.
`DEFAULT_SEED = 20260619` (date-stamped; no wall-clock anywhere).

## 7. Calibration-evaluation interface

The estimator contract (`lattice.evaluate.IVIMQuantileEstimator`):

```python
def predict_quantiles(signals: (n, n_b), q_levels) -> (n, n_params, n_levels)
```

`lattice.evaluate.to_scorer_inputs(cohort, q_pred, q_levels)` returns exactly the
keyword arguments of `caliper.metrics.score_quantiles`, so scoring on Lattice is
one line. Lattice imports no scorer itself (one-way dependency); the canonical
scorer is Caliper (coverage / ECE / sharpness). Tiny dependency-free
`interval_coverage` / `mean_sharpness` helpers exist only for smoke checks.

## 8. Self-consistency gates

- Clean round-trip (bi-exp): max relative recovery error `7.99e-10` (gate `<1e-2`).
- Continuity residuals: `0.000e+00` (exact families), `3.9e-09` (gamma, k=1e8).
