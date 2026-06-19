# Lattice — a UQ-calibration digital reference object (DRO) for IVIM

*Research software (MIT). A reusable **resource**, not a research result.*

Lattice is a **digital reference object** for *uncertainty-quantification
calibration* in IVIM (intravoxel incoherent motion) diffusion MRI. It packages
three things that, together, let any UQ method be benchmarked on a common,
reproducible, ground-truthed footing:

1. **Reference parameter distributions** — physiologically-grounded priors over
   IVIM ground truth `(D, D*, f)` (a `uniform` prior and a correlated
   `realistic` Gaussian-copula prior).
2. **Alternative-model signal generators** — five clean-room forward models
   (bi-exponential plus four *misspecification* families) that each reduce to the
   bi-exponential model at a continuity limit, so you can test calibration under
   controlled model error.
3. **A standardized calibration-evaluation interface** — one estimator contract
   and an adapter that hands a cohort's ground truth and a method's predicted
   quantiles straight to a quantile calibration scorer (canonically
   [Caliper](../Caliper/): coverage / ECE / sharpness).

Everything here is **synthetic and PHI-free**, fully reproducible from an integer
seed, and **depends on no publication** — the cohorts and generators are solid
now. (Only an eventual *citable* release waits on the companion papers' DOIs; see
[Citable release](#citable-release).)

## Why Lattice exists (positioning)

Lattice is deliberately distinct from its neighbours — see
[`docs/POSITIONING.md`](docs/POSITIONING.md):

| Neighbour | What it is | How Lattice differs |
|-----------|-----------|---------------------|
| **OSIPI TF2.4** (`Fashion/` upstream) | Accuracy-focused fitting benchmark + reference DRO for *point-estimate correctness* | Lattice targets **UQ calibration** (coverage/sharpness under misspecification), not accuracy |
| **Gauge's cohort** | Synthetic cohort used *internally* for one paper's results | Lattice is a **packaged, versioned, reusable** DRO with a documented API |
| **Caliper** | The calibration **scorer** (the ruler) | Lattice is the **data** (the reference object) the scorer consumes |

The distinct contribution is the **alternative-model generator families + a
standardized calibration-evaluation API over a reference parameter
distribution** — not present, as a packaged resource, anywhere else in the
family.

## Install

```bash
pip install -e .                 # core: numpy only
pip install -e ".[selfcheck]"    # + scipy, for the NLLS round-trip self-check
pip install -e ".[external-data]"# + requests, for optional OSIPI download-on-demand
```

## Quickstart

```python
import lattice

# A bi-exponential cohort (the baseline reference object)
cohort = lattice.make_cohort("biexp", n=2000, snr=50, seed=lattice.DEFAULT_SEED)
cohort.params      # (2000, 3) ground truth, columns (D, D*, f)
cohort.signals     # (2000, 22) noisy multi-b signals
cohort.bvalues     # the 22-point b-value schedule

# A misspecification cohort: stretched-exponential perfusion
alt = lattice.make_cohort("stretched", n=2000, extra={"beta": 0.7})

# Score any UQ method (must implement predict_quantiles -> (n, n_params, n_levels))
q_pred = my_method.predict_quantiles(cohort.signals, lattice.DEFAULT_QUANTILE_LEVELS)
payload = lattice.to_scorer_inputs(cohort, q_pred)   # ready for caliper.metrics.score_quantiles
```

## Cohort families and their continuity limits

Each non-baseline family reduces to bi-exponential at a continuity limit. Three
reduce **exactly** (verified residual `0.000e+00`); gamma reduces asymptotically.

| Family | Signal (perfusion term) | Extra params | Continuity limit → bi-exp | Verified residual |
|--------|-------------------------|--------------|---------------------------|-------------------|
| `biexp` | `f·e^{-bD*}` | — | identity | `0.000e+00` |
| `dispersion_gamma` | `f·(1+bμ/k)^{-k}` | `k` | `k→∞` (asymptotic) | `3.9e-09` (k=1e8) |
| `dispersion_lognormal` | `f·E[e^{-bD*}]`, `D*~LogN(μ,cv)` | `cv` | `cv=0` (exact) | `0.000e+00` |
| `stretched` | `f·e^{-(bD*)^β}` | `β` | `β=1` (exact) | `0.000e+00` |
| `triexp` | `f(1-g)e^{-bD*}+fg·e^{-bD*₂}` | `Dstar2_mult, g` | `g=0` (exact) | `0.000e+00` |

## Self-consistency (the gates)

- **Clean round-trip** — fit a noise-free bi-exp cohort by NLLS and recover the
  generating `(D, D*, f)`. Verified **max relative error `7.99e-10`** (gate
  `< 1e-2`).
- **Continuity** — every family at its continuity limit matches the
  bi-exponential cohort built from identical ground truth (table above).

Reproduce: `python examples/continuity_demo.py` and `pytest -q`.

## Layout

- `lattice/generators.py` — the five clean-room forward models + noise.
- `lattice/cohort.py` — the `Cohort` schema, priors, `make_cohort`, continuity helper.
- `lattice/evaluate.py` — the estimator contract + `to_scorer_inputs` adapter.
- `lattice/selfcheck.py` — NLLS round-trip self-validation (scipy; `selfcheck` extra).
- `lattice/osipi.py` + `scripts/fetch_osipi.py` — optional OSIPI download-on-demand + provenance.
- `lattice/publication.py` — citation gate for an eventual citable release (OFF by default).
- `examples/`, `docs/`, `tests/`.

## One-way dependency

Caliper and the companion papers **consume** Lattice; Lattice **imports nothing
back**. The `lattice` core package has no dependency on any sibling project.
Scoring "via Caliper" is demonstrated only in `examples/evaluate_with_caliper.py`,
where Caliper is an *optional, example-only* import. See
[`docs/CLEANROOM.md`](docs/CLEANROOM.md).

## Citable release

The DRO is usable now. An eventual JOSS/Zenodo citable release is gated **only**
on its *citations* — the companion papers' DOIs (Fashion, Gauge, Minos) — never
on the DRO content. `lattice.publication` renders those citations as honest
`@unpublished` placeholders until a real `paper_doi` is filled in.

## License

MIT — see [`LICENSE`](LICENSE). All in-tree data is synthetic; any external
(OSIPI) data is fetched on demand under its own CC-BY-4.0 license and never
redistributed here.
