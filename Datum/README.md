# Datum — the IVIM calibration benchmark

**Datum** is a benchmark for *uncertainty-quantification calibration* in
intravoxel incoherent motion (IVIM) diffusion MRI. It answers a single, fixed
question for any method that reports an error bar on IVIM parameters: **is that
error bar honest?** — scored on **Fashion's calibration ruler** (the coverage /
ECE / sharpness standard), as packaged read-only by **Caliper**, over a synthetic
data substrate.

> ⚠️ **PROVISIONAL.** Datum is built on Fashion's calibration ruler, which is **in
> review at *NMR in Biomedicine*** (retooled, boundary-railing-first). The retooled
> Fashion demotes the ruler to a **scoped secondary** (ground-truth/synthetic only)
> reported under the **honest CRLB** convention — Datum's substrate is synthetic, so
> the scope holds. Every reference number Datum produces by scoring through that
> ruler is flagged **PROVISIONAL** and must not be presented as a final reference
> value until the ruler locks (DOI assigned). See [`ASSUMPTIONS.md`](ASSUMPTIONS.md)
> and [`datum/manifest.py`](datum/manifest.py); re-validate with
> `python revalidate.py`.

## What Datum is (and what it is not)

A *benchmark* = a fixed task + curated baselines + reference numbers + a way to
submit a new method. Datum is distinct from the three things it sits next to:

| | What it is | Datum adds |
|---|---|---|
| **Caliper** | the calibration *ruler* + an internal demo sweep (explicitly *not* a citable benchmark) | a **frozen, versioned task**, a **curated multi-paradigm baseline panel**, **reference numbers with bootstrap CIs**, and a **submission/scoring interface** |
| **Lattice** | a data *substrate* / DRO (Datum's primary substrate) | the *task + baselines + numbers + interface* layered on that substrate |
| **OSIPI** | fitting algorithms + a DRO scored on **point accuracy** (bias/RMSE) | scoring on **uncertainty-calibration honesty** via Fashion's ruler-as-standard |

Datum also doubles as the concrete artifact behind **Fashion's "ruler-as-standard"
differentiation from Casali**: Casali (2026) is a model-based UQ *method*; Fashion's
ruler is the *standard* that measures any method's honesty, and Datum is where that
standard is applied as a leaderboard.

## The task (frozen, versioned — `datum/task.py::CURRENT_TASK`, v2)

- **Substrate:** the **Lattice IVIM DRO** (`lattice.make_cohort`, `seed = 20260619`),
  train/cal/test — the intended substrate, now built (v2 swapped it in for the v1
  Gauge bootstrap, which is kept runnable). An **OSIPI DRO** is the external-validation
  substrate. All synthetic and PHI-free.
- **Contract:** a method predicts quantiles `(n, 3, L)` for `(D, D*, f)` at fixed
  `quantile_levels`.
- **Metrics (Fashion's ruler):** coverage, coverage-gap, ECE, sharpness, pinball,
  interval-score — **marginal and per-D\* tercile** (the high-D\* identifiability
  wall), with bootstrap CIs on the load-bearing numbers.
- **Baseline panel (all reused, none reinvented):** NLLS+Gaussian σ (the known
  under-coverer), NLLS+residual-bootstrap, segmented reference, split-conformal,
  CQR, Mondrian-CQR, and MAF posterior quantiles.

## Layout

- `datum/` — `task.py` (the frozen spec), `substrate.py` (read-only Gauge/OSIPI/Lattice
  adapters), `ruler.py` (read-only adapter over `caliper.metrics`), `baselines.py`
  (the curated panel registry), `run.py` (the benchmark runner), `ci.py` (bootstrap
  CIs), `convert.py` (Gauge↔Caliper convention), `submit.py` (the submission/scoring
  interface), `osipi_fetch.py` (external-validation fetch), `manifest.py` +
  `provisional.py` (assumption pins + PROVISIONAL stamping), `_paths.py` (bootstrap).
- `examples/submit_demo.py` — score a new method on Datum (runnable).
- `docs/` — [`benchmark.md`](docs/benchmark.md) (the task), [`submitting.md`](docs/submitting.md)
  (the contract), [`ruler-as-standard.md`](docs/ruler-as-standard.md) (the Casali
  framing), [`release.md`](docs/release.md) (citable-release path, gated).
- `results/` — `REFERENCE.md` + `reference_numbers.csv` (PROVISIONAL) + OSIPI provenance.
- `tests/` — import/manifest/task/convert/run/submit/osipi gates. `revalidate.py` —
  one-command re-validation. `ASSUMPTIONS.md` — SOLID-vs-PROVISIONAL + pinned ruler.

## Quickstart

```bash
# Inside the ResearchProj monorepo (Caliper/ and Gauge/ are siblings of Datum/):
pip install -e Datum              # numpy only; Caliper/Gauge resolved via sibling bootstrap
python Datum/revalidate.py        # check pins + prove substrate→ruler pipeline resolves
python -m datum.run               # produce the (PROVISIONAL) reference numbers
python -m datum.osipi_fetch       # (optional) fetch the OSIPI DRO for external validation
python Datum/examples/submit_demo.py   # score a new method on Datum (worked example)
pytest Datum/tests                # gates
```

## Reference numbers (PROVISIONAL)

Produced by [`results/REFERENCE.md`](results/REFERENCE.md) /
[`results/reference_numbers.csv`](results/reference_numbers.csv) via
`python -m datum.run`. **Scored on Fashion's in-review ruler — PROVISIONAL, no
tuning.** The headline (Lattice DRO, D\* coverage at nominal 0.90): raw error bars
under-cover D\* (NLLS-Gaussian gap −0.089, segmented −0.649, MAF −0.038); conformal
restores marginal coverage (split +0.005, CQR −0.007, MAF+CQR −0.001, Mondrian
−0.003). The high-D\* identifiability wall persists for marginal split-conformal
(coverage 0.80) and is recovered by CQR/Mondrian only by inflating width (~19→185).
The OSIPI DRO external validation reproduces the story on an independent synthetic
phantom (D\* gap −0.130 raw → −0.002 split, +0.012 CQR).

## Status

CP1–CP3 complete: subrepo embedded (own clean synthetic-only history, mirroring
Minos), registered in the monorepo README, read-only imports resolve, assumptions
manifest pins the ruler; the curated baseline panel has been run through the ruler on
the **Lattice DRO** (primary; + OSIPI DRO external validation, Gauge bootstrap kept
runnable) to produce reference numbers with bootstrap CIs; and a
**submission/scoring interface** (`datum.submit`), worked
example, docs, and the ruler-as-standard / Casali framing are in place — **every
ruler-derived number PROVISIONAL** until the ruler locks. A citable JOSS/Zenodo
release is documented ([`docs/release.md`](docs/release.md)) but **not executed** —
gated on Fashion's ruler locking, exactly like Caliper.

## License

MIT — see [`LICENSE`](LICENSE).
