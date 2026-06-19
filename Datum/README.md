# Datum — the IVIM calibration benchmark

**Datum** is a benchmark for *uncertainty-quantification calibration* in
intravoxel incoherent motion (IVIM) diffusion MRI. It answers a single, fixed
question for any method that reports an error bar on IVIM parameters: **is that
error bar honest?** — scored on **Fashion's calibration ruler** (the coverage /
ECE / sharpness standard), as packaged read-only by **Caliper**, over a synthetic
data substrate.

> ⚠️ **PROVISIONAL.** Datum is built on Fashion's calibration ruler, which is **in
> review at *MRM***. Every reference number Datum produces by scoring through that
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
| **Lattice** | a data *substrate* (not built yet) | the *task + baselines + numbers + interface* layered on a substrate |
| **OSIPI** | fitting algorithms + a DRO scored on **point accuracy** (bias/RMSE) | scoring on **uncertainty-calibration honesty** via Fashion's ruler-as-standard |

Datum also doubles as the concrete artifact behind **Fashion's "ruler-as-standard"
differentiation from Casali**: Casali (2026) is a model-based UQ *method*; Fashion's
ruler is the *standard* that measures any method's honesty, and Datum is where that
standard is applied as a leaderboard.

## The task (frozen, versioned — `datum/task.py::TASK_V1`)

- **Substrate:** Gauge's seeded synthetic IVIM cohort (`seed = 20260613`), train/cal/test.
  Lattice replaces it when built; an OSIPI DRO is wired for external validation.
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
  (the curated panel registry), `manifest.py` (assumption pins + provisional policy),
  `provisional.py` (PROVISIONAL stamping), `_paths.py` (sibling bootstrap for the
  read-only deps).
- `tests/` — import/manifest/task gates (CP1). `revalidate.py` — one-command
  re-validation. `ASSUMPTIONS.md` — the SOLID-vs-PROVISIONAL split and the pinned
  ruler version.

## Quickstart

```bash
# Inside the ResearchProj monorepo (Caliper/ and Gauge/ are siblings of Datum/):
pip install -e Datum            # numpy only; Caliper/Gauge resolved via sibling bootstrap
python Datum/revalidate.py      # check pins + prove substrate→ruler pipeline resolves
pytest Datum/tests              # CP1 gates
```

## Status

CP1 complete: subrepo embedded (own clean synthetic-only history, mirroring Minos),
registered in the monorepo README, read-only imports resolve, assumptions manifest
pins the ruler. **Reference numbers are the CP2 deliverable** and are PROVISIONAL
by construction. A citable JOSS/Zenodo release is documented but **not executed** —
gated on Fashion's ruler locking, exactly like Caliper.

## License

MIT — see [`LICENSE`](LICENSE).
