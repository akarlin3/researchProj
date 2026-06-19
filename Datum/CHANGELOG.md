# Changelog

All notable changes to Datum are recorded here. Datum follows the monorepo's
provisional-until-the-ruler-locks discipline; see `ASSUMPTIONS.md`.

## [0.1.0] — 2026-06-19 — CP1: subrepo + README + manifest

### Added
- Datum package scaffolding: `task.py` (frozen `TASK_V1` spec), `substrate.py`
  (read-only Gauge cohort + OSIPI DRO + Lattice swap-in stub), `ruler.py`
  (read-only adapter over `caliper.metrics`), `baselines.py` (curated 7-method
  panel registry), `manifest.py` + `provisional.py` (assumption pins + PROVISIONAL
  stamping), `_paths.py` (sibling-path bootstrap for the read-only deps).
- `ASSUMPTIONS.md` pinning Fashion's ruler (v0.1.0 @ `f078802`, in review at MRM)
  and the substrate (Gauge cohort seed 20260613; OSIPI DRO DOI
  10.5281/zenodo.14605039; Lattice not built).
- `revalidate.py` one-command re-validation; CP1 import/manifest/task tests.

### Notes
- Embedded into the ResearchProj monorepo mirroring Minos: own clean synthetic-only
  history merged with `--allow-unrelated-histories`.
- No reference numbers yet — those are the CP2 deliverable and are PROVISIONAL until
  Fashion's ruler locks.
