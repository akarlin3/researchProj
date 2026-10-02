# Forge — Monte Carlo dose-simulation feasibility benchmark

*No manuscript — a timing/feasibility study, not a paper.* Forge benchmarks Monte Carlo
dose simulation for an MR-Linac (MR-guided linear accelerator) workflow to decide
whether generating a 10,000-case dataset is feasible before a December 2026 deadline.
It uses the TOPAS engine with a 7 MV flattening-filter-free photon beam in a 1.5 T
field (Elekta Unity-class), and validates the **Electron Return Effect** (ERE) — the
magnetic-field-induced dose enhancement at tissue boundaries that is load-bearing for
MR-Linac commissioning.

See the root [README.md](../README.md#forge--monte-carlo-dose-simulation-feasibility-benchmark)
for the full writeup, including the proton-contamination correctness finding.

## Layout
- `forge/` — `geom.py` (phantom + 7 MV FFF photon-spectrum generator), `benchmark.py`, `benchmark_fidelity.py`, `check_ere.py`, `run_case.py`; `decks/` holds TOPAS input templates.
- `cases/` — 12 randomized benchmark cases (voxel phantoms + decks) with `manifest.json`.
- `RESULTS_mc_floor.md` — smoke-test floor (runway verdict: YES, 0.28 extrapolated weeks for 10k cases vs. 22.5-week runway).
- `RESULTS_mc_floor_photon.md` — the proton-contamination discovery and corrected photon-only timing.
- `RESULTS_mc_floor_fidelity.md` — superseded fidelity timing (kept for provenance).

## Reproduce
```bash
python3 forge/benchmark.py       # CPU MC timing benchmark (mock data if TOPAS is absent)
python3 forge/check_ere.py       # Electron Return Effect validation
```
See `RESULTS_mc_floor.md` for TOPAS/Geant4 installation guidance if the benchmark falls
back to mock run data.
