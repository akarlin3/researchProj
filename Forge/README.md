# Forge — Monte Carlo dose-simulation feasibility benchmark

*No manuscript — a timing/feasibility study, not a paper.* Forge benchmarks Monte
Carlo dose simulation for an MR-Linac (MR-guided linear accelerator) workflow to
decide whether generating a 10,000-case dataset is feasible before a December 2026
deadline. It uses the TOPAS engine with a 7 MV flattening-filter-free (FFF) photon
beam in a 1.5 T field (Elekta Unity-class), and validates the **Electron Return
Effect** (ERE) — the magnetic-field-induced dose enhancement at tissue boundaries
that is load-bearing for MR-Linac commissioning.

## Status — modality bug found and fixed; ERE gate PENDING on the TOPAS host

The original 10k-case generator drew a ~50/50 proton/photon mix (11 protons + 1
gamma in the committed 12-case sample), and even its photon arm was a
non-clinical 4–10 MeV monoenergetic gamma rather than a 7 MV FFF spectrum — the
dataset itself was wrong, not just the benchmark. `forge/geom.py` now emits
photons only, with a continuous 7 MV FFF bremsstrahlung spectrum (fluence-mean
2.09 MeV, matching the published Unity ~2.11 MeV commissioning figure) and
randomized field size/isocenter. This fix applies to the full 10k generator, not
just the 12-case benchmark sample.

The per-case timing re-time and the ERE gate against a cited reference band
(~15.4–17.9% exit-dose enhancement at 1.5 T, ±3 pp tolerance) both require
`OpenTOPAS`, a licensed engine not installable in this container; they are
**PENDING** execution on the TOPAS host (`forge/benchmark_fidelity.py`). The
prior "training-fidelity" median of 0.326 core-h/case is **invalid** — it was
proton-contaminated (7/8 timed cases were 120–180 MeV protons) — and the prior
ERE reading of +0.32% is implausibly small for 1.5 T (likely a monoenergetic
beam + single-voxel readout + coarse-voxel artifact), so the gate is
**UNCALIBRATED**, not PASS/FAIL, until re-measured on the corrected photon setup.

## Layout

- `forge/` — `geom.py` (phantom + 7 MV FFF photon-spectrum generator), `benchmark.py`, `benchmark_fidelity.py`, `check_ere.py`, `run_case.py`; `decks/` holds TOPAS input templates.
- `cases/` — 12 randomized benchmark cases (voxel phantoms + decks) with `manifest.json`.
- `RESULTS_mc_floor.md` — smoke-test floor (50k histories, no σ analysis; not a fidelity number).
- `RESULTS_mc_floor_fidelity.md` — the (superseded, proton-contaminated) training-fidelity timing.
- `RESULTS_mc_floor_photon.md` — the proton-contamination discovery, the photon-beam fix, and the current PENDING re-time/ERE gate.

## Reproduce on the TOPAS host

```bash
python3 forge/geom.py                 # regenerate the photon dataset/decks
python3 forge/benchmark_fidelity.py   # N=8 photon re-time + ERE gate -> RESULTS_mc_floor_photon.md
```
