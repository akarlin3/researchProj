# Ferry — grounding Matrix on a real RT dataset (honest scope)

Ferry is a **real-data substrate adapter** for Matrix's closed loop. It swaps Matrix's
*synthetic* anatomy/dose twin for **real anatomy + real dose geometry** from a public
radiotherapy dataset, to pre-empt the *"shown only on a pure synthetic twin"* objection for a
pre-PhD Matrix publication. It is an **interface-swap only**: a `GroundedTwin` drops into the
existing `run_iteration` engine, so **`loop.py` is byte-unchanged**.

> ## Scope ceiling — read this first
>
> Ferry grounds **anatomy** and **dose geometry** on real data. **Perfusion / IVIM stays
> synthetic** — there is no scanner, so there is **no real diffusion data**. A grounded result
> means *"the loop closes on **real geometry**"* and **nothing more**. It is **not** a
> real-IVIM result, **not** a validated clinical loop, and makes **no clinical claim.**

## The dataset (CP0 selection: clinical-leaning venue)

| | |
|---|---|
| dataset | **TCIA Pancreatic-CT-CBCT-SEG** (Version 2, 2022-08-23) |
| DOI | `10.7937/TCIA.ESHQ-4D90` |
| license | **CC BY 4.0** (verified; gate-free, attribution only) — see [`matrix/ferry/LICENSE_DATASET.md`](matrix/ferry/LICENSE_DATASET.md) |
| what is used | one patient's **RTSTRUCT** (contours) + **RTDOSE** (delivered 3-D dose grid) |
| clean IP | **no data blobs committed**; loaded by script via the TCIA NBIA API into a git-ignored cache |

The **venue sets the realism bar.** This is a *clinical-leaning* grounding (real patient
pancreatic anatomy + a real clinical dose plan), chosen because it matches the twin's
abdominal/pancreatic target and its RTDOSE grid maps directly onto the loop's per-voxel dose.
The cost of that choice is a **higher realism bar**: a clinical reviewer will (correctly) note
that perfusion is still synthetic. A methods/DRO venue (e.g. TROTS/CORT liver, dose-influence
matrices) would lower that bar but with weaker anatomical relevance; the substrate adapter is
venue-agnostic — only the loader's dataset changes.

## What becomes real vs stays synthetic

| loop quantity | source | notes |
|---|---|---|
| `labels` (anatomy) | **REAL** | RTSTRUCT: target VOI → TUMOR; small bowel + stomach/duodenum → OAR; rest → NORMAL |
| `dose` (dose geometry) | **REAL** | RTDOSE 3-D grid, rescaled into `[dose_min, dose_max]` preserving spatial geometry |
| `D, D*, f` (IVIM) | **SYNTHETIC** | seeded priors per label — **no scanner, no real perfusion** |
| `scan()`, `snr_map`, `lowsnr` | **SYNTHETIC** | the same seeded acquisition + degradation model as the twin |
| `highdstar` | **SYNTHETIC** | the same seeded high-D* identifiability sub-region |

The synthetic layers use the **identical mechanism and RNG order** as `Twin.build`, so the
*only* differences between the synthetic twin and the grounded twin are the two real fields.
That isolation is what lets CP2 attribute behaviour changes to *real geometry* and nothing else.

## What a reviewer CAN and CANNOT conclude

**CAN conclude (SOLID):** Matrix's closed loop **closes on real geometry**. On a real
pancreatic target and a real delivered-dose grid the loop runs end-to-end and reproducibly;
the trust gate holds the *action* rate on untrustworthy voxels at 0 while trustworthy voxels
stay free; trusted tumour perfusion falls under treatment (bootstrap CI excludes 0) and the
treatment decision converges (`n_treat → 0`). The harness is not an artefact of the synthetic
twin's tidy geometry. (Gate: `verify_ferry_cp2.py`; numbers in
[`results/RESULTS_FERRY_CP2.md`](results/RESULTS_FERRY_CP2.md).)

**CANNOT conclude:** anything about **real perfusion / IVIM**. The `(D, D*, f)` ground truth,
the scan, and the SNR map are synthetic. Ferry does **not** validate the IVIM estimator, the
calibration ruler, or the trust gate against real diffusion data — there is none. Every
calibration / trust / convergence number still rides on the synthetic perfusion layer and the
placeholder components (NOT-Fashion / NOT-Minos / NOT-Forge).

## The residual real-diffusion gap (only a scanner closes it)

Ferry closes the *geometry* gap (real anatomy + real dose). It does **not** close the
*diffusion* gap: turning any of Matrix's perfusion-dependent readings into a clinical claim
requires **real IVIM acquisition** — a scanner producing real multi-b diffusion data on the
same anatomy. That is Keystone's real-time / offline modes (scanner + Forge), explicitly out
of scope here. Ferry's honest contribution is exactly one rung of the ladder: *the loop
mechanics survive real geometry*, leaving only real diffusion to a future scanner study.

## Findings under real geometry (flagged, not failures)

The grounded run surfaced behaviour the synthetic twin cannot show (full detail in
`results/RESULTS_FERRY_CP2.md`):

- **F1 — action-suppression is not outcome-protection.** On the synthetic twin (and on real
  anatomy with a flat-baseline dose), "holding" an untrusted voxel leaves its perfusion
  untouched (held drop `0.000`). On **real dose geometry** the held untrusted tumour
  perfusion still drops (`≈0.148`, CI excludes 0): the real *delivered* dose already
  devascularises it. The trust gate suppresses new **actions**, not dose already delivered —
  a substantive insight for adaptive loops, visible only once dose geometry is real.
- **F2 — a NOT-Forge placeholder artefact.** The strict "TREAT ⇒ dose strictly increases"
  warrant breaks on a real non-uniform prescription (a hot voxel TREATed toward the
  placeholder's *absolute* boost target decreases) and on re-TREAT (Δ=0). This is a property
  of the analytic NOT-Forge placeholder — present even in the synthetic baseline at a finer
  grid — and Forge's real geometry-aware engine resolves it. It is not a loop failure.
- **F3 — scale/shape.** Real anatomy is larger and irregular; both trusted and untrusted
  tumour populations remain non-empty and the loop's qualitative behaviour is preserved.

## Run it

```bash
PROT=/opt/homebrew/Caskroom/miniforge/base/envs/proteus/bin/python   # the proteus env
# CP1 — the drop-in proof (loop.py byte-unchanged; contract; reproducible; end-to-end)
PYTHONPATH=Matrix $PROT Matrix/verify_ferry_cp1.py
# CP2 — the grounded closed-loop run + side-by-side vs the synthetic baseline (needs network)
PYTHONPATH=Matrix $PROT Matrix/verify_ferry_cp2.py
# tests (offline mechanism + one network-gated real-data test)
PYTHONPATH=Matrix $PROT -m pytest Matrix/tests/test_ferry.py -q
```

The loader downloads ~30 MB (one patient's RTSTRUCT + RTDOSE) from TCIA/NBIA on first run and
caches the derived grids locally (git-ignored). No patient data is ever committed.
