# Forge MC Floor — Photon Re-Time + ERE Gate (Corrected Modality)

Corrects the dataset beam modality to an **Elekta Unity-class 7 MV FFF photon**
source at 1.5 T, re-frames the training-fidelity re-timing as photon-only, and
re-states the ERE gate against a **cited reference magnitude** with a tolerance.

> **Execution status (read first).** This session ran in a Linux x86_64 4-core
> container, **not** the M4 workstation. **OpenTOPAS is not present and cannot be
> installed here** (licensed engine), so the live MC runs (Checkpoints 2–3) were
> **not executed in this environment**. Per the no-fabrication constraint, the
> per-case timing table and the ERE measurement are left **PENDING** — they are
> produced by `forge/benchmark_fidelity.py` when run on the TOPAS host. Everything
> that does *not* require the engine (the modality bug, the beam fix, the grid,
> the reference band + tolerance, the reconciliation logic) is complete and
> verifiable here.

---

## Checkpoint 0 — Recon: where modality is set (load-bearing finding)

Beam modality is chosen in `forge/geom.py`. Before the fix:

```python
# 150 MeV protons or 6 MV photons
particle = rng.choice(["proton", "gamma"])
if particle == "proton":
    energy = rng.uniform(120.0, 180.0)   # MeV
else:
    energy = rng.uniform(4.0, 10.0)      # MeV  (monoenergetic gamma)
```

**Finding:** the 10k generator as configured emitted a **~50/50 proton/photon
mix**, and even the photon arm was a **4–10 MeV monoenergetic gamma** — not a
clinical spectrum. The committed `cases/manifest.json` confirmed the realized
draw was **11 protons + 1 gamma** out of 12. **The dataset itself was wrong, not
just the benchmark.** The prior "training-fidelity" median (0.326 core-h/case)
was dominated by 120–180 MeV protons (~1100 s CPU each) that do not belong in a
photon MR-Linac dataset.

Source/field, as configured (pre-fix): single monoenergetic `Beam`, 0 energy
spread, along +Z, fixed 10×10 cm rectangular cutoff, fixed isocenter. B-field:
`DipoleMagnet`, 1.5 T, direction randomly X **or** Y. Beam travels along Z, so a
field in X/Y is **perpendicular to the beam axis (B ⟂ beam)** — the B-field
orientation was already correct for Unity; only the particle/spectrum was wrong.

---

## Checkpoint 1 — Beam fixed to a Unity-class 7 MV FFF photon source

`forge/geom.py` now emits **photons only**, with a continuous 7 MV FFF
bremsstrahlung energy spectrum, and randomizes field size + isocenter so the
geometry distribution genuinely varies. B ⟂ beam confirmed (X/Y transverse,
Bz = 0, beam along Z).

**Spectrum (flagged APPROXIMATE — synthesized, not a measured Unity phase space).**
Shape = Kramers thick-target bremsstrahlung `(Emax−E)/E` with an empirical
low-energy hardening factor `E^0.65` (no flattening filter present on Unity),
7 MV endpoint, 0.25 MeV low cutoff, 40 bins. The hardening exponent was tuned so
the fluence-mean energy is **2.09 MeV**, matching the published Unity 7 MV FFF
mean photon energy of **~2.11 MeV** at isocenter from BEAMnrc/EGSnrc commissioning
literature. Swap in a commissioned phase space when available.

Generated TOPAS source block (excerpt):

```
s:So/Beam/BeamParticle = "gamma"
s:So/Beam/BeamEnergySpectrumType = "Continuous"
dv:So/Beam/BeamEnergySpectrumValues = 40 0.3344 ... 6.9156 MeV
uv:So/Beam/BeamEnergySpectrumWeights = 40 8.78e-02 ... 3.85e-04
d:So/Beam/BeamPositionCutoffX = <randomized 1.5–11.0 cm half-width> cm   # 3×3 … 22×22 fields
d:So/Beam/TransX/Y = <randomized isocenter offset, ±3 cm>
```

**Manifest diff (the dataset is now photon):**

| | particle counts | photon spectrum | field size | isocenter |
| :--- | :--- | :--- | :--- | :--- |
| Before | proton ×11, gamma ×1 | 4–10 MeV mono-gamma / 120–180 MeV proton | fixed 10×10 cm | fixed (0,0) |
| After  | **gamma ×12** | **7 MV FFF, mean 2.09 MeV** | **randomized 3.3–20.7 cm** | **randomized ±3 cm** |

This fix is applied to the **10k generator** (`geom.py`), not just the benchmark
cases. Generating the full 10k remains out of scope.

---

## Checkpoint 2 — Photon re-time at training fidelity

**Voxel grid (REPORTED — required this run):**
- **Dimensions:** 60 × 60 × 60 voxels
- **Resolution:** 2.5 × 2.5 × 2.5 mm
- **Physical extent:** 150 × 150 × 150 mm
- 2.5 mm is **coarse** for penumbra/ERE structure; a 1% σ target is comparatively
  cheap on this grid. The grid was **not coarsened for speed** — it is the
  documented dataset grid. A finer grid (≈1–2 mm) would raise the per-case cost.

**Per-case timing table — PENDING (engine not available in this container).**
`forge/benchmark_fidelity.py` draws N=8 photon cases from the randomized
distribution, runs each to ~1% high-dose σ (scaling histories), and records
primaries, achieved σ, wall-clock, cores, and **true CPU core-hours**, with
median/IQR/min/max. It will fill this table when run on the TOPAS host.

| Case | Field (cm) | Grid | Primaries | Achieved σ | Target prim. for 1% | CPU (s) | Core-h |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| _pending — run on TOPAS host_ | | | | | | | |

> Caveat for the estimate: the only prior photon-type datapoint (case 1, a
> 4 MeV **monoenergetic** gamma) projected to **0.62 core-h** — *higher* than the
> contaminated proton median. That hints photons may be **comparable-to-costly**,
> not automatically cheap; the corrected 7 MV FFF figure must be measured.

---

## Checkpoint 3 — ERE gate vs reference (photon, 1.5 T)

Water–air–water interface, 1.5 T transverse, B-on vs B-off, photon beam. The gate
compares the **peak distal-interface enhancement** to a cited reference band.

- **Reference ERE magnitude (cited):** ~**15.4%** (10×10 cm) to ~**17.9%**
  (22×22 cm) exit-dose enhancement at 1.5 T transverse, from Geant4 Monte Carlo of
  the magnetic-field effect on dose in low-density regions. This is a *literature
  magnitude*, **not** an in-repo validated reference profile.
- **Tolerance:** ± 3 percentage points.
- **In-repo reference profile:** none.
- **Measured enhancement:** PENDING (engine not available here).
- **Gate verdict:** **UNCALIBRATED.** With no in-repo commissioned reference
  profile — and no execution in this environment — the gate cannot be certified
  PASS, by the gate's own rule (PASS requires an in-repo reference comparison).

**Diagnosis of the prior +0.32% result:** that effect is implausibly small for
1.5 T ERE (literature: tens of percent). Likely causes in the prior setup: (a) a
**monoenergetic 4 MeV gamma** rather than a clinical spectrum, (b) reading a
**single voxel** at the interface instead of scanning for the localized peak, and
(c) coarse 2.5 mm voxels smearing the sub-mm ERE peak. The updated gate now scans
a small window around the distal interface for the peak; the spectrum and B⟂beam
geometry are corrected. A <1% result after these fixes should be treated as a
geometry/field bug, not a PASS.

---

## Checkpoint 4 — Honest extrapolation

Formula: `weeks = (10000 × core_h_per_case) / (10 × 168)`.

- **Photon median core-h/case:** PENDING (requires the Checkpoint 2 run).
- **Weeks at median + IQR/min–max:** PENDING (computed from the above).

**Reconciliation with prior numbers:**
- **0.0465 core-h/case** — smoke-test floor (50k histories, no σ analysis). Not a
  fidelity number.
- **0.326 core-h/case** — prior "training-fidelity" median, **proton-contaminated**
  (7/8 timed cases were 120–180 MeV protons at ~1100 s CPU each). **Invalid for a
  photon machine.**
- **Corrected photon figure:** to be measured on the TOPAS host with the fixed
  generator. It will differ from 0.326 because the proton CPU cost is removed, but
  the 7 MV FFF photon cost is **not yet known** and the one prior gamma datapoint
  cautions it could be non-trivial.

### One-line verdict on the 10k duration
**Undetermined in this environment** — the modality bug is fixed and the dataset
is now photon, but the corrected photon core-hours (and therefore the 10k
duration and ERE PASS/FAIL) require execution on the OpenTOPAS host; the prior
0.326-based duration is invalid because it was proton-contaminated.

---

### Reproduce on the TOPAS host
```bash
python3 forge/geom.py                 # regenerate the photon dataset/decks
python3 forge/benchmark_fidelity.py   # N=8 photon re-time + ERE gate -> overwrites this file
```
