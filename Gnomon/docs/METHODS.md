# Methods (the complete write-up Fashion was flagged for lacking)

> This is the methods section Gnomon exists to produce — **complete from line one**,
> covering every item Fashion's review flagged. Sections marked `[CP2]`/`[CP3]` are
> filled by the rebuild run; the *specifications* are fixed now so nothing is
> reverse-engineered after seeing results. Numbers, when added, trace to
> `results/reproduction.json`.

## 1. Datasets (every dataset named with a stable ID)

### 1.1 Synthetic substrate — Lattice DRO
- Source: `lattice.make_cohort` (sibling, MIT, synthetic-only), imported read-only.
- Ground truth `(D, Dstar, f)` drawn from Lattice's documented physiological priors;
  seed = `manifest.MASTER_SEED` (20260620), fully reproducible from the seed.
- **Headline 9-cell design:** 3 ground-truth conditions × SNR ∈ {10, 20, 40} × 200
  noise realizations (Rician), matching Fashion's headline set (README.md:58-59).

### 1.2 Open real data — OSIPI abdomen
- **OSIPI TF2.4**, Zenodo record **14605039** (`OSIPI_TF24_data_phantoms.zip`),
  CC-licensed; the in-vivo **abdomen** acquisition (`Data/abdomen.*`).
- Fetched on demand into a gitignored `download/` with a provenance manifest
  (record id, SHA-256 `2a53054d…b3e`). Never redistributed in-tree. `gnomon/osipi.py`.
- **Acquisition (read from the archive):** 144×144×21 volume, **104 measurements**
  over 12 unique b-values {0,10,20,30,40,50,75,100,150,250,400,600} s/mm² (Philips 3T).
- **ROI selection (documented):** the archive's `Data/mask_abdomen_homogeneous.nii.gz`
  homogeneous-tissue mask = **1932 voxels**; each voxel b0-normalized (÷ mean b=0).
  Fashion's "1618 high-SNR ROI voxels" is this ROI under an *unstated* SNR cut — a
  completeness gap. Gnomon reports the D\* railing rate on the **full homogeneous ROI**
  and on a b0-SNR>25 subset (b0-SNR = mean/std over the 15 b=0 repeats), so the
  load-bearing number (the rate) carries its selection sensitivity rather than being
  pinned to an unreproducible voxel count.

## 2. Acquisition (b-value schemes, explicit)

- **clinical-sparse (8 b):** `(0, 10, 20, 40, 80, 200, 400, 800)` s/mm² — Gnomon's
  documented clean-room choice for the *synthetic* design (Fashion's prose did not
  list its 8-b values).
- **dense (16 b):** `(0,10,20,30,50,75,100,150,200,300,400,500,600,700,800,1000)`
  s/mm² — quoted verbatim from Fashion (REVIEWER_RESPONSE.md:126-127).
- **real OSIPI abdomen:** the native 104-point / 12-unique-b scheme above (used as-is
  for the T1 railing fit).

## 3. Forward model

IVIM bi-exponential `S(b)/S0 = (1-f)·exp(-b·D) + f·exp(-b·(D+Dstar))` (Le Bihan
1988), reimplemented from physics (`gnomon/forward.py`). Estimators fit in **scaled**
params `(S0, D3, f, Ds3)` with `D = D3·1e-3`, `Dstar = Ds3·1e-3` so every fit variable
is O(1). Analytic Jacobian (verified vs finite differences). Self-consistency gates:
clean-signal NLLS round-trip recovers truth; continuity `f=0 ⇒ mono-exponential`.

## 4. Estimators & fitting (full detail)

### 4.1 NLLS + boundary-railing
- Box-constrained four-parameter fit (`scipy.optimize.least_squares`, trust-region
  reflective, analytic Jacobian, `max_nfev=400`). **Box** (scaled): `S0∈[0.5,1.5]`,
  `D3∈[0.2,3.0]`, `f∈[0,0.5]`, `Ds3∈[3,150]` (D/Dstar bounds follow Fashion's stated
  NPE prior range). **Init** `(1.0, 1.0, 0.1, 20.0)`.
- Covariance: `σ²·pinv(JᵀJ)` (known σ=1/SNR for the synthetic Laplace; residual-based
  for the real-data fit), each SD capped at the box span so a railed/unidentified D\*
  yields a finite-but-pathological interval.
- **Railing:** parameter railed iff `|x̂ − bound|/(upper−lower) < rail_tol`;
  `rail_tol = 1e-3` primary, `1e-2` sensitivity. D\* railing **rate** is target T1.

### 4.2 Laplace posterior
- Gaussian at the MAP (NLLS fit) with the CRLB covariance (known σ); symmetric SD
  interval → T3a.

### 4.3 MCMC posterior
- Clean-room per-voxel random-walk Metropolis, **vectorized across all voxels**, over a
  Gaussian likelihood (known σ=1/SNR) with a uniform-box prior (the NLLS box).
- Proposal: isotropic Gaussian RW, per-voxel per-param step = `0.6 ×` the CRLB SD;
  init at the NLLS MAP. **burn 1500, keep 2000, thin 2**, seed `MASTER_SEED+3`;
  acceptance ~0.25 (reported per run). From one chain: **SD interval** (T3b) and
  **2.5/97.5 quantile interval** (T3c).

### 4.4 MAF amortized posterior (NPE)
- Conditional autoregressive flow (5 affine AR layers, per-dim MLP conditioners,
  hidden 64, dim order reversed between layers; standard-normal base) over scaled
  `(D3,f,Ds3)` standardized by the prior; context = standardized signal. **Sim budget
  80 000** `(θ,x)` pairs, θ~uniform-box, SNR~U(10,40); Adam lr 1e-3, **40 epochs**,
  batch 512, seed `MASTER_SEED+5`. Inference: 1500 posterior draws/voxel → quantiles.
  Optional `[flow]` extra (torch). No sbi/nflows — written from scratch.

## 5. CRLB assumption (the reviewer-flagged item) `[CP2]`

The Laplace/asymptotic covariance is a **Gaussian** approximation to the posterior.
It is **weakest exactly where D\* is most skewed** (low SNR), which is where the
central under-coverage claim sits. Stated as a limitation; the direction of the bias
(conservative — it *understates* how overconfident the Gaussian interval is) is noted
so no claim rests on the approximation being tight.

## 6. Calibration ruler (re-derived independently) `[CP2]`

Coverage, ECE (mean |empirical − nominal| over a level grid), sharpness (mean
interval width), pinball loss, interval score — from published definitions, numpy
only, no Caliper import. Plus a D\*-tercile conditional-coverage probe.

## 7. Uncertainty on the numbers `[CP2]`

Bootstrap CIs (percentile/BCa, seeded from `manifest.BOOTSTRAP`) on every
load-bearing number; directional gaps (T4) pass only if the CI excludes 0.

## 8. Reproduction & verdict `[CP3]`

`gnomon/reproduce.py` runs the rebuild, compares to the frozen manifest, and writes
`results/reproduction.json` + the verdict (REPRODUCES / DOES NOT REPRODUCE). One
command: `bash reproduce.sh`.

## 9. Claims scope

Claims are held to what the rebuild supports. Where Fashion was read as overextending
(e.g. single-subject in-vivo generalization), Gnomon scopes the corresponding
statement to its evidence and records the limitation here.
