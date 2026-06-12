# Proteus

Structure-first mining of the **dark proteome** for **divergent PET-hydrolase
candidates**. Rather than starting from sequence homology to known PETases — which
by construction can only re-find close relatives — Proteus structurally triages
uncharacterized proteins, folds the survivors, then gates on the *catalytic
architecture* of a serine hydrolase (α/β-hydrolase fold, Ser-His-Asp triad,
oxyanion hole, a suitable substrate cleft). Sequence clustering is used only to
**dereplicate** the input corpus, never to homology-gate it: **dereplicate ≠
homology-gate**, so the divergent tail we actually care about stays in play.

## Compute split: LOCAL (M4 Air) vs GCE burst

The pipeline narrows the corpus **locally and cheaply** before spending any GPU
time, so only the S2 shortlist is ever folded.

- **LOCAL — M4 MacBook Air (Apple Silicon, arm64, MPS/CPU, no CUDA):** the cheap
  triage + analysis half. MMseqs2 dereplication (S0), ProstT5 seq→3Di (S1),
  Foldseek fold-class triage (S2), catalytic-geometry and cleft analysis (S4–S5),
  fpocket, and small AutoDock Vina runs.
- **GCE burst (scaffolded in [`gce/`](gce/), NOT installed here):** the heavy fold
  half. ESMFold batch folding (S3). This project has **no GPU quota**, so the fold
  runs on **CPU** (a C4 VM) — slow + RAM-heavy, but the narrowed shortlist is small;
  flip `compute.gce_burst.accelerator` to move to GPU once quota is granted. Only the
  **S2 shortlist FASTA** is shipped up; folded structures are pulled back for S4–S5.
  See [`gce/sync.md`](gce/sync.md).

Folding is offloaded, so there is **no hard local compute blocker**.

## Pipeline stages (S0–S5)

| Stage | File | Where | What it does |
|---|---|---|---|
| **S0** | `src/proteus/s0_dereplicate.py` | LOCAL | MMseqs2 clustering to **dereplicate** (collapse near-identical seqs). NOT a homology gate. |
| **S1** | `src/proteus/s1_tokenize.py` | LOCAL | ProstT5 translates sequence → Foldseek 3Di alphabet (cheap structural proxy; MPS, cpu fallback). |
| **S2** | `src/proteus/s2_foldclass_triage.py` | LOCAL | Foldseek triage vs the **α/β-hydrolase fold CLASS** (unseeded — fold, not template). Emits the shortlist FASTA. |
| **S3** | `src/proteus/s3_fold.py` (local dry-run) + `gce/run_fold.py` (on-box) | **GCE** | ESMFold batch fold; mean-pLDDT filter; length-chunk long seqs. **Local = dry-run only** (validate FASTA + emit job manifest; never fold on MPS). The on-box runner folds on **CPU** (no GPU quota), resumable + pLDDT-gated. |
| **S4** | `src/proteus/s4_geometry.py` | LOCAL | Catalytic geometry gate on returned models: Ser-His-Asp triad + oxyanion hole. |
| **S5** | `src/proteus/s5_cleft_filter.py` | LOCAL | Cleft metrics A–E (fpocket), scored **anchored to the positive controls**. |
| P4 | `src/proteus/docking/` | LOCAL (small) | AutoDock Vina docking of the PET-mimic BHET into the catalytic cleft (box on the S4 Ser OG). Large/GPU docking bursts to GCE. |

### Ingestion & orchestration

| Module | What it does |
|---|---|
| `proteus.fetch_corpus` | Resolve `corpus.sources` (UniProt queries / URLs) into `data/raw` shards + a provenance manifest. |
| `proteus.corpus` | Assemble + length-filter the raw shards into one corpus FASTA (the S0 input); sanitises ids so MMseqs2/Foldseek and the parsers agree. |
| `proteus.pipeline` | One command: corpus → S0 → S1 → S2 → S3 manifest, with the narrowing-funnel report. |
| `proteus.launch` | Drive the GCE burst (validate manifest → stage-up → create → fold → stage-down → delete). **Dry-run by default.** |
| `proteus.screen` | Resume after folding: S4 geometry → S5 cleft → control-anchored score → ranked PETase-like hits. |
| `proteus.calibrate` | Calibrate S4+S5 on the controls (separation verdict, operating point) + held-out divergent-positive recovery. |

## Run the local narrowing pipeline

The front door (`proteus.corpus`) assembles + length-filters the raw FASTA shards
in `data/raw` into one corpus; `proteus.pipeline` then chains **corpus → S0 → S1 →
S2 → S3 job manifest** in one command and prints the narrowing funnel:

```bash
# (optional) fetch the configured corpus.sources (UniProt queries / URLs) into data/raw:
PYTHONPATH=src python -m proteus.fetch_corpus
# ...or drop dark-proteome FASTA(.gz) shards in data/raw/ (corpus.fasta_glob) directly, then:
PYTHONPATH=src python -m proteus.pipeline --out data/interim
#   corpus N (length-filtered) → S0 reps → S2 shortlist → S3 manifest
# ship data/interim/s2_shortlist.fasta + s3_job_manifest.json to GCE (gce/sync.md)
```

S3 folds on the GCE burst box. When the models come back, **resume locally**:

```bash
PYTHONPATH=src python -m proteus.screen \
  --folded structures/folded --struct-dir structures \
  --out data/processed/s4s5_candidates        # S4 geometry + S5 cleft, control-anchored
PYTHONPATH=src python -m proteus.docking \
  --candidates data/processed/s4s5_candidates.json --models structures/folded \
  --out data/processed/docking                # Vina dock BHET into the PETase-like hits
```

(S2 needs pre-built reference DBs — set `s2_foldclass_triage.references[].db`.)

## Controls

`controls/references.csv` is the locked control set:

- **Positives** (anchor the cleft scoring): IsPETase (6EQE), LCC_WT (4EB0); LCC_ICCG
  (6THS) is the S165A inactivation trap.
- **Negatives** (share the fold/triad but not PET activity): CalB (1TCA), AChE (1EA5),
  CRL (1CRL), Est2 (1EVQ).
- **Recovery** — held-out **divergent positives** that test generalization without
  re-anchoring: PET46 (8B4U, archaeal), Cut190 (4WFI), TfCut2 (4CG1). GuaPA and MG8
  remain **unresolved** (no reachable sequence/structure) and need a GCE fold to
  enter calibration — the three structures above stand in for them.

```bash
python controls/fetch_controls.py     # downloads PDBs -> structures/, writes MANIFEST.json (+sha256)
```

Structure rows are pulled from RCSB with a sha256 manifest. The PET-mimic docking
ligand is committed at `controls/ligands/bhet.pdbqt` (BHET, prepared with Open Babel).

## Validation

`envlog/validation-run.md` records the first end-to-end run on a known-answer corpus:
the screen separates PET-hydrolases from non-PET serine hydrolases (8/15 hits, margin
1.4), a single widened operating point recovers all three held-out divergent positives
(bacterial + actinomycete + archaeal) at precision 1.0, and BHET docking is shown to be
confirmatory, not discriminating. It is the face-validity evidence for
"dereplicate ≠ homology-gate" before pointing the pipeline at a real dark corpus.

## Reproduce the environment (on the M4)

The local env targets **osx-arm64 (Apple Silicon), MPS/CPU, no CUDA**, built with
Miniforge:

```bash
CONDA_SUBDIR=osx-arm64 mamba env create -f environment.yml
conda activate proteus
python -c "import torch; print('mps', torch.backends.mps.is_available())"   # expect True
python controls/fetch_controls.py
pytest tests/test_smoke.py -v          # every LOCAL tool must show positive output → GREEN
# then capture the real locks:
pip freeze        > requirements-lock.txt
conda env export  > envlog/conda-env-resolved.yml
```

Resolution order per tool: native osx-arm64 conda → Homebrew (arm64) → x86 via
`CONDA_SUBDIR=osx-64` under Rosetta 2 (last resort).

Pinned/resolved versions live in the lockfiles — treat them as **generated, not
hand-edited**:

- `requirements-lock.txt` — `pip freeze` of the resolved env.
- `envlog/conda-env-resolved.yml` — `conda env export` of the resolved env.
- `envlog/recon-report.md` — host recon (OS/arch/RAM/MPS) + GO/WARNING verdict.
- `envlog/env-failures.md` — per-tool resolution source + any deviations.

> ⚠️ **This scaffold was generated in a Linux x86_64 cloud container, not on the
> M4.** MPS and the osx-arm64 conda solve therefore could not be exercised here, so
> the committed lockfiles are **placeholders pending generation on the Mac** — run
> the commands above to produce the real freeze and a GREEN env. The host-agnostic
> piece that *was* validated: `s3_fold.py --dry-run` emits a valid job manifest, and
> the smoke suite collects/skips cleanly. Details in `envlog/env-failures.md`.

## Reproducibility / random seed

The global `random_seed` lives in [`config/proteus.yaml`](config/proteus.yaml)
(`random_seed: 1729`) alongside paths, device (`mps|cpu`), corpus selection, and
per-stage thresholds. **Every stochastic step** (S0 MMseqs2, S3 ESMFold/torch on
GCE, docking search) must read the seed from there — do not hardcode seeds in
stage code.

## Tests & CI

The suite is **positive-output** (asserts real artifacts, not just "no exception") and
**tool-gated** (a test skips cleanly if its tool/weights/structures are absent, so it
never falsely fails). Install the local toolchain — the **same script** CI and the web
SessionStart hook use — then run ruff + pytest:

```bash
bash scripts/setup_tools.sh                 # MMseqs2, Foldseek, fpocket, Open Babel,
                                            # Python deps + ruff + vina, control structures.
PROTEUS_FETCH_WEIGHTS=1 bash scripts/setup_tools.sh   # also fetch the ~2.4 GB ProstT5 weights
export PATH="$HOME/.proteus-tools/bin:$PATH"
ruff check src/ tests/
python -m pytest -q
```

- **CI** ([`.github/workflows/ci.yml`](.github/workflows/ci.yml)) runs ruff + pytest on
  every PR/push, reusing `scripts/setup_tools.sh` (no weights → S1/S2/pipeline skip).
- **Claude Code on the web**: a `SessionStart` hook (`.claude/`) auto-installs the full
  toolchain (incl. weights) so tool-gated tests run in-session.

## Layout

```
config/proteus.yaml      paths, corpus (+ sources), thresholds, device (mps|cpu), RANDOM SEED
controls/                locked control set + fetch script + MANIFEST (+sha256); ligands/bhet.pdbqt
src/proteus/             stages s0–s5; corpus, fetch_corpus, pipeline, screen, launch, calibrate; docking/, utils/
gce/                    burst fold image (Dockerfile.fold), on-box run_fold.py, sync notes (sync.md)
tests/                   positive-output tests per stage + smoke (per-tool) + end-to-end funnel
scripts/setup_tools.sh   shared toolchain installer (CI + SessionStart hook)
.github/workflows/ci.yml ruff + pytest on PRs
.claude/                 Claude-Code-on-the-web SessionStart hook
envlog/                  recon report, env-failures, resolved-version snapshots
data/{raw,interim,processed}/   gitignored corpora/artifacts
```

## Manual steps left to the operator

1. **GuaPA / MG8 sequences (unresolved)** — GuaPA (archaeal PETase, Acosta et al.
   2025) repo `marcottelab/GuaPA` is unreachable; MG8 (saliva-metagenome PETase,
   Eiamthong et al., *Angew. Chem.* 2022) has no PDB and its SI accession is not
   located. If you obtain either sequence, drop the FASTA in `data/raw/`, fold it on
   GCE, and add the structure to calibration. Until then the divergent-positive
   structures PET46/Cut190/TfCut2 stand in (see Controls).
2. **ADFRsuite** — x86-only; install separately (Scripps) under Rosetta 2 if needed
   for receptor prep, or prefer Open Babel / Meeko. Not pip/conda installable on arm64.
3. **Provision a GCE box for S3** — build/push `gce/Dockerfile.fold` (pin the
   ESMFold toolchain on the first build — CPU torch, no CUDA), then drive the burst with
   `proteus.launch` (or the manual steps in `gce/sync.md`).
4. Run install + `pytest` **on the M4** to reach a GREEN env and generate the real
   lockfiles (`requirements-lock.txt`, `envlog/conda-env-resolved.yml`).
