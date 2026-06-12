# Proteus

Structure-first mining of the **dark proteome** for **divergent PET-hydrolase
candidates**. Rather than starting from sequence homology to known PETases — which
by construction can only re-find close relatives — Proteus structurally triages
uncharacterized proteins, folds the survivors, then gates on the *catalytic
architecture* of a serine hydrolase (α/β-hydrolase fold, Ser-His-Asp triad,
oxyanion hole, a suitable substrate cleft). Sequence clustering is used only to
**dereplicate** the input corpus, never to homology-gate it: **dereplicate ≠
homology-gate**, so the divergent tail we actually care about stays in play.

## Compute split: LOCAL (M4 Air) vs VAST.AI burst

The pipeline narrows the corpus **locally and cheaply** before spending any GPU
time, so only the S2 shortlist is ever folded.

- **LOCAL — M4 MacBook Air (Apple Silicon, arm64, MPS/CPU, no CUDA):** the cheap
  triage + analysis half. MMseqs2 dereplication (S0), ProstT5 seq→3Di (S1),
  Foldseek fold-class triage (S2), catalytic-geometry and cleft analysis (S4–S5),
  fpocket, and small AutoDock Vina runs.
- **VAST.AI burst (Linux + CUDA — scaffolded in [`vast/`](vast/), NOT installed
  here):** the GPU-heavy half. ESMFold batch folding (S3), Chai-1 refinement, and
  GPU docking (GNINA/DiffDock). Only the **S2 shortlist FASTA** is shipped up;
  folded structures are pulled back for S4–S5. See [`vast/sync.md`](vast/sync.md).

There is **no hard local GPU blocker** precisely because folding is offloaded.

## Pipeline stages (S0–S5)

| Stage | File | Where | What it does |
|---|---|---|---|
| **S0** | `src/proteus/s0_dereplicate.py` | LOCAL | MMseqs2 clustering to **dereplicate** (collapse near-identical seqs). NOT a homology gate. |
| **S1** | `src/proteus/s1_tokenize.py` | LOCAL | ProstT5 translates sequence → Foldseek 3Di alphabet (cheap structural proxy; MPS, cpu fallback). |
| **S2** | `src/proteus/s2_foldclass_triage.py` | LOCAL | Foldseek triage vs the **α/β-hydrolase fold CLASS** (unseeded — fold, not template). Emits the shortlist FASTA. |
| **S3** | `src/proteus/s3_fold.py` | **VAST** | ESMFold batch fold; mean-pLDDT filter; length-chunk long seqs. **Local = dry-run only** (validate FASTA + emit job manifest; never fold on MPS). |
| **S4** | `src/proteus/s4_geometry.py` | LOCAL | Catalytic geometry gate on returned models: Ser-His-Asp triad + oxyanion hole. |
| **S5** | `src/proteus/s5_cleft_filter.py` | LOCAL | Cleft metrics A–E (fpocket), scored **anchored to the positive controls**. |
| P4 | `src/proteus/docking/` | LOCAL (small) | AutoDock Vina docking of PET-mimic ligands into ranked clefts. Large/GPU docking bursts to Vast. |

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
# ship data/interim/s2_shortlist.fasta + s3_job_manifest.json to Vast (vast/sync.md)
```

S3 folds on the Vast burst box. When the models come back, **resume locally**:

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

`controls/references.csv` is the locked control set (positives: IsPETase, LCC_ICCG;
reference: LCC_WT; recovery: GuaPA, MG8; negative: CalB). Fetch the structures:

```bash
python controls/fetch_controls.py     # downloads PDBs -> structures/, writes MANIFEST.json (+sha256)
```

Structure rows are pulled from RCSB with a sha256 manifest. Sequence rows
(GuaPA, MG8) print manual-fetch instructions — see "Manual steps" below.

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
Vast, docking search) must read the seed from there — do not hardcode seeds in
stage code.

## Layout

```
config/proteus.yaml      paths, corpus, thresholds, device (mps|cpu), RANDOM SEED
controls/                locked control set + fetch script + MANIFEST (+sha256)
src/proteus/             pipeline stages S0–S5, docking/, utils/
vast/                    Vast.ai burst fold image (Dockerfile.fold) + sync notes (sync.md)
tests/test_smoke.py      positive-output test per local tool + MPS sanity + S3 dry-run
envlog/                  recon report, env-failures, resolved-version snapshots
data/{raw,interim,processed}/   gitignored corpora/artifacts
```

## Manual steps left to the operator

1. **GuaPA sequence** — archaeal PETase (Acosta et al. 2025); fetch the protein
   sequence from the Marcotte Lab GitHub repo (`marcottelab/GuaPA`), save FASTA to
   `data/raw/`, record the commit hash.
2. **MG8 sequence** — saliva-metagenome PETase (Eiamthong et al., *Angew. Chem.*
   2022); locate the accession/sequence in the paper's SI, deposit FASTA in
   `data/raw/`.
3. **CalB / 1TCA** — confirm PDB **1TCA** is *Candida antarctica* lipase B
   (Uppenberg et al. 1994), the lid-bearing lipase used as the structural negative.
4. **ADFRsuite** — x86-only; install separately (Scripps) under Rosetta 2 if
   needed for receptor prep, or prefer Meeko. Not pip/conda installable on arm64.
5. **Provision a Vast.ai box for S3** — build/push `vast/Dockerfile.fold`, then ship
   the S2 shortlist up and fold (`vast/sync.md`).
6. Run install + `pytest tests/test_smoke.py` **on the M4** to reach a GREEN env and
   generate the real lockfiles.
