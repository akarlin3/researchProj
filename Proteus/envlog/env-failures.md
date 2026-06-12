# Proteus — env resolution status, failures & deviations

Per Checkpoint 2: "Any tool that won't resolve: record in `envlog/` with the
reason and continue — don't abort the env over one package."

**Target host:** M4 MacBook Air (osx-arm64, MPS/CPU, no CUDA).
**Where this scaffold was generated:** Linux x86_64 cloud container (NOT the Mac).

## The one structural fact that gates everything below

The install + verify step (Checkpoint 2) and the positive-output smoke suite
(Checkpoint 3) **must run on the M4** to be meaningful. They could not be executed
here because this container is Linux x86_64 with no conda/mamba and no Apple
Silicon / MPS. So the lockfiles (`requirements-lock.txt`,
`envlog/conda-env-resolved.yml`) are **placeholders pending generation on the Mac**
— not a real freeze. Fabricating osx-arm64 pins from a Linux host would be wrong
(it would pin the wrong-platform wheels), so we explicitly did not.

## Per-tool resolution plan & expected source (resolve on the Mac, then lock)

| Component | Expected source on osx-arm64 | Status / note |
|---|---|---|
| **PyTorch (MPS)** | conda-forge osx-arm64 (`pytorch`), else pip arm64 wheel | resolve on Mac; assert `torch.backends.mps.is_available()`. CPU/MPS only — **no CUDA**. |
| **transformers / sentencepiece** | pip | ProstT5 tokenizer. Some `transformers` releases mis-route ProstT5's SentencePiece model through an incompatible converter; smoke test falls back to loading `spiece.model` via `sentencepiece` directly. Pin a transformers whose T5 converter handles ProstT5, or rely on the fallback. |
| **Foldseek** | bioconda osx-arm64 → **Homebrew** (`brew install foldseek`) | bioconda arm coverage varies; Homebrew is the reliable arm64 fallback. Record which resolved. |
| **MMseqs2** | bioconda osx-arm64 → **Homebrew** (`brew install mmseqs2`) | same fallback chain as Foldseek. |
| **fpocket** | bioconda → **Homebrew** (`brew install fpocket`) → source build | most likely Homebrew or source on arm64. |
| **AutoDock Vina** | conda-forge osx-arm64 (`vina`) + `meeko` (pip) | Meeko preferred for ligand/receptor prep. |
| **ADFRsuite** | **NOT installed — x86-only** | Scripps distribution has no arm64 build. If required, run under Rosetta 2 (`CONDA_SUBDIR=osx-64`) or prefer Meeko. Flagged, not blocking. |
| **numpy/pandas/biopython/biotite/pyyaml/tqdm/pytest** | conda-forge osx-arm64 | standard, expected clean. |
| **ESMFold / fair-esm / Chai-1 / GNINA / DiffDock** | **n/a locally — Vast-only** | Intentionally absent from the local env. Live in `vast/Dockerfile.fold`; pinned there on first build. |

## What WAS validated here (host-agnostic, pure-Python)

- `s3_fold.py --dry-run` — validates a FASTA and emits a valid job manifest
  (`run_location: vast`, per-seq sha256, fold params from config). ✅ runs & passes.
- Smoke suite **collects and runs**; the S3 dry-run test passes; every
  tool-dependent test SKIPS cleanly (tools absent on this container) rather than
  erroring. ✅
- `controls/fetch_controls.py` reachability: RCSB download of 6EQE returned the
  exact byte count recorded in `controls/MANIFEST.json` (614952 bytes). ✅

## S0/S1 stage decisions (CP0 audit for the dereplicate + tokenize PR)

These were resolved while implementing S0 (dereplicate) and S1 (tokenize). They
are recorded here because the stage code comments point back to this file.

### S1 backend — Foldseek-native ProstT5 (CHOSEN), transformers as fallback
`foldseek createdb --help` exposes **`--prostt5-model STR`** (and
`foldseek databases ProstT5` to fetch the weights), i.e. this Foldseek build does
native AA→3Di generation. We use it as the **primary S1 backend**:

- It runs the ProstT5 translation on **CPU**, so on the M4 it sidesteps both the
  MPS question and the arm64 `transformers`/`sentencepiece` fragility noted above.
- Its output is already a **valid Foldseek query DB** (the `_ss` 3Di sub-DB),
  which S2 consumes directly — no separate DB-build step.
- We also export an inspectable 3Di FASTA from the DB (`lndb` the headers onto
  `_ss`, then `convert2fasta`).

The `transformers` ProstT5 path (`s1_tokenize.tokenize_transformers`) is kept as a
**fallback for older Foldseek builds without the flag**; it honours
`s1_tokenize.{batch_size,device}` with a cpu fallback. Weights resolve from
`paths.prostt5_weights` → `PROTEUS_PROSTT5_MODEL` → auto-download via
`foldseek databases ProstT5` cached under `paths.models` (gitignored; ~2.4 GB
`prostt5-f16.gguf`). The 3Di model itself is logged on each run.

### MMseqs2 has no clustering RNG seed
`mmseqs easy-cluster --help` exposes **no `--seed`** — clustering is a
deterministic greedy set-cover, and the only stochastic knob is `--shuffle`
(input-order shuffle of the DB). So S0 does **not** fabricate a nonexistent seed
flag: it pins determinism with **`--shuffle 0`** and still reads + logs the single
global `random_seed` from config for provenance. (If a future seedable MMseqs2
step is added, it reads the same `random_seed`.)

### S2 fold-class reference — USE BOTH (decided; informs the next PR)
The single open decision flagged for the S2 PR — *what reference defines the
α/β-hydrolase "fold class" S2 searches against* — is resolved: **use BOTH
references, not one.**

- **Curated ESTHER / representative α/β-hydrolase set** — high-precision
  fold-class anchors (known α/β-hydrolase superfamily members).
- **Broad AF-DB / PDB search + post-filter** — keeps the divergent dark tail
  that a curated set alone would miss; the post-filter removes off-fold hits.

S2 searches the S1 query DB against both, then **unions** the survivors (curated
for precision, broad+post-filter for recall). This stays true to the unseeded
fold-CLASS intent: we match architecture, never specific PETase templates.

**Implemented + empirically validated** (`s2_foldclass_triage.py`). On the mini
corpus the Foldseek search of ProstT5-predicted 3Di separates the fold class
decisively at the configured thresholds (evalue≤0.01, bits≥50):

| query rep | best e-value | bits | verdict |
|---|---|---|---|
| IsPETase, LCC_WT, CalB, AChE, CRL, Est2 (α/β-hydrolases) | ≤ 4.5e-7 | ≥ 270 | **shortlisted** |
| decoy_allalpha, decoy_random (non-fold) | ≥ 6.6 | ≤ 4 | **dropped** |

i.e. a ~6-order-of-magnitude e-value gap. The 4 NEGATIVE controls
(CalB/AChE/CRL/Est2) are NOT PETases yet ARE shortlisted — proof S2 gates on
FOLD, not PETase homology; they are separated later at S4/S5. The reference DBs
used in production (curated ESTHER/representative set; broad AF-DB/PDB subset) are
host-local and pre-built — config `db` paths left blank to avoid guessing.

## Net env status (on this Linux container): YELLOW — scaffold only

Not GREEN, and cannot be from here: no Apple Silicon to exercise MPS, no osx-arm64
conda solve. GREEN is reachable **only on the M4** after `mamba env create`, a
clean import/run of each LOCAL tool, and a smoke suite where every LOCAL tool shows
positive output. There is **no hard GPU blocker** — folding (S3) is offloaded to
Vast.ai by design (`vast/Dockerfile.fold`), so MPS/CPU is sufficient locally.
