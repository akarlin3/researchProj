# Proteus

Structure-first mining of the **dark proteome** for **divergent PET-hydrolase
candidates**. Rather than starting from sequence homology to known PETases — which
by construction can only re-find close relatives — Proteus folds and structurally
triages uncharacterized proteins, then gates on the *catalytic architecture* of a
serine hydrolase (α/β-hydrolase fold, Ser-His-Asp triad, oxyanion hole, a suitable
substrate cleft). Sequence clustering is used only to **dereplicate** the input
corpus, never to homology-gate it: dereplicate ≠ homology-gate, so the divergent
tail we actually care about stays in play.

## Pipeline stages (S0–S5)

| Stage | File | What it does |
|---|---|---|
| **S0** | `src/proteus/s0_dereplicate.py` | MMseqs2 clustering to **dereplicate** the corpus (collapse near-identical seqs). NOT a homology gate. |
| **S1** | `src/proteus/s1_tokenize.py` | ProstT5 translates sequence → Foldseek 3Di alphabet (cheap structural proxy). |
| **S2** | `src/proteus/s2_foldclass_triage.py` | Foldseek triage vs the **α/β-hydrolase fold CLASS** (unseeded — fold, not template). |
| **S3** | `src/proteus/s3_fold.py` | ESMFold batch on the local RTX 4080; mean-pLDDT filter; length-chunk long seqs. |
| **S4** | `src/proteus/s4_geometry.py` | Catalytic geometry gate: Ser-His-Asp triad + oxyanion-hole geometry. |
| **S5** | `src/proteus/s5_cleft_filter.py` | Cleft metrics A–E (fpocket), scored **anchored to the positive controls**. |
| P4 | `src/proteus/docking/` | AutoDock Vina (CPU) docking of PET-mimic ligands into ranked clefts. |

## Compute split: local vs Vast.ai

- **Local (CPU + single RTX 4080, 16 GB):** S0–S5 — ESMFold folding, Foldseek,
  MMseqs2, fpocket, and CPU Vina docking.
- **Vast.ai burst (P3/P4):** Chai-1 cofolding and GPU docking (GNINA / DiffDock)
  run as interruptible/on-demand burst jobs — see [`vast/`](vast/). These are
  **not** installed in the local env; only their config is scaffolded.

## Controls

`controls/references.csv` is the locked control set (positives: IsPETase, LCC_ICCG;
reference: LCC_WT; recovery: GuaPA, MG8; negative: CalB). Fetch the structures:

```bash
python controls/fetch_controls.py          # downloads PDBs -> structures/, writes MANIFEST.json
```

Structure rows are pulled from RCSB with a sha256 manifest. Sequence rows
(GuaPA, MG8) print manual-fetch instructions — see "Manual steps" below.

## Reproduce the environment

The local env targets a **CUDA-enabled RTX 4080 on Linux/WSL2**. On that host:

```bash
mamba env create -f environment.yml        # set pytorch-cuda to your recon CUDA version first
conda activate proteus
pytest tests/test_smoke.py -v               # every LOCAL tool must show positive output → GREEN
```

Pinned/resolved versions live in the lockfiles — treat them as generated, not
hand-edited:

- `requirements-lock.txt` — `pip freeze` of the resolved env.
- `envlog/conda-env-resolved.yml` — `conda env export` of the resolved env.
- `envlog/recon-report.md` — host recon (OS/GPU/CUDA) + GO/BLOCKER verdict.

> ⚠️ **This scaffold was generated on a CPU-only host with no GPU/CUDA** (see
> `envlog/recon-report.md`). The committed lockfiles reflect a **CPU best-effort**
> resolution: ESMFold imports but its folding is unvalidated, and `torch` is the
> CPU build. The env is **YELLOW/PARTIAL**, not GREEN. Re-run the install + smoke
> suite on the RTX 4080 host to produce the real CUDA-matched lock and a GREEN env.
> Per-component status is in `envlog/env-failures.md`.

## Reproducibility / random seed

The global `random_seed` lives in [`config/proteus.yaml`](config/proteus.yaml)
(`random_seed: 1729`) alongside all paths, corpus selection, and per-stage
thresholds. **Every stochastic step** (S0 MMseqs2, S3 ESMFold/torch, docking
search) must read the seed from there — do not hardcode seeds in stage code.

## Layout

```
config/proteus.yaml      paths, corpus, thresholds, RANDOM SEED
controls/                locked control set + fetch script + MANIFEST
src/proteus/             pipeline stages S0–S5, docking/, utils/
vast/                    Vast.ai burst image (Chai-1 + GNINA) + launch notes
tests/test_smoke.py      one positive-output test per local tool
envlog/                  recon report + resolved version snapshots
data/{raw,interim,processed}/   gitignored corpora/artifacts
```

## Manual steps left to the operator

1. **GuaPA sequence** — archaeal PETase (Acosta et al. 2025); fetch the protein
   sequence from the Marcotte Lab GitHub repo (`marcottelab/GuaPA`), save FASTA to
   `data/raw/`, record the commit hash.
2. **MG8 sequence** — saliva-metagenome PETase (Eiamthong et al., *Angew. Chem.*
   2022); locate the accession/sequence in the paper's SI, deposit FASTA in
   `data/raw/`.
3. **CalB / 1TCA** — confirmed: PDB **1TCA** is *Candida antarctica* lipase B
   (Uppenberg et al. 1994), the lid-bearing lipase used as the structural negative.
4. **ADFR suite** — install separately (Scripps) on the RTX 4080 host for receptor
   prep; not pip/conda installable.
5. Re-run install + `pytest tests/test_smoke.py` on the GPU host to reach a GREEN env.
