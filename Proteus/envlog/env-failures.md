# Proteus — env resolution failures / deviations (CPU best-effort run)

Recorded per Checkpoint 2: "If any package will not resolve on this hardware,
record the failure in `envlog/` and continue." Host = CPU-only ephemeral cloud
container (Checkpoint 0 BLOCKER, see `recon-report.md`).

| Component | Status | Reason / deviation |
|---|---|---|
| **conda / mamba** | not used | Neither installed on host; `environment.yml` not solved. Used a pip venv + static binaries instead. |
| **PyTorch + CUDA** | **CPU build only** | No GPU/CUDA on host. Installed `torch 2.12.0+cpu`. CUDA-matched build cannot be resolved (no CUDA version to match). |
| **ESMFold (fair-esm)** | import OK, **inference SKIPPED** | `esm` imports fine, but `esmfold_v1()` inference needs a CUDA GPU. Smoke test skips by design. |
| **Foldseek** | OK | Installed `linux-avx2` static binary from mmseqs.com (bioconda equivalent). |
| **MMseqs2** | OK | Installed `linux-avx2` static binary from mmseqs.com. |
| **fpocket** | OK | Not in apt; built 4.2.2 from source (gcc 13.3.0). Default serial `make` (the `-j` parallel build races). |
| **AutoDock Vina** | OK | `vina` 1.2.7 python bindings via pip; scorer instantiates. |
| **Meeko** | OK | Needed extra deps `rdkit` + `gemmi` beyond a bare `pip install meeko`. |
| **ADFR suite** | not installed | Scripps distribution, not pip/conda installable. Deferred to real host (receptor prep). |
| **ProstT5 tokenizer** | OK (with caveat) | `transformers 5.10.2` mis-routes ProstT5's SentencePiece model through a tiktoken converter (raises on `spiece.model`). Worked around by loading the same asset via `sentencepiece` directly. On the real host, pin a transformers version whose T5 converter handles ProstT5, or load via sentencepiece. |
| **Chai-1 / GNINA / DiffDock** | n/a (by design) | Vast.ai burst targets (P3/P4); only `vast/` config scaffolded. |

## Net env status: YELLOW / PARTIAL

5/6 local smoke tests show positive output (Foldseek, MMseqs2, ProstT5, fpocket,
Vina). ESMFold is SKIPPED for lack of a GPU. Per the green-without-output rule,
the env is **not GREEN** on this host and cannot be until ESMFold folds on a real
CUDA GPU. Re-run on the RTX 4080 host to reach GREEN.
