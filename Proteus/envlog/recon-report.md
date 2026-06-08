# Proteus — Checkpoint 0 Recon Report

**Generated:** 2026-06-08
**Host:** `vm` (Claude Code on the web — ephemeral cloud container)

## Probes

| Probe | Result | Verdict |
|---|---|---|
| OS / arch | Ubuntu 24.04.4 LTS, `x86_64`, kernel 6.18.5 | Linux native (not WSL — OK) |
| GPU | `nvidia-smi` absent; no `/dev/nvidia*`; no `nvidia` kmod; no NVIDIA/VGA on `lspci` | **NO NVIDIA GPU** |
| CUDA toolkit | `nvcc` absent; no `libcuda`/`libcudart` in `ldconfig` | **NO CUDA** |
| Python | 3.11.15 (system) | OK |
| conda / mamba | both absent | install needed (Miniforge) |
| pip | 24.0 | OK |
| Free disk | `/dev/vda`: 31 GB free of 252 GB | tight for corpora |
| Existing torch | not installed | clean |

## Verdict: BLOCKER (hard)

Named hard blocker hit: **"no NVIDIA GPU detected for the local folding env."**
Four independent signals confirm the absence of GPU/CUDA hardware — this is not a
driver/PATH glitch.

### Consequences for the local stack
- **ESMFold (S3)** cannot run meaningfully (GPU-only in practice).
- **PyTorch + CUDA** has no CUDA version to match (Checkpoint 2 premise void).
- **Smoke tests** for ESMFold / CUDA-torch can never reach positive-output GREEN here.
- Per the green-without-output rule, this env **cannot be marked GREEN** on this host.

### Out of scope (unaffected)
- Chai-1 / GNINA / DiffDock are Vast.ai burst targets (P3/P4) — not installed locally
  in any case; only their config is scaffolded.

## Decision taken
User selected **"Scaffold + CPU best-effort"** after the blocker was reported:
- Full repo scaffold is created.
- CPU-only / non-GPU tools are attempted (MMseqs2, Foldseek, fpocket, Vina, ProstT5 CPU,
  standard scientific Python).
- **ESMFold and CUDA-enabled PyTorch are intentionally skipped** on this host.
- Env status is **YELLOW / PARTIAL**, never GREEN, by design.

## Remedy to reach a real GREEN env
1. Run on the intended Linux/WSL2 host with the RTX 4080 + NVIDIA driver (>= R550 for
   CUDA 12.x) visible to `nvidia-smi`, then re-run all checkpoints.
2. Or re-provision this session with a GPU-backed compute policy:
   https://code.claude.com/docs/en/claude-code-on-the-web
3. Install conda/mamba (Miniforge) on whichever host — `environment.yml` and the bioconda
   tools assume a conda solver.
