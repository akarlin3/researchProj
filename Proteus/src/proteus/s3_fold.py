"""S3 — Fold survivors with ESMFold on the local RTX 4080.

Batch-folds the S2 survivors, applies a mean-pLDDT filter, and length-chunks
long sequences to fit VRAM. GPU-only in practice — on a host without CUDA this
stage cannot run (see envlog/recon-report.md).

Reads device / plddt_min / chunk_size / max_recycles from config (s3_fold) and
seeds torch from random_seed.
"""

# TODO(P3): load ESMFold (fair-esm), set torch manual_seed(random_seed), batch-fold
#           with length chunking, write PDBs to structures/ and per-model pLDDT arrays,
#           drop models below plddt_min.
