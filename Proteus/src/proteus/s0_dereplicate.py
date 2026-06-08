"""S0 — Dereplicate the input corpus with MMseqs2.

This is DEREPLICATION, not a homology gate: we collapse near-identical
sequences (>= s0_dereplicate.min_seq_id) to a representative set so downstream
folding isn't wasted on duplicates. We deliberately do NOT filter by similarity
to known PETases — that would discard the divergent dark-tail candidates we are
hunting for.

Reads thresholds and the random seed from config/proteus.yaml.
"""

# TODO(P1): wrap `mmseqs createdb / cluster / createtsv`, honor min_seq_id,
#           coverage, cov_mode from config; pass --seed from random_seed; emit a
#           representative FASTA + cluster TSV into data/interim.
