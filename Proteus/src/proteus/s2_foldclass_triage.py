"""S2 — Fold-CLASS triage with Foldseek (unseeded).

Match candidates against the alpha/beta-hydrolase fold CLASS, not against
specific PETase templates. The point is to retain anything with the right
*architecture* regardless of sequence-level homology to known enzymes —
keeping the divergent tail in play.

Reads target_fold / evalue / min_bits from config/proteus.yaml (s2_foldclass_triage).
"""

# TODO(P2): build a Foldseek 3Di db from S1 output, search vs an alpha/beta-hydrolase
#           fold-class profile/db, keep hits above evalue/min_bits; write survivors list.
