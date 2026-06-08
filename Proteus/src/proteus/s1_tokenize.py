"""S1 — Tokenize sequences to 3Di with ProstT5.

ProstT5 (Rostlab/ProstT5) translates amino-acid sequence into Foldseek's 3Di
structural alphabet without folding, giving a cheap structural representation
for the fold-class triage in S2.

Reads device/batch settings from config/proteus.yaml (s1_tokenize).
"""

# TODO(P1): load ProstT5 via transformers (weights path from config.paths.prostt5_weights),
#           run AA->3Di translation in batches, write 3Di strings to data/interim.
