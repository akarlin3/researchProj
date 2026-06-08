"""S5 — Cleft filter: metrics A-E, control-anchored scoring.

Run fpocket on triad-positive models and score the active-site cleft on five
metrics (volume, hydrophobicity, depth, polarity, druggability). Scores are
anchored RELATIVE to the positive controls (e.g. IsPETase) from
controls/references.csv rather than to absolute cutoffs, so the bar tracks
what real PET hydrolases look like.

Reads fpocket_min_pockets / metrics / control_anchor from config (s5_cleft_filter).
"""

# TODO(P4): run fpocket per model, extract A-E metrics for the catalytic cleft,
#           z-score / rank against the control anchor, emit a ranked candidate table.
