# FLAGS.md — needs-compute / needs-data / needs-decision (friction remediation)

Separates what was *done* this run from what remains. Done items are in
CHANGES.md; claim-weakening outcomes are in ADVERSE_RESULTS.md.

## Done this run (no flag)
- CP0 runnability table (`paper/FRICTION_DIAGNOSIS.md`).
- CP2 bias-aware floor (`npe/run_cp2_bias_aware_floor.py`, `npe/cp2_bias_aware_floor.csv`) + manuscript reframe.
- CP3 NLLS railing audit (`npe/run_cp3_railing_audit.py`) + S4/main-text both-ways reporting.
- CP6 graveyard positioning + 6 verified references.
- CP1 budget sweep (`npe/run_cp1_budget_sweep.py`) — **executed** in the isolated venv.
- CP4/CP5 decorrelation + in-silico control (`npe/run_cp45_decorrelate_control.py`) — **executed** on the canonical model.

## needs-external-data (CC cannot produce)
- **HC2 / CS2 — real in-vivo parameter ground truth.** The in-silico control (CP5)
  separates posterior-width overconfidence from forward misfit *in simulation*, and
  an RNPE-style criticism is run on the in-vivo signal, but a true in-vivo
  parameter-recovery test needs a **phantom or repeated-acquisition dataset with
  known/repeatable parameters**. No such dataset is in the repo or cited. The in-vivo
  limb remains a held-out-signal / OOD check, not a parameter-recovery validation.
- **HC6 in-vivo 3-way b-split decorrelation.** The simulation decorrelation (CP4) is
  done with ground truth. Repeating it on the *in-vivo* brain data with a disjoint
  third b-subset needs the raw signals (`download/Data/brain.*`, Zenodo 14605039) +
  a trained `.pt`; the committed CSVs do not store the posterior-predictive samples
  needed to reconstruct it. Downloadable but not executed this run.

## needs-compute (runnable; not all run this round)
- 5M-budget point of the CP1 sweep ("if feasible" in the brief) — not run
  (~hours of additional CPU); the {50k…1M} sweep is sufficient for the verdict.
- Per-protocol recalibration experiment (the proposed safeguard) — not implemented;
  it inherits Ovadia's shift-degradation caveat (now stated, CP6) and is future work.

## needs-decision (human / author)
- Hero & Fessler (1993) exact volume/pages for typeset (`note` flag in `refs.bib`).
- Whether to regenerate Figures 3/4 and the supplement figures from freshly trained
  models for the camera-ready (current committed figures predate this env; the
  committed CSV numbers are unchanged and the new analyses are additive).
- Reconciling the in-vivo abdominal `figS4` voxel count (1618 here vs "1618 high-SNR
  ROI voxels"/"54.7%" in text) at typeset — consistent in this run.

## Reproduce
Isolated env: `python -m venv --system-site-packages $VENV` (from proteus python),
`pip install sbi==0.26.1`; run with `KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=4
PYTHONPATH=.`. Canonical NPE: `train_npe.py --mode set --budget 500000 --epochs 200
--log-dstar --seed 0`.
