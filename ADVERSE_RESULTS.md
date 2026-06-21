# ADVERSE_RESULTS — close/open-friction

Every experiment whose outcome did not simply confirm a claim is recorded here
verbatim, per the honesty-first contract. Confirmatory outcomes are in CHANGES.md.

Legend: **ADVERSE** = weakens a claim; **TEMPERED** = claim survives but must be
narrowed; **CONFIRMATORY (concession required)** = survives but a concession must
be stated.

---

## Railing is a SPECIFIC but NOT SENSITIVE flag of D\* unrecoverability — **TEMPERED**

Known-truth recovery (`phantom_recovery.json`, B3, over the trained prior): railed
voxels carry a 3.4× larger median absolute D\* recovery error than non-railed
(0.0307 vs 0.0090 at SNR 20), so railing **does** flag a genuinely harder subset.
But as a *classifier* of "large recovery error" the railed flag has **AUC ≈ 0.56
(SNR 20) / 0.52 (SNR 40)** — only marginally above chance. Reason: D\* is weakly
identified *broadly*, not only where the fit rails, and a voxel can have a bad
interior D\* without railing (railing under-flags). 

**Implication:** the paper must present railing as a *visible, specific* failure
flag, **not** a complete detector of D\* miscalibration. This is consistent with the
prior friction finding that a self-consistency gate's AUC against true parameter
error collapses toward chance. Stated as such in `sextant.tex` Sec. Robustness.

---

## In-silico railing MAGNITUDE under the full prior is ~30%, not 54.7% — **TEMPERED (regime, not misfit)**

`misspecification_isolation.json`: well-specified railing averaged over the full
trained prior is 29.8% (SNR 10), well below the real 54.7%. The gap is **not**
forward-model misfit (misspecification moves railing by ≤6.2 pp) — it is **regime
concentration**: restricted to the low-perfusion / high-D\* corner the curated real
ROI occupies, the well-specified rate rises to **50.7%** (SNR 10), bracketing 54.7%
(corner n=75, so noisy). The honest reading: the *mechanism* (weak-identifiability
railing) reproduces exactly in pure simulation with known truth; the *headline
magnitude* is a property of the curated ROI's regime, not a universal constant.
The paper does not claim the 54.7% magnitude is reproduced in unconditioned
simulation; it claims the mechanism is, and the corner brackets it.

---

## Model criticism finds REAL forward/noise misspecification in vivo (9.1%) — **CONFIRMATORY (concession required)**

`model_criticism.json`: the truth-free held-out-b criticism flags 9.1% of real
abdominal voxels as model-misspecified, a real excess over the 4.7% well-specified
baseline (and in the same direction as the prior NPE-based 12.8% on brain). So the
bi-exponential model **is** mildly misspecified on real tissue — this is conceded.
It does **not** undermine the railing claim, because the misspecification is
*anti-correlated* with railing (criticised|railed 3.3% < criticised|non-railed
16.2%): the railed voxels are the well-fit-but-unidentified ones. Concession
required: the in-vivo bi-exp fit carries ~9% genuine misspecification, reported as
a property of the data, separate from the railing.

---

## Model-criticism sensitivity is misspecification-type-dependent — **TEMPERED**

The self-calibrating (averaging-robust) criticism statistic is calibrated (4.7% on a
well-specified null) and detects **structural shape** misspecification (tri-exp
tissue tail: strong 6.2%, mild 5.8%) but is **robust to / insensitive to**
amplitude/noise deviations the bounded bi-exp fit can absorb (diffusional kurtosis
4.5%, baseline offset 2.1%). This is the correct behaviour for a *shape* check but
means the 9.1% real-data flag is a lower bound on total misspecification (amplitude
deviations are not counted). Stated honestly rather than presented as a universal
misspecification detector. (An external-SNR-normalised variant is sensitive to
amplitude in silico but is confounded by multi-acquisition averaging on the real
scan — hence the self-calibrating statistic is preferred and its scope disclosed.)

---

## HC4/CS3 budget sweep — **CONFIRMATORY (not adverse)**

`research_debt/budget_sweep/`: below-floor D\* overconfidence is flat across a 20×
budget range (claimed-ratio span 0.007; below-floor span 0.004) → budget-invariant,
not starvation. The honesty gate would have required qualifying the "intrinsic"
reading had the effect weakened at high budget; it did not. BANKED, not in the NMR
paper.
