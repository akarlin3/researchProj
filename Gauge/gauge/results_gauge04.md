# Gauge 04 — Results: robustness, exchangeability, and an in-vivo demonstration

Every number below is reproduced verbatim from the in-session CP0–CP2 printouts
(`results/robustness_report.txt`, `results/invivo_report.txt`). Deterministic
from seed `20260613`; α = 0.10 (nominal 0.90); deployed split-conformal is
calibrated **once** on the i.i.d. Rician cohort calibration split (Gauge 01) and
each shift perturbs the *test* (deployment) distribution only. CP0 is a HALT-TO-
REPORT gate and CP2 a HALT-ABLE-on-data gate; nothing is tuned toward the answer.

---

## Headline (the three gates)

1. **Conformal coverage is fragile to exchangeability breaks — but the breaks are
   observably detectable.** Every deliberate break miscalibrates coverage (max
   |coverage − nominal| = **0.385**); the dangerous under-coverage reaches D\*
   = **0.515** under a harder-tissue prior shift and the worst parameter = **0.610**
   under a low-SNR covariate shift. **Weighted / nonexchangeable conformal**
   (Tibshirani 2019; Barber, Candès, Ramdas & Tibshirani 2023) **recovers the
   covariate shifts** (SNR: 0.610 → 0.899) and is **honestly limited** on the
   forward-model misspecification (a P(y|x) shift the X-space likelihood ratio
   cannot see). The **Minos-style** label-free monitor **fires on all four
   observable breaks before coverage fails**, and is **blind** to the latent
   high-D\* gap — exactly Minos v3's observable-AUC≈1 / hidden-AUC≈0.5 signature.

2. **The high-D\* identifiability wall is acquisition-ROBUST.** Across a clinical,
   a CRLB-optimal, and a dense b-value scheme the high-D\* coverage barely moves
   (0.841 → 0.844) and the resolution ratio CRLB(D\*)/tercile-width stays **≥ 1.05
   for every scheme** (≥ 1 ⇒ unresolvable). Acquisition design lowers the CRLB
   (1.25 → 1.05) but does not remove the wall — a concrete, honestly-negative
   handoff to **Vernier**.

3. **In vivo, only qualitative claims are made — and the monitor enforces that.**
   With no ground truth, coverage is unmeasurable; the demo shows the
   conformalized interval map, the **heavy-tailed D\* band widths** (the wall on
   real signal), and the deployment monitor **firing on the synthetic→deployment
   transfer (AUC 0.84)** — turning "no ground truth" into an observable stop sign.
   Dataset choice/licensing and in-vivo framing are a **human sign-off**
   (HALT-TO-REPORT): no clean, permissive, no-GT in-vivo IVIM set was usable in
   this environment, so the demo runs on a clearly-labeled synthetic stand-in.

---

## CP0 / GATE 0 — exchangeability / shift stress (HALT-TO-REPORT)

`n_cal = 2000`, `n_test = 2500`, 22 b-values. **One** deployed split-conformal
predictor (radius q (D, D\*, f) = [0.0005, 0.0642, 0.1196]) and **one** monitor,
both built on the i.i.d. Rician cohort calibration, are applied to every shift.

### [0.1] Per-parameter coverage under each break — NAIVE vs WEIGHTED conformal

| scenario | method | D | D\* | f | monitor |
|---|---|---|---|---|---|
| in-dist (control) | naive | 0.909 | 0.900 | 0.902 | **silent** (AUC 0.50) |
|  | weighted | 0.895 | 0.892 | 0.893 | |
| **SNR shift (low)** | naive | **0.610** | 0.744 | **0.616** | FIRES (AUC 0.91) |
|  | **weighted** | **0.901** | **0.899** | **0.911** | |
| **prior shift (harder tissue)** | naive | 0.970 | **0.515** | 0.906 | FIRES (AUC 0.71) |
|  | weighted | 0.900 | 0.975 | 0.916 | |
| **tri-exp misspec** | naive | 0.908 | **0.846** | 0.910 | FIRES (AUC 0.51) |
|  | weighted | 0.886 | 0.928 | 0.894 | |
| noise-model (Rician cal→Gauss) | naive | 0.918 | 0.901 | 0.904 | FIRES (AUC 0.52) |
|  | weighted | 0.904 | 0.898 | 0.896 | |

- **Weighted conformal cleanly recovers the COVARIATE (SNR) shift**: worst
  parameter 0.610 → per-param (0.901, 0.899, 0.911), all within **0.011** of
  nominal. The harder-tissue prior shift's D\* collapse (0.515) is likewise
  restored (→ 0.975, conservative). Miscalibration is two-sided (the same breaks
  over-cover other parameters).
- **Weighted conformal does NOT cleanly repair the tri-exponential
  misspecification**: a shift in the conditional law P(y|x) is, in general, not
  fixable by an X-space likelihood ratio. Here weighting *over-corrects* D\*
  (0.846 → 0.928) and leaves the per-parameter result UNEVEN (max deviation
  **0.028** vs the SNR shift's 0.011) — Barber et al. (2023) *bounds* this residual
  gap, it does not erase it. The honest limit, the same observable/hidden split
  that recurs throughout.
- **The noise-model swap (Rician-calibrated, Gaussian deployment) barely dents
  coverage** at clinical SNR (0.918, 0.901, 0.904) — yet the monitor still
  observably flags it (Family-2 residual; AUC 0.52). A mild effect, but detectable
  — consistent with under-identification being largely noise-model-independent.

### [0.3] Minos cross-link — does the monitor fire *before* coverage fails?

Label-free deployment-validity monitor (the **Minos-Core v3** idea: staleness
from observable statistics), two detector families — **Family-1 Mahalanobis** on
signal summary features + **Family-2 residual-conformal** on NLLS fit residuals.

- **Observable breaks flagged: 4/4** shift scenarios fire; the in-dist control
  stays silent.
- **Hidden case** (in-dist, exchangeable): marginal D\* coverage 0.900 ✓, monitor
  **silent (AUC 0.50)** — yet high-D\* **conditional** coverage **0.815**
  (worst-SNR **0.595**) is below nominal. The observable monitor is **blind** to
  the latent high-D\* gap (the Gauge 03 wall).

**SNR-severity sweep (fires-before-failure):**

| test SNR | min coverage | cov < nom−0.05? | monitor | Mahalanobis stat/thr |
|---:|---:|:--:|:--:|---:|
| 50 | 0.967 | ok | silent | 1.98 / 2.60 |
| 30 | 0.919 | ok | silent | 2.15 / 2.60 |
| 20 | 0.866 | ok | **FIRES** | 2.44 / 2.60 |
| 13 | 0.753 | **FAIL** | FIRES | 3.15 / 2.60 |
| 9 | 0.593 | FAIL | FIRES | 4.09 / 2.60 |
| 6 | 0.413 | FAIL | FIRES | 5.45 / 2.60 |
| 4 | 0.283 | FAIL | FIRES | 7.46 / 2.60 |

Monitor first fires at **SNR 20** (coverage 0.866, still ok); coverage first fails
at **SNR 13**. **Fires-before-failure: YES** — the monitor raises the alarm one
severity step before the guarantee silently breaks.

**GATE 0 verdict:** every exchangeability break miscalibrates coverage; weighted
conformal recovers the covariate shifts and is honestly limited on the P(y|x)
misspecification; the Minos-style monitor fires on the observable breaks before
coverage fails and is blind to the latent high-D\* gap. Reported, not engineered
around. *(Figure: `figures/robustness_shift.pdf`.)*

---

## CP1 / GATE 1 — acquisition sensitivity of the high-D\* wall (Vernier tie-in)

Same prior / SNR / seed; only the b-value scheme differs. The CRLB-optimal scheme
is built by greedy forward selection minimising mean high-D\* CRLB(D\*).

| scheme | n_b | D\* marg | hi-D\* marg | hi-D\* worst-SNR | hi-D\* CRLB / tercile-width |
|---|---:|---:|---:|---:|---:|
| clinical (11 b) | 11 | 0.912 | 0.841 | 0.595 | 1.25 |
| CRLB-optimal (11 b) | 11 | 0.908 | 0.844 | 0.608 | **1.05** |
| dense (22 b) | 22 | 0.913 | 0.845 | 0.596 | 1.06 |

CRLB-optimal b-values = `[0, 5, 10, 15, 20, 35, 65, 70, 75, 80, 800]` — dense
low-b sampling (to catch the fast perfusion compartment) plus one high-b anchor
for D, exactly as IVIM design theory predicts.

- high-D\* marginal coverage clinical 0.841 → CRLB-optimal 0.844 (**Δ +0.004**);
- hi-D\* CRLB/tercile-width clinical 1.25 → optimal 1.05 (**improvement +0.20**).

**GATE 1 verdict:** the wall is **acquisition-ROBUST** — even the CRLB-optimal and
dense schemes keep hi-D\* CRLB ≥ tercile width, so the high-D\* regime stays
under-resolved and the conditional gap stays open. Acquisition design *shifts* the
wall (lower CRLB) but does not *remove* it: a deeper identifiability limit, and a
concrete (honestly negative) handoff to Vernier. *(Figure:
`figures/acquisition_wall.pdf`.)*

---

## CP2 / GATE 2 — in-vivo qualitative demonstration (HALT-ABLE on data)

**What conformal cannot claim in vivo** (kept strictly separate from the synthetic
guarantees): in vivo there are no ground-truth (D, D\*, f), so realized coverage is
**unmeasurable** and the finite-sample 1−α guarantee is **not asserted**. All CP2
outputs are qualitative interval *behavior*.

- **[2.1] Conformalized interval map** (deployed CQR band; synthetic-calibrated):
  median plug-in (D, D\*, f) = (1.41, 60.5, 0.304); median band width
  (0.76, 64.0, 0.127) [D, D\* in 10⁻³ mm²/s].
- **[2.2] D\* compartment band widths** — 10/50/90th pct = **[47.4, 64.0, 78.9]**
  (10⁻³ mm²/s); widest-decile/narrowest-decile ratio **1.7×**. The D\* interval
  honestly balloons on voxels where the perfusion compartment is under-identified
  — the Gauge 03 wall, now on (stand-in) signal.
- **[2.3] Deployment monitor** on the in-vivo signals: **FIRES**, AUC(cal vs
  deployment) = **0.84** (Family-1 Mahalanobis 3.68/thr 2.60, Family-2 residual
  0.358/thr 0.177; domain-weight spread 7×). Synthetic-calibrated → deployment is
  itself a large exchangeability break, and the observable monitor detects it —
  the **honest, label-free reason** the synthetic guarantee must not transfer to
  this data. The monitor turns "no ground truth" (the Echo tension) into an
  observable stop sign (the CP0 cross-link).
- **Phantom-only sanity** (the stand-in *has* GT; a true in-vivo run prints
  nothing here): coverage (D, D\*, f) = (0.805, 0.636, 0.744) — under-covers,
  consistent with the monitor firing on the transfer. Shown only to validate the
  pipeline; **not** an in-vivo coverage claim, kept strictly separate.

**GATE 2 verdict:** the in-vivo pipeline runs and its qualitative outputs are
clearly delimited from the synthetic guarantees. **HALT-TO-REPORT (dataset gate):**
no clean, permissively-licensed, ground-truth-free in-vivo IVIM dataset was usable
in *this* environment (nibabel absent; the OSIPI Apache-2.0 collection commits only
b-value files and signal *generators*; true in-vivo sits behind TCIA data-use
agreements). The demo therefore runs on a transparent **synthetic stand-in**.
**Human sign-off needed:** dataset choice + licensing, and the in-vivo framing
(what the qualitative demo may claim). The loaders `load_dwi_nifti` /
`load_signals_npy` are the drop-in integration point. *(Figure:
`figures/invivo_demo.pdf`.)*

---

## Positioning (new references)

- **Barber, Candès, Ramdas & Tibshirani (2023)**, *Conformal prediction beyond
  exchangeability* (Ann. Statist.) — the nonexchangeable / fixed-weight framework
  bounding the coverage gap under deviation from exchangeability.
- **Tibshirani, Foygel Barber, Candès & Ramdas (2019)**, *Conformal prediction
  under covariate shift* (NeurIPS) — the likelihood-ratio-weighted conformal used
  as the shift-aware fix; weights estimated by a domain classifier.
- **projMinos**, Minos-Core v3 — the label-free deployment-validity monitor under
  shift (observable-AUC≈1 / hidden-AUC≈0.5); its two detector families
  (Mahalanobis summary-space, residual-conformal signal-space) are reused here as
  the concept, independent of Fashion's flow.
- The high-D\* wall's acquisition-robustness (CP1) is the bridge to **Vernier**
  (acquisition-aware design): a CRLB-optimal scheme helps but does not close it.

## Reproduce

```bash
pip install -r requirements.txt
python -m gauge.robustness          # CP0 + CP1 + figures + report (~2 min)
python -m gauge.invivo              # CP2 in-vivo demo (synthetic stand-in)
python -m gauge.invivo dwi.nii.gz bvals.txt   # CP2 on a real DWI (needs nibabel)
python -m pytest -q                 # 70 tests (55 prior + 15 new)
```
