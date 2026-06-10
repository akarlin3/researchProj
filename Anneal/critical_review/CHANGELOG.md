# Project Anneal — PRE Critical-issue review (C1, C2): findings changelog

**Mode:** run-first, report-only. No manuscript (`.tex`) edits were made. Every number below
traces to printed output from analysis scripts in `critical_review/` (re-runnable). All Weibull
fits use the audited profile-likelihood censored-MLE fitter `anneal-hazard/src/survival.py`
(`fit_all`).

**Independent verification:** every headline number was re-derived by an adversarial verification
pass that used *different* implementations — `lifelines.WeibullFitter` (where available), two
from-scratch scipy censored-Weibull MLEs (profiled + joint), and an MLE-free Kaplan–Meier
Weibull-plot slope — plus refutation/sensitivity probes. **All three checkpoints CONFIRMED**;
independent fits agree with the author fitter to <1e-6 in k; no probe broke any claim (several
strengthened them). Scripts in `critical_review/verification/`.

Two distinct models are involved (see `critical_review/` and CP0 diagnosis):
- **C1** → driven two-population mean-field Sakaguchi–Kuramoto "absorption" campaign at A=0.5
  (main checkout data: `absorption_results/`, `reduced_results/`).
- **C2** → nonlocally-coupled **ring** (Wolfrum–Omel'chenko), this worktree
  (`anneal-hazard/`), death = ρ_std(t)<0.04 held 50, β∈{0.110…0.130}, N∈{32,64,128,192,256}.

---

## C1 — Aging on a deterministic system (A=0.5): **SURVIVED**

**Worry:** at A=0.5 the flow is deterministic, so pooled Weibull k_abs>1 might be the spread of
fixed transient times across the seed ensemble (IC heterogeneity), not a rising state-level hazard.

**Test:** stratify the 200-seed ensemble per N into 4 quantile bins by the collective IC and
refit a right-censored Weibull within each stratum (`critical_review/cp1_stratified.py`).

**Result — within-stratum k stays > 1 with CI excluding 1, unanimously:**

| binning variable | cells | k range | CI excludes 1 |
|---|---|---|---|
| pooled per-N (provenance) | 7 | 2.22 – 3.03 | 7/7 — **reproduces the manuscript k_abs exactly** |
| \|Δφ₀\| (collective phase gap) | 28 | 1.88 – 5.12 | **28/28** |
| reduced-model predicted lifetime t_capture (decisive) | 28 | 2.00 – 5.00 | **28/28** |
| R̄_incoh,0 | 28 | 1.88 – 4.83 | **28/28** |

- Pooled per-N k_abs = **2.27, 2.22, 2.49, 2.67, 2.44, 2.72, 3.03** (N=4→64) — matches the paper's
  §`sec:stats` values bit-for-bit (fitter provenance).
- Conditioning on **t_capture** (the reduced ODE's full IC→lifetime map) is the strongest test:
  even within a narrow predicted-lifetime band the measured lifetimes keep a Weibull k≥2.0 shape
  (CI excludes 1) in every cell; LRT vs exponential p ≤ 1e-7 throughout.
- Within-bin CV(τ_abs) is substantial (**0.21–0.57**, not a degenerate spike) and tracks the
  fitted Weibull k quantitatively (k≈2→CV≈0.52, k≈5→CV≈0.22), confirming a genuine aging-shaped
  within-bin lifetime spread rather than heterogeneity.
- **Ratchet cross-check (consistent):** per-pass hazard increment rises 0.031→0.052→0.078→0.091
  with N; ratchet_frac→1.0 by N≥32 — the same IC-independent "hazard rises as the chimera ages."

**Verdict (spec decision rule):** within-stratum k>1 with CIs excluding 1 across N → **aging is
trajectory-level**, surviving conditioning on the collective initial condition. C1 holds.

*Independent verification (CONFIRMED):* reproduced with lifelines + a from-scratch profile-likelihood
MLE + an MLE-free KM Weibull-plot slope (k∈2.07–2.96), triangulating the shape from three statistical
principles. Refutation probes (4/6/8/10 bins; binning by Rsync0/Rincoh0/dphi0; joint 2-D & 3-D IC
bins; pooled 20-bin) never drive k→1 except in isolated <22-run cells (small-sample noise); within
narrow strata k *grows* with N — the opposite of an IC-confound signature.
*Honest refinement:* at A=0.5 the censored fraction is **exactly 0%** (all 1400 runs absorb before
t_max=2000), so these are effectively *uncensored* MLEs — the censoring machinery is not exercised
here, which makes the k>1 result cleaner rather than weaker.
Artifacts: `cp1_stratified/cp1_stratified.json`, `cp1_strata_table.csv`, `cp1_k_vs_bin.{png,pdf}`.

---

## C2 — Ring lifetime inversion: **SURVIVED (and strengthened)**

### CP2 — is the dying object a canonical chimera? **VALIDATED**

**Worry:** the ρ_std<0.04 criterion might time a non-canonical decay channel, contaminating k(N)/k∞.

**Test:** the spatial field ρ_k(t,·) was never stored, so all 300 seeds per cell (β=0.110 & 0.130
at N=256) were **re-run with field dumping** (deterministic; τ reproduces the ensemble **bit-for-bit,
0/600 mismatches**), and the pre-death state classified (contrast gate ρ_std>0.08; head count =
dominant spatial Fourier mode of ρ_k). `critical_review/cp2_batch.py`, `cp2_fields.py`.

| cell | canonical single-arc | multi-headed chimera | degenerate channel | k unfiltered | k chimera-death filt. | k single-arc filt. |
|---|---|---|---|---|---|---|
| β=0.110 N=256 | 83.7% | 14.3% | **2.0%** | 1.220 [1.12,1.32] | 1.241 [1.14,1.35] | 1.241 [1.13,1.36] |
| β=0.130 N=256 | 89.3% | 8.0% | **2.7%** | 1.466 [1.35,1.58] | 1.481 [1.36,1.60] | 1.462 [1.34,1.59] |

- Unfiltered k reproduces the published VERDICT values (1.220, 1.466) **exactly**.
- **97–98% are genuine chimera deaths** (single + multi-head). The degenerate "dissolve-into-limbo"
  channel C2 feared is **real but rare (2–3%)**.
- Removing it **does not weaken k — it slightly strengthens it** (CIs still exclude 1); strict
  single-arc-only refit is also k>1. → object **validated**, k∞ claim not contaminated.

*Independent verification (CONFIRMED):* 0/600 τ mismatches re-derived with an independent detector;
k reproduced by two from-scratch MLEs (agree to <7e-7). Gate sweep c∈{0.06,0.08,0.10,0.12}: degenerate
fraction rises (β=0.130: 0.7/2.7/9.0/19.0%) but the filtered k moves *away* from 1 (1.470→1.569, all
CI>1); dropping the weakest-structured decile *raises* k; even the worst structured-life decile alone
fits k>1. The aging signal is broad-based, not an artifact of the rare degenerate channel.
Artifacts: `cp2_validation/cp2_summary.json`, `cp2_runs_b*.csv`, `cp2_representative.{png,pdf}`.

### CP3 — is the trend / hazard shape an artifact of the criterion? **ROBUST**

**Test:** re-detect death on the *same* β=0.130 trajectories (all 5 N) under alternative criteria,
changing only the label (`critical_review/cp3_criterion.py`). N∈{32,64,128} from the original saved
traces; N=192 regenerated; N=256 from the CP2 re-run (all reproduce τ bit-for-bit).

**Median lifetime DECREASES with N under every criterion** (ratio τ̃(256)/τ̃(32)):

| criterion | N=32 | 64 | 128 | 192 | 256 | ratio |
|---|---|---|---|---|---|---|
| original ρ_std<0.04 | 872 | 422 | 404 | 349 | 345 | 0.40 |
| mean-coh ρ_mean>0.78 | 810 | 379 | 362 | 310 | 276 | 0.34 |
| struct-loss ρ_std<0.08 | 828 | 383 | 372 | 314 | 308 | 0.37 |
| struct-loss ρ_std<0.10 | 815 | 372 | 356 | 303 | 289 | 0.35 |

→ the lifetime **inversion is criterion-robust** (real physics of the near-boundary regime, not a
detector artifact). Provenance: original median-τ = 872/422/404/349/345 matches the published
SUMMARY exactly.

**Weibull k(N) is k>1 under every ρ_std-based criterion:**

| criterion | N=32 | 64 | 128 | 192 | 256 | all CI>1 |
|---|---|---|---|---|---|---|
| original ρ_std<0.04 | 1.33 | 1.55 | 1.68 | 1.68 | 1.47 | ✅ (== VERDICT) |
| struct-loss ρ_std<0.08 | 1.28 | 1.45 | 1.61 | 1.56 | 1.42 | ✅ |
| struct-loss ρ_std<0.10 | 1.26 | 1.41 | 1.59 | 1.50 | 1.40 | ✅ |
| mean-coh ρ_mean>0.78 (raw) | 1.26 | 1.44 | **0.88** | 1.13 | **0.62** | ❌ (artifact) |
| mean-coh, artifact-corrected | 1.26 | 1.44 | 1.62 | 1.54 | 1.57 | ✅ |

- The mean-coherence-level criterion's dip <1 at large N is a **diagnosed censoring artifact**:
  6.7% of N=256 collapses go to a spatially-uniform **twisted** state (ρ_k≈0.50, ρ_std≈0 → genuinely
  dead) rather than the high-coherence sync plateau (~0.83). A ρ_mean>0.78 rule misses these →
  spuriously censors 20 runs at T_max=12000 → drags k to 0.62. **Excluding the mis-censored runs
  restores k=1.573 [1.44,1.71]**, consistent with all other criteria.
- This actually demonstrates *why* the ρ_std spatial-homogeneity criterion is the appropriate
  detector: it captures both sync- and twisted-collapse channels.

**Verdict:** the decreasing-median-with-N trend and the k>1 rise-then-saturate hazard structure are
**robust to the death criterion**; the one criterion that appears to break k>1 does so via a
diagnosed artifact, not real memorylessness. C2 holds.

*Independent verification (CONFIRMED):* all 20 mean-coherence-censored N=256 runs independently
confirmed twisted (final ρ_mean=0.502, max ρ_mean ever <0.79) **and** genuinely ρ_std-dead at
t=642–2298 — so they are real deaths the rule mis-censors, not long-lived survivors. Re-labeling them
as deaths at their true ρ_std<0.04 time gives k(N) all>1; censoring is exactly 0% under every ρ_std
criterion (so the median decrease cannot be a censoring artifact); and an independent KM
log-cumulative-hazard slope = [1.51,1.96,2.01,2.08,1.82] (all>1) corroborates increasing hazard
without the Weibull MLE.
Artifacts: `cp3_criterion/cp3_criterion.json`, `cp3_criterion.{png,pdf}`.

---

## Bottom line

- **C1 SURVIVED** — mean-field aging at A=0.5 is trajectory-level (within-stratum k>1, CI excl 1, 84/84
  cells; survives conditioning on the reduced-model lifetime prediction).
- **C2 SURVIVED & STRENGTHENED** — the ρ_std<0.04 criterion kills a genuine chimera in 97–98% of
  runs; the lifetime inversion and k>1/saturation are criterion-robust; a contamination channel exists
  but is rare and does not move k.

No headline claim failed its analysis. (Per the spec, had one failed, this file would report the
numbers and stop without any manuscript edit.)
