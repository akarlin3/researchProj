# Echo — gates and locked thresholds

Pre-coding verification, in the spirit of `Minos/VERIFICATION.md`. Thresholds are **locked
before any run** (set at CP0) and are **not tuned** afterwards. Each gate is a hard halt;
reaching a named fallback (Lethe / Reverb) is a valid verdict, not a failure.

## CP0 — legitimacy gate (PASSED)
Echo's statistic must (a) not equate repeatability with accuracy/coverage, and (b) be
distinct from and beyond Gauge §4.2.2's rank check. Both proven analytically and re-checked
empirically by the CP1 self-test. → **PASS** (Echo viable; not routed to Lethe at CP0).

## CP1 — scaffold + method self-test (this checkpoint)
**Gate:** the subrepo resolves the same way Minos's does (`git log -- Echo/` shows Echo's
own commits); README registers Echo; package imports; the read-only Caliper import resolves;
nothing tainted. **Method self-test** (`scripts/run_harness.py`, SOLID) must satisfy:
- a correctly measurement-scaled 90% interval recovers the analytic ≈0.755 test–retest
  coverage (|error| < 0.01 at n=20000);
- **bias invariance**: a large systematic bias leaves coverage unchanged (|Δ| < 0.01) —
  the precision-not-accuracy guarantee;
- **scale sensitivity**: coverage tracks the analytic law across width scales 0.5–2.0;
- **distinctness from Gauge**: a pure width rescale leaves Spearman fixed (|Δ| < 1e-9) while
  moving coverage (> 0.05).

## CP2 — data gate (hard halt)
**Gate:** suitable public repeatability data is present for Echo's regime. ACRIN-6698
TrT0/TrT1 (n≈76) demonstrably exists (Gauge used it). `scripts/fetch_invivo.py --check`
confirms availability; `--from-gauge` materialises ROI-mean test–retest pairs
(download-on-demand, provenance only). **HALT → Reverb** (the synthetic test–retest harness
in `harness.py`) only if the real data is genuinely unsuitable for the scale regime.

## CP3 — validation (hard, honest gate) → RENDERED: **LETHE**
**Outcome (seed 20260613, n=76):** coverage_D = **0.263 [BCa 95% 0.158, 0.355]** vs the 0.755
target; R(D) = 0.247 (interval ~4× too narrow); under-scaled across the whole SNR grid →
**LETHE** under the locked gate below. D* over-scaled (coverage 0.797, R=19.5) by the
identifiability wall. See `LETHE.md`. CP4 is therefore not executed (PASS-only).

Run `scripts/run_validation.py` on the real data; report coverage + z-dispersion with BCa
95% CIs. **Locked gate on the well-identified parameter D:**
- **PASS:** coverage_D BCa-95%-CI **excludes 0.65** (not under-scaled) **and excludes 1.00**
  (has discriminative content); point estimate ∈ **[0.72, 0.96]**.
- **→ Lethe** if: coverage_D < 0.65 (interval can't capture scan–rescan scatter);
  coverage_D ≥ 0.99 or CI includes 1.0 (saturated, no scale content — overclaim/conservative);
  or the coverage signal is statistically redundant with Gauge's Spearman (collapse).
- **D\*** (under-identified): pre-registered identifiability **consistency probe**, not
  pass/fail — expect degraded/different behavior mirroring Gauge's D/D\* split.
- **No tuning:** α = 0.10 fixed; seeds fixed; thresholds above frozen.

## CP4 — result + manuscript (PASS only)
Build `paper/echo.tex` (ebgaramond + microtype). Every number traces to a seeded printout
via `paper/consistency.py` → `numbers.tex`. Every Fashion/Minos/Gauge-dependent claim carries
the PROVISIONAL marker. Honest scoping: precision ≠ accuracy; the explicit delta vs Gauge
§4.2.2; positioned in the no-ground-truth / QIBA repeatability literature.
