# Echo — repeatability as a *scale* check on conformal intervals

**Status: speculative gated build — VERDICT RENDERED. CP0 legitimacy PASSED; CP1 scaffold
up; CP2 data gate PASSED (real ACRIN-6698, n=76); CP3 → LETHE (constrained validation).**
On real data the conformal interval is ~4× too narrow to cover real test–retest variation of
D (coverage 0.263 [0.158, 0.355] vs 0.755 target) — a valid honest-limitation verdict, not a
failure. See [`LETHE.md`](LETHE.md) and [`results/RESULTS_VALIDATION.md`](results/RESULTS_VALIDATION.md).

Echo asks one ground-truth-free question about a deployed conformal IVIM interval:

> Is the interval the right **size** to capture a measurement's own irreproducibility?

It answers with **test–retest interval coverage** — does one scan's parameter estimate
fall inside the *other* scan's deployed conformal interval — reported per parameter with a
**BCa bootstrap CI**, on public same-day scan–rescan data (ACRIN-6698, n≈76).

## What makes Echo legitimate (the two hard constraints)

1. **Precision, not accuracy — provably.** Writing each estimate as
   `est = θ_true + bias + ε`, the test–retest discrepancy `Δ = est_B − est_A = ε_B − ε_A`
   **cancels bias exactly**. Echo's statistic is therefore invariant to any systematic
   error common to both scans (demonstrated in `echo_repeat/harness.py` and `tests/`). Echo
   certifies an interval is correctly **sized to measurement noise**; it is **blind to
   accuracy/bias** and makes **no ground-truth coverage claim**. A perfectly
   measurement-scaled 90% interval is *expected* to show ≈76% test–retest coverage, **not
   90%** — the derivable gap between accuracy-coverage and repeat-coverage is exactly why
   the two cannot be conflated.

2. **Distinct from and beyond Gauge — provably.** Gauge's published check (paper §4.2.2;
   the "§3.7" of the protocol) is a Spearman **rank** correlation: does the interval *width*
   widen where scan–rescan scatter is larger (D `r=+0.60`, D\* null)? Echo measures **scale**
   (a coverage rate). These are mathematically independent: rescaling every width by a
   constant leaves Spearman invariant but moves Echo's coverage arbitrarily. Gauge asks
   *"does the band widen where noise is larger?"*; Echo asks *"is the band the right size to
   capture that noise?"* — the question Gauge explicitly declined.

If on real data the coverage signal collapsed to Gauge's rank check, saturated (no scale
content), or under-scaled, Echo routes to **Lethe** (honest-limitation regime) — a valid
verdict. See `VERIFICATION.md` for the locked PASS/FAIL thresholds.

## Layout

```
echo_repeat/
  statistic.py    the core: test-retest coverage, standardized-residual scale check,
                  analytic reference, Spearman (for contrast), numpy-only BCa bootstrap
  harness.py      synthetic test-retest generator + the CP1 method self-test
                  (also the Reverb fallback spec)
  invivo.py       IVIM forward + segmented plug-in fit + Caliper-conformal deployer +
                  real test-retest signal loader
  provenance.py   download-on-demand provenance manifest writer (mirrors Gauge's posture)
  _paths.py       read-only import chokepoint (Caliper SOLID; Gauge/Fashion/Minos PROVISIONAL)
scripts/
  run_harness.py     CP1 method self-test (SOLID) -> results/RESULTS_HARNESS.*
  fetch_invivo.py    CP2 download-on-demand fetch (reuses Gauge's data-handling template)
  run_validation.py  CP3 real-data validation, locked gate -> results/RESULTS_VALIDATION.*
paper/            CP4 manuscript (ebgaramond + microtype), built PASS-only
tests/            unit tests for the statistic, harness, and fit
ASSUMPTIONS.md    pinned Fashion/Minos/Gauge inputs + SOLID/PROVISIONAL split
PROMOTION.md      promotion (PASS) / Lethe-fold / Reverb-fold paths
VERIFICATION.md   the CP gates and locked thresholds
reproduce.sh      one-command re-validation (CP1 -> CP2 -> CP3 -> CP4)
```

## Reproduce

```bash
pip install -e .            # numpy-only core
bash reproduce.sh          # CP1 self-test always; CP2/CP3 run iff data present; CP4 PASS-only
```

## Data & IP

Echo's own repo and history are **synthetic + open only**. Real repeatability data
(ACRIN-6698, CC-BY-4.0, DOI 10.7937/tcia.kk02-6d95) is **download-on-demand**: no pixel
data is committed, only a provenance manifest. Echo redistributes nothing and mirrors
Gauge's in-vivo data posture exactly.

## License

MIT (see `LICENSE`). Echo imports Caliper (MIT) read-only and depends on Gauge/Fashion/Minos
by read-only import; those dependencies are **PROVISIONAL** (in review) — see `ASSUMPTIONS.md`.
