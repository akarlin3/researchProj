# Echo — assumptions manifest

Echo's validation method (the test–retest *scale* statistic) is **SOLID** — it depends on
no upstream paper and is machine-checked by `scripts/run_harness.py`. Echo's *real-data
result and its framing* depend on inputs that are **PROVISIONAL** (in review). This file
pins every such input and states what becomes invalid if it changes.

## 0. SOLID vs PROVISIONAL split

| Component | Depends on Fashion/Minos/Gauge? | Status |
|---|---|---|
| The statistic (`statistic.py`) + analytic reference | No — self-contained, machine-checked | **SOLID** |
| Method self-test (`run_harness.py`, CP1) | No — pure synthetic | **SOLID** |
| Conformal deployer (`invivo.py`, via Caliper) | Caliper only (MIT, in-tree, stable) | **SOLID** (ruler) |
| Real-data validation (`run_validation.py`, CP3) | Yes — Caliper ruler on real signals; Gauge fetch/baseline | **PROVISIONAL** |
| Manuscript (`paper/echo.tex`, CP4) | Yes | **PROVISIONAL** (caveated inline) |

## 1. CALIPER (SOLID — the ruler)

| key | pinned value | source | role for Echo |
|---|---|---|---|
| `caliper.version` | 0.1.0 | `Caliper/pyproject.toml` | package version |
| `caliper.license` | MIT | `Caliper/LICENSE` | clean reuse |
| `caliper.api.conformal` | `caliper.conformal.conformal_offset(scores, alpha)` / `SplitConformalResidual` | `Caliper/caliper/conformal.py` | split-conformal interval offsets |
| `caliper.api.metrics` | `caliper.metrics.empirical_coverage` | `Caliper/caliper/metrics.py` | coverage primitive (cross-check) |

Caliper is in-tree, MIT, and not in review, so the ruler itself is **not** a provisional
risk. If Caliper's conformal offset definition changed, the deployed interval *widths*
would change and CP3 would need a re-run — `invivo.build_deployer` carries a numpy-only
fallback equal to the textbook split-conformal offset so the method is never silently
different.

## 2. GAUGE (PROVISIONAL — baseline + data template)

| key | pinned value | source | role for Echo |
|---|---|---|---|
| `gauge.seed` | `DEFAULT_SEED = 20260613` | `Gauge/gauge/cohort.py` | synthetic determinism |
| `gauge.data.template` | `scripts/fetch_invivo.py`, `gauge/invivo.py` | Gauge | download-on-demand fetch/provenance template |
| `gauge.data.dataset` | ACRIN-6698 TrT0/TrT1, CC-BY-4.0, DOI 10.7937/tcia.kk02-6d95 | `Gauge/results/invivo_real_provenance.json` | the public test–retest data |
| `gauge.baseline.D` | width vs \|ΔADC\| Spearman **r=+0.60**, 95% CI [0.42, 0.72] (n=76) | Gauge paper §4.2.2 | the RANK check Echo must beat/be-distinct-from |
| `gauge.baseline.Dstar` | **r=−0.17** (null, CI [−0.39, 0.05]) | Gauge paper §4.2.2 | identifiability split Echo reproduces as a probe |

**The Gauge assumption:** Gauge's §4.2.2 result and dataset land as submitted. If the
baseline number or the dataset arm changes in revision, Echo's *delta-vs-Gauge* framing must
be re-stated. The *data* itself (ACRIN-6698) is public and stable regardless.

## 3. FASHION (PROVISIONAL — the posterior/ruler the widths derive from)

**Paper status (PINNED):** in review at **NMR in Biomedicine** (retooled,
boundary-railing-first; resubmitted from MRM); no DOI. Source: `Fashion/paper_retool/`,
`Gnomon/handoff/CLAIMS_LEDGER.md`.

| key | pinned value | source | role for Echo |
|---|---|---|---|
| `fashion.version` | 0.1.0 | `Fashion/pyproject.toml` | package version |
| `fashion.zenodo` | 10.5281/zenodo.20649669 | `Fashion/README.md` | citable snapshot |
| `fashion.role` | calibrated IVIM posterior / **scoped, ground-truth-only** ruler that conformal widths are built on | Fashion | the uncertainty whose *scale* Echo checks |

**The Fashion assumption:** the calibrated posterior that defines a "well-scaled" interval
lands as submitted. If Fashion's ruler changes, the *interpretation* of what a correctly
sized interval is changes, and CP3's reading is re-stated.

**Retool note — reinforces Lethe's thesis.** The retooled Fashion now openly owns a
**bounded, conditional limit**: under the honest CRLB its calibration ruler is a
*scoped secondary* (ground-truth/synthetic only — it cannot be applied to a real
scan), and even on ground truth it under-covers D\* *conditionally* in the high-D\*
tercile (honest CRLB 0.63 [0.60, 0.67]); the dramatic marginal 0.30/0.67 severity is
**dropped**. This **reinforces** Echo/Lethe's central claim that ground-truth-free
trust is bounded: the ruler's own paper now agrees its calibration certificate is
scoped and conditional — exactly the kind of bound Lethe records. Echo's own result
lives on a *disjoint axis* (real-data test–retest **precision**, coverage 0.263 vs an
analytic 0.755 target); it never quoted Fashion's coverage figures, so nothing here is
re-numbered — only the framing aligns. No overclaim is introduced: Fashion's bound is
on *synthetic* conditional coverage; Lethe's is on *real* repeatability scale; the two
bounds are independent and reported as such.

## 4. MINOS (PROVISIONAL — the decision/trust lens)

| key | pinned value | source | role for Echo |
|---|---|---|---|
| `minos.role` | decision-value-of-calibration lens | `Minos/` | frames *why* a correctly-scaled interval matters for trust |
| `minos.zenodo` | (pending) | Minos | citation once published |

**The Minos assumption:** Minos's decision-calibration framing lands as submitted; Echo
cites it to motivate the trust claim, not to derive a number.

## 5. DATA SOURCE

- **Synthetic** (`harness.py`, `invivo.synthetic_cohort`) — seeded, open; the SOLID half.
- **Real** ACRIN-6698 / I-SPY2 breast DWI, same-day test–retest (TrT0/TrT1, n≈76) —
  CC-BY-4.0, download-on-demand, **provenance manifest only** (`results/invivo_provenance.json`),
  pixel/array data git-ignored. No `pancData3` / MSK / clinical data in tree or history.

**Default for the SOLID self-test: synthetic only. IP gate passes — nothing private required.**

## 6. Re-validation contract

When Fashion / Minos / Gauge publish (or revise to final):
1. Update the `*.version` / `*.zenodo` / DOI / baseline rows above to the published artifact.
2. Run `bash reproduce.sh` (one command): CP1 self-test, CP2 data check, CP3 validation,
   CP4 consistency.
3. If all green, PROVISIONAL flags may be cleared (see `PROMOTION.md`). If any fails, the
   dependent result is invalidated — fix it before clearing.

Provenance at build time: monorepo HEAD `73d588e`; Caliper/Fashion/Gauge/Minos versions as
pinned above.
