# Retool hand-off — Gnomon → the retooled Fashion paper

This is the package the retooled Fashion manuscript consumes. It is **spine-agnostic**:
it commits to validated *components* and a *claims ledger*, **not** to which claim
leads. It is equally usable for a **ruler-first** paper (the quantile/flow calibration
story leads) or a **boundary-railing-first** paper (the real-data identifiability
signature leads). Choosing the spine — and the actual redraft (incorporating Sextant) —
is the next, separate step.

**Verdict (CP3): PARTIAL.** The diagnostic and the mechanism reproduce; two *marginal*
headline numbers do not, for a documented and resolved reason (the railed-voxel SD
convention). Nothing here is asserted without a seeded run with bootstrap CIs.

---

## 1. Contents

| Artifact | What it is |
|---|---|
| **Clean reference implementation** — `gnomon/` | Independent rebuild (no Fashion/Caliper code): `forward.py`, `nlls.py` (+railing), `bayes.py` (Laplace + MCMC), `flow.py` (MAF/NPE), `metrics.py` (ruler), `bootstrap.py`, `osipi.py`, `reproduce.py`, `reframe.py`. |
| **Claims ledger** — [`handoff/CLAIMS_LEDGER.md`](handoff/CLAIMS_LEDGER.md) | Every Fashion claim → KEEP / REFRAME / DROP → run-evidence + CI → reworded text; primary-claim candidates flagged (no order fixed). |
| **Reframed conditional table** — [`handoff/conditional_coverage.json`](handoff/conditional_coverage.json) (+ §3 below) | Per-D\*-tercile D\* coverage for Laplace SD & MCMC SD under **both** SD conventions, bootstrap CIs. The honest replacement for the marginal 0.30/0.67. |
| **Complete methods** — [`docs/METHODS.md`](docs/METHODS.md) | Every Huang-flagged item incl. both SD conventions; completeness checklist §10. |
| **Keep-set + verdict** — [`results/reproduction.json`](results/reproduction.json), [`VERDICT.md`](VERDICT.md) | The reproduced numbers (T1/T3c/T4) and the CP3 verdict. |
| **Divergence appendix** — §4 below + [`scripts/divergence_diagnostic.py`](scripts/divergence_diagnostic.py) | Why the headline changed, documented and resolved. |

## 2. Keep-set (reproduced; carry as-is, with CIs)

- **K1 — NLLS D\* boundary-railing, real OSIPI abdomen:** **54.2% [52.0, 56.4]**
  (vs claimed 54.7%); 56.2% at rail_tol 1e‑2; 58.7% for b0‑SNR>25.
- **K2 — quantile interval restores *marginal* D\* coverage:** **0.90 [0.89, 0.91]**
  (D 0.95, f 0.96). The interval *shape*, not a wider SD, is the fix.
- **K3 — MAF flow beats railed NLLS:** D\* coverage **0.98 vs 0.76**, ECE **0.069 vs
  0.121**, sharpness **0.112 vs 0.181**; all three gaps' CIs exclude 0.

## 3. The reframe — per-true-D\*-tercile D\* coverage (central 0.95), both conventions

Replaces the marginal 0.30 / 0.67. Honest CRLB is **recommended** (see §4); the floored
convention is shown only to explain the original headline.

| Estimator · convention | low D\* | mid D\* | **high D\*** | pooled |
|---|---|---|---|---|
| Laplace SD · **honest** | 0.91 [0.89,0.94] | 0.86 [0.83,0.89] | **0.63 [0.60,0.67]** | 0.80 [0.78,0.82] |
| Laplace SD · floored | 0.87 [0.85,0.90] | 0.75 [0.72,0.78] | 0.41 [0.37,0.45] | 0.68 [0.65,0.70] |
| MCMC SD · **honest** | 0.95 [0.93,0.97] | 0.95 [0.94,0.97] | **0.81 [0.78,0.84]** | 0.90 [0.89,0.92] |
| MCMC SD · floored | 0.92 [0.90,0.94] | 0.85 [0.83,0.88] | 0.65 [0.61,0.69] | 0.81 [0.79,0.83] |
| MCMC **quantile** (recommended) · honest | 0.93 [0.90,0.95] | 0.97 [0.95,0.98] | 0.81 [0.78,0.84] | 0.90 [0.89,0.91] |

Reading: under the honest CRLB a **real, reproducible** under-coverage survives in the
**high-D\*** tercile (Laplace 0.63, MCMC-SD 0.81); even the quantile fix leaves a
residual high-D\* gap (0.81) — the identifiability wall. The dramatic marginal severity
is a property of the floored convention, not of the data.

## 4. Divergence appendix — why the headline changed (documented & resolved)

Fashion's marginal Laplace 0.30 / MCMC-SD 0.67 did not reproduce (clean rebuild: 0.80 /
0.90). Two undocumented choices reconstruct them, both **run** not asserted
(`scripts/divergence_diagnostic.py`, `gnomon/reframe.py`):

1. **Cohort regime** — pooling over a prior-spanning cohort dilutes a failure that is
   concentrated in high D\* (honest Laplace pooled 0.80 vs high-tercile 0.63).
2. **Railed-voxel SD convention** — flooring an unidentified D\*'s SD ("overconfident by
   design") drops honest Laplace pooled 0.80 → **0.68** ≈ Fashion's MCMC-SD 0.67. The
   honest CRLB does the opposite (wide where uninformative).

**Resolution:** report under the **honest CRLB**, state the convention explicitly
([`docs/METHODS.md` §5](docs/METHODS.md)), and present coverage **conditionally**
(§3). The finding stands; only the manufactured marginal severity is removed.

## 5. Using this package under either spine

- **Boundary-railing-first:** lead with **K1** (real-data 54.2% railing) as the
  identifiability signature; the conditional table (§3) and the flow superiority (K3)
  become the calibration consequence. Sextant's in-vivo railing replications slot in
  alongside K1.
- **Ruler-first:** lead with **K2/K3** (shape-correct/amortized intervals fix marginal
  calibration) as the method; K1 + §3 become the motivating failure and its
  conditional residual.

The package fixes neither order. The claims ledger marks K1 *and* K2/K3 as primary-claim
candidates; pick the lead downstream.

## 6. Reproduce every number — one command

```bash
KMP_DUPLICATE_LIB_OK=TRUE PYTHONPATH=Gnomon python Gnomon/scripts/build_handoff.py
#  -> results/reproduction.json (keep-set + verdict) and
#     handoff/conditional_coverage.json (reframe, both conventions)
# gates only:  python -m pytest Gnomon/tests -q     # 16/16
```

All hand-off prose numbers trace to those two JSON files (seed 20260620).

## 7. Clean-IP (unchanged)

Synthetic substrate = **Lattice** (read-only sibling import); real data = **OSIPI**
Zenodo 14605039 (CC-BY-4.0), download-on-demand, gitignored, sha256-verified. No
proprietary/clinical data (no `pancData3` / MSK) in tree or history. **Caliper's ruler
import remains forbidden** and is enforced (`gnomon/_paths.py` lattice-only + static
AST test). Gnomon is **not** a standalone paper — it is the clean technical core that
feeds the retool.
