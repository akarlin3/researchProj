# Huang coverage map — every critique → where/how the retool resolves it

Reviewer "Huang" returned the earlier Fashion manuscript at *Magnetic Resonance
in Medicine* on **methods**, in three buckets: **internal inconsistencies**,
**incompleteness** (under-specified dataset IDs, training/fitting detail, the
Cramér–Rao approximation, the SD conventions), and **overextended claims** (a
calibration statement is only as trustworthy as the simulator behind its
reference truth; single-subject generalization; the dramatic marginal headline).

This table maps each complaint to the exact place in the retooled,
boundary-railing-first manuscript (`manuscript.tex`) and the claims ledger
(`Gnomon/handoff/CLAIMS_LEDGER.md`) where it is resolved. It is both the internal
completeness check (CP3 gate) and the factual basis for the cover letters (CP4).

Section references are to `manuscript.tex`; ledger codes (K/R/D) are the
disposition from the Gnomon CP4 hand-off.

---

## A. Internal inconsistencies

| # | Huang complaint | How the retool resolves it | Where | Ledger |
|---|---|---|---|---|
| A1 | The marginal $D^*$ under-coverage headline (0.30 Laplace / 0.67 MCMC-SD) is not internally reproducible; it depends on undocumented choices. | **Dropped as a headline.** An independent clean-room rebuild gives an *undramatic* honest marginal (0.80 / 0.90). The dramatic figure is shown only as a labelled reconstruction under a named floored SD convention, never asserted as a result. The kept, reproducible statement is the per-tercile conditional table. | §Results/“The convention reconstructs…”; Table (conditional coverage) | **D1** drop; **R1** reframe |
| A2 | The quantile "fix" is presented as fully restoring coverage, inconsistent with a surviving high-$D^*$ failure. | The fix is scoped to the **marginal** level (0.90); the **conditional** residual high-$D^*$ wall (0.81) is stated explicitly in the same paragraph — no contradiction. | §Results/secondary ruler ("interval shape…") | **K2** keep, **R2** reframe |
| A3 | "Below-CRLB-floor overconfidence" is stated alongside an honest CRLB without reconciling the two. | Reconciled: it is the **same axis as the SD convention**. The floored convention reconstructs pooled 0.68 ≈ the original 0.67; the honest CRLB does the opposite. Both are documented in one place. | §Methods/CRLB(b); §Results/convention | **R3** reframe |

## B. Incompleteness

| # | Huang complaint | How the retool resolves it | Where | Ledger |
|---|---|---|---|---|
| B1 | Dataset identifiers under-specified — which data, which ROI, what defined the "1618 high-SNR voxels." | Every dataset named with a stable ID: OSIPI Zenodo **14605039** (CC-BY-4.0, sha256-verified), homogeneous mask = **1932 voxels**; TCGA-LIHC DOI **10.7937/K9/TCIA.2016.IMMQW8UQ**. The load-bearing number is the **rate** with its selection sensitivity, not the unreproducible 1618 count. | §Methods/Datasets | K1 |
| B2 | Training/fitting detail missing (NLLS box/init/solver; MCMC sampler; the amortized network's architecture and training). | Stated in full: NLLS box/init/solver/railing threshold; MCMC burn 1500 / keep 2000 / thin 2 / seed; NPE 5-layer flow, 80k sims, 40 epochs, batch 512. | §Methods/Forward model and estimators | K2, K3 |
| B3 | The Cramér–Rao (CRLB) approximation is not stated. | Stated explicitly: the Gaussian/CRLB covariance is weakest exactly where $D^*$ is most skewed (low SNR / high $D^*$); the MCMC results do not rely on it. | §Methods/CRLB(a) | R1 |
| B4 | **The railed-voxel SD convention — the headline-driving choice — is undocumented.** | Both conventions are documented: **honest CRLB** (wide where uninformative; the default) vs **overconfident floor** (narrow; "overconfident by design"). The convention that reconstructs the original headline is named and quantified. | §Methods/CRLB(b); completeness checklist (Table) | R1, R3 |

## C. Overextended claims

| # | Huang complaint | How the retool resolves it | Where | Ledger |
|---|---|---|---|---|
| C1 | The calibration ruler is overextended: empirical coverage is only as trustworthy as the assumed noise/forward model behind its reference truth. | **The contribution is re-aimed.** The primary claim is now boundary-railing, read directly off the optimiser with **no ground-truth or noise-model assumption** — it cannot be overextended. The ruler is demoted to a scoped secondary, explicitly undefined on real scans. | §Introduction; §Results/primary; §Results/secondary ruler | K1 primary |
| C2 | Single-subject in-vivo generalization is overclaimed. | No population claim is made from the single open subject; a within-subject voxel statistic with a voxel bootstrap is reported, and generality is supported by **cross-cohort agreement** (full ROI + independent TCGA-LIHC liver, an unrelated organ/scanner/site). | §Results/primary; §Discussion/Limitations | K1 |
| C3 | The dramatic marginal coverage severity is overextended. | Reported **conditionally** (per true-$D^*$ tercile) under the honest CRLB; the marginal headline is dropped. | §Results/secondary; Table | D1, R1 |
| C4 | Claims are made beyond the evidence (OOD self-consistency gate; estimator timing; in-vivo brain held-out-$b$). | These are **named as not claimed** — explicitly out of scope of the evidence presented and left to future work. | §Discussion/"What is and is not claimed" | OUT OF SCOPE (not dispositioned) |

---

## Internal check (CP3 gate)

- **Every Huang bucket is covered:** A (3/3 inconsistencies), B (4/4 completeness
  items, incl. the decisive SD convention), C (4/4 overextension points).
- **No claim outside the ledger:** every cell above maps to a KEEP/REFRAME/DROP
  code or to an explicitly not-claimed OUT-OF-SCOPE item.
- **The decisive item (B4)** — the railed-voxel SD convention Huang flagged as
  missing — is now documented in Methods §CRLB(b) *and* shown to be the lever that
  reconstructs the original headline, which simultaneously closes A1 and A3.
- **Numbers-gate** (`consistency.py`) passes: 126 macros, every number traced to a
  seeded Sextant/Gnomon result JSON; the dropped 0.30/0.67 marginal is guarded
  against (asserted absent as an honest coverage macro).
