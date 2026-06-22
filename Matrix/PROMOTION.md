# PROMOTION.md — swapping a placeholder for the real component

> **Documented, not executed.** This records exactly how each consumed component's labelled
> placeholder (NOT-Fashion / NOT-Minos / NOT-Forge) is replaced by the real component once it
> lands, and how Matrix is re-validated against it. **None of this is run now** — Matrix runs on
> placeholders until a real component is available and its gate is green.

## The one rule

The loop (`matrix/loop.py` stages) is written against three interfaces and **never** references a
concrete component. Promotion is therefore always the same two-step move, and it **never edits a
loop stage**:

1. Implement an adapter behind the existing interface (same method signature).
2. Point `matrix/loop.py :: Interfaces` at the adapter instead of the placeholder.

`verify_cp3.py` check 2 already proves this is sufficient: it swaps a different dose engine into
`Interfaces.dose_engine` and the full loop runs unchanged.

## Per-component promotion

| component | placeholder (now) | adapter (to build) | interface (unchanged) |
|---|---|---|---|
| **Fashion** ruler | `PlaceholderRuler` | `FashionRulerAdapter` → `uq.bayesian.mcmc_uncertainty` + `uq.calib.coverage/ece` | `Ruler.calibrate(mu, raw_sigma, truth=None)` |
| **Minos** gates | `PlaceholderTrustGate`, `PlaceholderActionGate` | `MinosGateAdapter` → `minos.gate.{gate_fires,votg,gated_actions}`, `minos.decision.bayes_action` | `TrustGate.trustworthy`, `ActionGate.act` |
| **Forge** dose | `PlaceholderDoseEngine` | `ForgeDoseAdapter` → Forge MC engine (deferred 2027) | `DoseEngine.replan(current_dose, action, state, cfg)` |

## Precondition gate (must be green before clearing a PROVISIONAL flag)

1. The real component is available (published / built) and its pins in `ASSUMPTIONS.md` are updated.
2. Its adapter is wired into `Interfaces` in place of the placeholder.
3. `bash Matrix/reproduce.sh` is **all-green** with the adapter in place — in particular CP2/CP4 still
   close and behave sensibly with the *real* ruler/gates, and CP3 with the *real* dose engine.

If a gate fails, the real component genuinely changed a dependent result. Fix the result; do not
promote a stale number.

## Flags that clear

- Remove the **PROVISIONAL** flag for that component from `ASSUMPTIONS.md` §0 and update its row to
  the published artifact (DOI / tag / commit).
- Replace the placeholder's `label` ("NOT-Fashion …") wherever it is surfaced (gate provenance,
  `RESULTS_CP4.md`) with the real component's citation.
- Matrix's **scope ceiling does not change**: even with all three real components wired in, Matrix
  remains a **synthetic-twin** harness. A clinical claim requires a real scanner + real data
  (Keystone's real-time/offline modes, or the Ferry grounding on a public RT dataset) — explicitly
  **out of scope** here.

## Ferry — the substrate swap (BUILT; real anatomy + dose geometry)

Ferry is the **fourth swap point** and the only one already grounded on real data. Unlike
Fashion/Minos/Forge (which swap a *consumed component*), Ferry swaps the **substrate** — the
anatomy/dose twin — behind the same one rule: the loop never references a concrete twin, so
`run_iteration(twin, …)` takes a `GroundedTwin` and **`loop.py` is untouched** (proven
byte-for-byte in `verify_ferry_cp1.py`).

| swap | placeholder/default | adapter (built) | interface (unchanged) |
|---|---|---|---|
| **Ferry** substrate | synthetic `Twin.build(cfg)` | `GroundedTwin.from_substrate` ← TCIA RTSTRUCT+RTDOSE (`matrix/ferry/`) | the `Twin` contract consumed by `run_iteration` |

**What Ferry promotes, and what it does NOT.** Ferry promotes the substrate from *synthetic
geometry* to *real anatomy + real dose geometry* — clearing the "shown only on a pure synthetic
twin" objection. It does **not** clear the perfusion/IVIM ceiling: perfusion stays synthetic
(no scanner). So a Ferry-grounded result clears exactly one flag — *geometry is real* — and the
**real-IVIM / clinical flag stays set** until a scanner study lands. Update `ASSUMPTIONS.md` §6
pins (dataset/DOI/license/patient) if the grounding dataset changes; the synthetic-twin gates
(CP1–CP4) are unaffected and must stay green.
