# STUB_LEDGER.md — what each placeholder stands in for, and what moves when it lands

Matrix runs on three labelled placeholders (NOT-Fashion / NOT-Minos / NOT-Forge) behind
fixed interfaces. This ledger records, for each: the current stubbed value/behaviour, exactly
which reported result depends on it, and the **re-check that fires on drop-in** — in
particular whether the headline number moves and whether the **F1 honest negative** is
affected (it is not; the reason is stated). Promotion mechanics are in `PROMOTION.md`; the
publication HOLD is in `RELEASE.md` / `release.json`.

The re-check on every drop-in is one command: wire the real adapter into
`matrix/loop.py :: Interfaces` (one line, never a loop stage) and run `bash reproduce.sh`.

---

## 1. Fashion — the calibration ruler  (release condition; `\PROV`)

| | |
|---|---|
| interface | `matrix/interfaces/ruler.py :: Ruler.calibrate(mu, raw_sigma, truth=None)` |
| placeholder | `PlaceholderRuler` — rescales the raw posterior spread by a per-parameter factor so the nominal interval hits nominal empirical coverage on the twin; reports coverage/ECE the way `uq.calib` does |
| drop-in | `FashionRulerAdapter` → `uq.bayesian.mcmc_uncertainty` + `uq.calib.coverage/ece` |
| **reported results that depend on it** | ruler coverage(f) `0.944 [0.903, 0.979]`, ECE_f `0.007`; and, through the calibrated σ_f it produces: the trust-gate firing (AUROC `1.00`), the action thresholds, and therefore **every CP4/Ferry suppression + convergence number** (gated `0.000` / ungated `0.422` / trusted `0.573`; trusted drop `0.176`; grounded drop `0.169`) |
| **on drop-in: does the headline move?** | **Possibly.** A real ruler may calibrate σ_f differently from the placeholder's variance-matching, which would shift coverage/ECE and (via σ_f) the suppression/convergence magnitudes. The interface is fixed, so the loop and gates are unchanged; only the numbers re-derive. |
| **F1 affected?** | **No.** F1 (held untrusted f-drop under real dose) depends on real delivered dose × the twin's response model, not on how σ_f is calibrated. The held voxels drop because dose was already delivered, regardless of the ruler. |
| **"loop closes" affected?** | **No** — SOLID; independent of any consumed component. |

## 2. Minos — the trust gate + action gate  (release condition; `\PROV`)

| | |
|---|---|
| interface | `matrix/interfaces/gates.py :: TrustGate.trustworthy`, `ActionGate.act` |
| placeholder | `PlaceholderTrustGate` (untrustworthy iff σ_f too wide, i.e. the low-SNR zone) + `PlaceholderActionGate` (Bayes-action-shaped interval rule on f; trust gate forces ESCALATE on untrustworthy voxels) |
| drop-in | `MinosGateAdapter` → `minos.gate.{gate_fires,votg,gated_actions}`, `minos.decision.bayes_action` |
| **reported results that depend on it** | trust-gate firing (AUROC `1.00`, fire `1.00`/`0.02`); suppression (gated `0.000`, ungated `0.422`, trusted `0.573`); convergence (drop `0.176`, TREAT `21→…→0`); the grounded suppression/convergence numbers |
| **on drop-in: does the headline move?** | **Possibly in magnitude, not in kind.** The trust-gate *concept* (suppress action where the error bar is untrustworthy) is exactly Minos's VoTG, so a faithful drop-in reproduces the qualitative behaviour; the exact rates may shift as the real gate's thresholds differ. |
| **F1 affected?** | **No (qualitatively); magnitude may shift slightly.** F1 needs only that untrusted voxels are *held* (no new action) yet still receive standing dose. A real Minos gate still holds the untrusted set (same VoTG), so held voxels remain exposed to delivered dose and still drop — the sign and CI-excludes-0 survive. The precise `0.148` could move if the trusted/untrusted partition changes, but the negative does not. |
| **"loop closes" affected?** | **No** — SOLID. |

## 3. Forge — the dose engine  (**NOT a release condition**; NOT-BUILT, deferred 2027)

| | |
|---|---|
| interface | `matrix/interfaces/dose.py :: DoseEngine.replan(current_dose, action, state, cfg)` |
| placeholder | `PlaceholderDoseEngine` — analytic prescription update (TREAT→boost toward target, SPARE→de-escalate, ESCALATE→hold); no geometry, no transport, no Electron Return Effect |
| drop-in | `ForgeDoseAdapter` → Forge's Monte-Carlo engine **when it exists** (deferred 2027) |
| **reported results that depend on it** | only the *physical realism* of the re-plan, and the **F2 placeholder-warrant artefact** (strict "TREAT ⇒ dose strictly increases" trips on a non-uniform real prescription / on re-TREAT) |
| **on drop-in: does the headline move?** | The F2 artefact resolves (a geometry-aware engine boosts relative to local dose). No headline calibration/trust/convergence number depends on the dose engine's internals. |
| **F1 affected?** | **No.** F1 is about dose *already delivered* by the standing prescription; the re-plan engine governs *future* dose only. |
| **release impact** | **None.** Forge is explicitly a non-condition in `release.json` — it never holds submission. Presented as drop-in future work. |

---

## Forward-citation finalization checklist (per DOI)

Both Fashion and Minos are cited as in-review forward-refs (`\bibitem{fashion}`,
`\bibitem{minos}` in `paper/matrix.tex`), with no fabricated DOI. On publication of each:

1. Set its flag `satisfied: true` and fill its real `doi` in `release.json`.
2. Update its pin (commit/version/DOI) in `ASSUMPTIONS.md` (§Fashion / §Minos).
3. Replace the `% PROVISIONAL-{fashion,minos}` `\bibitem` placeholder with the published
   citation + DOI; drop the `\PROV` markers on the numbers that rode on it **only after**
   step 5 confirms they still hold.
4. `python release_gate.py status` — both satisfied → RELEASE; the `SUBMISSION_HOLD` marker
   clears.
5. `bash reproduce.sh` with the real adapter wired in (per `PROMOTION.md`) — every gate must
   stay green. If a number moved, update it (it is regenerated from the seeded gate); if the
   F1 negative changed sign, stop — that would be a real finding, not a number to paper over
   (it should not, for the reasons above).

`Forge` needs no citation finalization: it is future work, not a dependency met.

**Central-finding survival (summary):** the F1 honest negative and "the loop closes" survive
all three drop-ins; only secondary magnitudes (coverage/ECE, suppression/convergence rates)
may re-derive when Fashion/Minos land, and the F2 artefact resolves when Forge lands.
