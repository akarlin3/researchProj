# Exploratory verdict — agent-derived TGKP slaved-cumulant scan (NOT the CP3 verdict)

**Register:** agent-derived closure, authorized downgrade (2026-06-10).
Protocol frozen at commit `fcf94b0` *before* any scoring run; interpretation
rules R1–R4 in `run_exploratory_tgkp.py`. The formal CP2/CP3 path (human
derivation, `f_corr.py`, `run_probe.py`) remains open and ungated by this
result. Manuscript untouched.

## Verdict: the bet did not land in this register (R1 — exploratory clean negative)

The TGKP slaved-two-cumulant drift correction, at data-anchored effective
noise intensity (σ² ∈ [4.4e-4, 4.4e-3] s⁻¹, from measured σ_HF·√N ≈ 0.047
and declared γ_R ∈ [0.1, 1] s⁻¹), produces a real, correctly-signed,
**cleanly N-independent** prolongation of **3–10%** — a factor of ~30 short
of the measured 3.2×.

| variant | σ² (s⁻¹) | per-N factor (8/16/32/64) | median | CV | verdict |
|---|---|---|---|---|---|
| V1 central | 1.33e-3 | 1.032 / 1.031 / 1.013 / 1.038 | 1.032 | 0.009 | FAIL (c1) |
| V2 γ_R low | 4.4e-4 | 1.010 / 1.017 / 1.246 / 1.016 | 1.016 | 0.094 | FAIL (c1) |
| V3 γ_R high | 4.4e-3 | 1.103 / 1.094 / 1.099 / 1.129 | 1.101 | 0.012 | FAIL (c1) |
| V4 per-N measured | 0.97–1.78e-3 | 1.023 / 1.029 / 1.013 / 1.051 | 1.026 | 0.014 | FAIL (c1) |
| V5 V1 + CLT noise | 1.33e-3 + c=0.05/√N | 0.941 / 0.989 / 1.027 / 0.978 | 0.984 | 0.031 | FAIL (c1) |

Conditions 2–4 (N-independence, breath-phase locking, rising k_cyc) **pass in
every variant**; only condition 1 (factor in [2.9, 3.5]) fails, everywhere.
The pre-registered secondary measured-pattern criterion fails in every
variant (S1 factor pattern off by ~3×). Gates: OA-forcing transcription
exact to 8.9e-16 vs `rhs_3d`; zero-correction determinism 0.23–0.25%;
harness previously shown to reject the known-wrong additive mechanism.

## What the scan establishes

1. **The mechanism class has the right *shape* but not the right *size*.**
   An N-free deterministic drift drag on the collective flow reproduces the
   phenomenology's signature effortlessly: N-independence (CV 0.009–0.014,
   far cleaner than any excluded stochastic mechanism), preserved phase
   locking, preserved rising hazard. Magnitude is the sole failure.
2. **Quantified gap.** V3 implies the realized drag on the corner escape
   rate (σ_spiral = +0.0124 s⁻¹) is ≈ σ²/4. Reaching 3.2× needs drag
   ≈ 0.69·σ_spiral, i.e. σ² ≈ 3e-2 s⁻¹ — roughly **25× the central
   data anchor** (equivalently c ≈ 0.2 vs measured 0.047, or γ_R ≈ 10 s⁻¹
   vs the declared 0.1–1 bracket). Per rule R3 no such variant was added;
   this estimate is reported, not run.
3. **The stochastic CLT part is confirmed irrelevant-to-harmful** (V5
   slightly shortens, consistent with Appendix B).

## Open-problem sharpening (draft, exploratory register)

> We additionally probed the one remaining principled candidate, a
> finite-size kinetic correction to the collective flow, using an
> exploratory two-cumulant (circular-cumulant) closure of the
> Tyulkina–Goldobin–Klimenko–Pikovsky type with the effective
> per-oscillator intensity anchored to the measured order-parameter
> fluctuations. At this order and anchoring, the correction is a
> deterministic, N-independent drag on the collective flow that reproduces
> the *qualitative* signature of the anomaly — N-independent prolongation
> with preserved breath-phase locking and rising cycle hazard — but at
> 3–10% magnitude, a factor of ~30 below the observed 3.2×; matching the
> observation would require an effective intensity ~25× the measured
> anchor. The discrepancy thus survives the kinetic-drift correction at
> physical order, and sharpens to: what mechanism supplies an
> N-independent, order-unity slowdown of the collective escape that
> standard system-size expansions render only at the few-percent level?

(Caveats for any human use of this paragraph: agent-derived closure;
slaved-κ truncation, O(σ⁴); γ_R amplitude-to-rate conversion bracketed,
not derived; dynamical-κ and ghost-specific boundary treatments untested.)

## Files

`kinetic_results/exploratory_tgkp/`: `runs_V*.jsonl` (5 × 800 runs, fixed
seeds), `score_V*.json` (frozen scorer output incl. secondary criterion),
`exploratory_tgkp_summary.json` (gates + verdicts). Regenerate:
`python3 tools/kinetic-probe/run_exploratory_tgkp.py` (~23 min, 9 workers).
