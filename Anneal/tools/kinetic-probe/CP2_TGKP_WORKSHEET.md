# CP2 worksheet — TGKP circular-cumulant route (basis chosen 2026-06-10)

The human (Avery) has selected the Tyulkina–Goldobin–Klimenko–Pikovsky
circular-cumulant expansion as the derivation basis for `f_corr`. This
worksheet pins down **exactly what remains to be supplied** and what the
harness will then do mechanically. It contains no candidate equations for our
system: the slots below are questions, not proposals, and the harness authors
will not fill them.

## What the published basis provides (verbatim, for orientation)

TGKP, PRL 120, 264101 (2018) (arXiv:1804.05326), Eq. (13) — per population,
with Z the OA order parameter, κ the second circular cumulant, h the complex
forcing (for Kuramoto-type coupling h is a linear combination of the
populations' Z's), σ² the intrinsic-noise intensity, h* = conj(h):

    Ż = (iΩ₀ − γ)Z + h − h*Z² − σ²Z − h*κ
    κ̇ = 2(iΩ₀ − γ)κ − 4h*Zκ − σ²(4κ + 2Z²)

Scalings: κ₂ ∼ σ²κ₁², |κ_m| ∼ σ^{2(m−1)}, truncation accuracy O(σ⁴). The
paper demonstrates the scheme on the two-symmetrically-coupled-population
chimera system (Abrams et al. class) with *injected* noise, and the closure
that makes the two-cumulant truncation a rigorous first-order beyond-OA
correction (for small noise) is Goldobin et al., Chaos 28, 101101 (2018).

## The four slots the human must fill (the derivation)

**Slot 1 — the effective noise intensity σ²_eff, and its N-dependence.**
The target system is deterministic with identical oscillators; TGKP's σ² is
intrinsic-noise intensity. What replaces it at finite N, derived from what,
and with what N-scaling? This is the load-bearing choice for the bet: the
measured prolongation is N-independent, and the no-tuning rule applies in
full — σ²_eff must be theory-fixed (no freedom), or declared as a range in
`F_CORR_META["undetermined_coefficients"]` and scored as variants with no
value selected.

**Slot 2 — the two-population forcing terms in our conventions.**
The forcing h entering each population's (Ż, κ̇) equations, written for our
coupling and lag conventions (μ = (1+A)/2, ν = (1−A)/2, α = π/2 − β, ω = 0,
A = 0.5, β = 0.05), and the same for the second population (Y, ν cumulant).
State whether the published two-population demo's form carries over verbatim
or is modified.

**Slot 3 — cumulant treatment: slaved or dynamical.**
Option A: adiabatically eliminate κ (slave it to Z via the σ²-scaling),
yielding a correction that is a pure function of the instantaneous collective
state — this fits the existing `f_corr(rho1, rho2, psi, N) -> (drift, B)`
interface directly. Option B: keep κ (and ν) as dynamical variables — the
state is then augmented beyond (ρ₁, ρ₂, ψ), and the harness driver will be
extended to integrate the augmented system (a mechanical change on our side;
say the word and specify the auxiliary initial conditions). Which one the
derivation justifies, and at what order, is the human's call.

**Slot 4 — stochastic vs deterministic character, and the boundary.**
Does the derived correction carry a stochastic term (a B matrix) in addition
to the deterministic drift, or is it drift-only (B = 0)? If stochastic:
Itô or Stratonovich (the harness driver is Euler–Maruyama, i.e. Itô —
conversion drift, if any, must be stated). And: is reflection at ρ ∈ [0,1]
(the Appendix B operator) the correct boundary treatment near the capture
boundary / homoclinic ghost for this derivation, or does the theory demand
something else?

## What the harness side will then do (mechanical, no judgment)

1. If the supplied form is in complex per-population variables (Z, Y, κ, ν),
   transcribe to (ρ₁, ρ₂, ψ) by Itô change of variables — shown back to the
   human for sign-off before any run (the Itô correction terms are where
   transcription errors hide).
2. Implement as `f_corr` (or the augmented-state driver, per Slot 3), with
   the zero-correction and bit-equivalence gates re-asserted.
3. Run CP3 once: `run_probe.py`, 4 N × 200 realizations, fixed seeds, scored
   against the frozen primary conditions and the pre-registered secondary
   measured-pattern criterion. No coefficient adjusted at any point.

## Standing context (from the verified literature sweep)

Every verified published correction in this class is N-dependent (LITERATURE.md);
for the bet to land, the supplied derivation must make the N-scaling cancel in
the prolongation factor itself. "The expansion yields no N-independent term at
physical order" is a valid Slot-1 answer and will be recorded as the clean
negative per CP3 path 2.
