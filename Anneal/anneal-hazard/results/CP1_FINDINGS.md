# CHECKPOINT 1 — engine validated, but a regime problem blocks the experiment

The integrator, order-parameter computation, and death detector all work and are
dt-converged. But validation surfaced a fundamental mismatch between the **model**
(two-population mean-field Sakaguchi–Kuramoto) and the **phenomenology the experiment
needs** (a metastable chimera that dies within T_max, with lifetime growing in N).

## Engine is correct and converged
Single run, A=0.50, N=128, β=0.10, seed=20260609 (`cp1_single_run.png`):
- dt vs dt/2: Δτ = 0.25 (0.2% of τ); Δ mean r_incoh = 4.3e-4 (0.05%). RK4 is well resolved.
- Mean-field reduction re-derived from the all-to-all identity; matches the spec exactly.

## Problem 1 — at β=0.10 there is NO chimera (it just synchronizes)
Trajectory scan, N=256, β=0.10, seed=20260609 (`diag_traj.png`): for **every** A in
{0.10,0.15,0.20,0.25,0.30,0.40} both populations fully synchronize (r1=r2≈1.000,
std 0.000) within ~100 t.u. No sustained incoherent plateau. β=0.10 is above the
chimera region.

## Where the chimera actually lives
(β,A) map, late-time mean min(r1,r2), N=256 (`*` = sustained intermediate r = chimera):

```
beta\A   0.05   0.10   0.20   0.30   0.40   0.50
0.02     1.00   0.82*  0.66*  0.73*  0.87*  1.00
0.05     1.00   1.00   1.00   0.74*  0.89*  1.00
0.08     1.00   1.00   0.68*  1.00   1.00   1.00
0.10     1.00   1.00   1.00   1.00   1.00   1.00
0.15+    1.00   1.00   1.00   1.00   1.00   1.00
```
Clean chimeras (r_incoh ≈ 0.67, the canonical Abrams value) appear only at small β
(≲0.08), around A≈0.2 (`diag_chimera.png`).

## Problem 2 — in the chimera band the chimera is STABLE (no deaths to analyze)
Metastability scans (5–6 seeds, T_max=2000, ε=0.03, hold=50):
- β=0.07, A=0.20, N ∈ {6,8,10,12,16,32,48,64,96,128,256}: essentially **0 deaths** —
  every realization survives to T_max. Even N=6 survives. Lifetime ≫ 2000 at all N.
- Approaching the A-boundary (β=0.05, A: 0.30→0.22) and the β-edge (A=0.20,
  β: 0.092→0.098): behavior is strictly **bimodal** — either eternal survival (censored
  at 2000) or near-instant death (τ ≈ 40–100). The fast deaths are **N-independent**,
  i.e. IC basin-misses (chimera never forms), NOT finite-N escapes from an established
  chimera. There is no window with median τ in [200,1800] or lifetime growing in N.

## Diagnosis
The two-population mean-field chimera is a genuine **stable attractor** within its
existence region (Abrams–Mirollo–Strogatz–Wiley 2008). Its collapse to sync is a
deterministic basin/bifurcation phenomenon, not a finite-N chaotic transient. The
"chimera death = chaotic transient whose lifetime grows ~exponentially with N"
phenomenology — the thing this survival/hazard study is built to measure — is the
hallmark of the **nonlocally-coupled ring** model (Wolfrum & Omel'chenko 2011), a
different system. As specified, there is no chimera-death hazard to measure here.

## Decision required (see chat) — candidate paths
1. Switch to the nonlocal-ring chimera (canonical lifetime-grows-with-N system;
   still O(N) via cosine-kernel order-parameter reduction). [recommended]
2. Keep the two-population model but add weak phase noise → Kramers escape with a
   well-defined, tunable hazard (memoryless vs structured becomes a real question).
3. Keep two-population deterministic; pivot to mapping the chimera→sync boundary and
   report the structural null (no finite-N hazard).
