# FIGURES_REPORT — Chaos (AIP) manuscript Figures 1 & 7

Two publication figures composed from existing merged work (the absorption re-campaign traces and the reduced-ODE bifurcation machinery). New code only under `tools/paper-figures/`; every artifact is deterministic and regenerable from `figures.config.json`. Prior `*_results/` are read-only.

## Figure 1 — annotated example trace

- **Chosen run:** N=16, A=0.5, β=0.05, **seed 100015** (deterministic; selected by CP0 from 85 N=16 candidates and pinned in config).
- **Labels (regenerated, reproduce the campaign bit-for-bit):** t_graze = 83.2 s, t_abs = 221.9 s, n_grazes_before_abs = 2, T_b = 20.5 s (8 breath peaks).
- **Determinism:** the trace is re-integrated from the logged seed with the shipped-identical RK4; the labeler returns t_graze=83.2 / t_abs=221.9, matching the absorption campaign row (83.2 / 221.9). ✔
- **Supervisor over-trigger:** the shipped 2-s detector fires 6 times on this single trajectory — at the first graze and at every subsequent graze that self-heals — before the absorbing crossing.
- **CP0 top-3 candidates (clarity-scored):**

  | seed | t_graze (s) | t_abs (s) | n_grazes | T_b (s) | reform dip | score |
  | --- | --- | --- | --- | --- | --- | --- |
  | 100015 ←chosen | 83.2 | 221.9 | 2 | 20.5 | 0.264 | 8.34 |
  | 100070 | 69 | 236.9 | 2 | 26 | 0.312 | 8.11 |
  | 100023 | 130.2 | 228.3 | 1 | 28.1 | 0.312 | 8.02 |

## Figure 7 — (β, A) stability diagram

- **Homoclinic self-check (β=0.05):** A_hc = 0.4096 vs expected 0.4096 (|Δ| = 9.55e-07, tol 0.002) → **PASS**. Bracket [0.4095655121654199, 0.4096325773458567].
- **Series vs numeric (max deviation over the check β):** saddle-node (Eq. 17) 1.36e-03; Hopf (Eq. 18) 8.71e-04. The series curves are drawn as the SN/Hopf boundaries and the numeric fold / trace-zero locations overplotted as open check marks; they agree to <1.5e-3 in A.
- **Homoclinic curve:** 12 β values in [0.03, 0.2] by escape-to-sync bisection (no degradation up to β=0.20).

  | β | A_hc | bracket width |
  | --- | --- | --- |
  | 0.030 | 0.4405 | 6.7e-05 |
  | 0.040 | 0.4244 | 6.7e-05 |
  | 0.050 | 0.4096 | 6.7e-05 |
  | 0.060 | 0.3962 | 6.7e-05 |
  | 0.070 | 0.3839 | 6.6e-05 |
  | 0.080 | 0.3728 | 6.6e-05 |
  | 0.100 | 0.3537 | 6.5e-05 |
  | 0.120 | 0.3387 | 6.4e-05 |
  | 0.140 | 0.3277 | 6.3e-05 |
  | 0.160 | 0.3208 | 6.1e-05 |
  | 0.180 | 0.3186 | 5.9e-05 |
  | 0.200 | 0.3220 | 5.7e-05 |

- **Takens–Bogdanov point:** (0.2239, 0.3372) — the SN, Hopf, and homoclinic curves meet.
- **Operating corners (both at β=0.05):**
  - A=0.2: stable chimera — never-absorbers
  - A=0.5: post-homoclinic — transient (this work)
  A=0.2 sits in the stable-chimera band (below the homoclinic); A=0.5 sits above A_hc=0.4096 in the post-homoclinic region (sync globally attracting) — the regime studied in this work.

## File inventory (`paper_figures/`)

- `FIGURES_REPORT.md`
- `cp0_audit.json`
- `fig1.pdf`
- `fig1.png`
- `fig1_caption.txt`
- `fig1_trace.json`
- `fig7.pdf`
- `fig7.png`
- `fig7_caption.txt`
- `fig7_curves.json`

## Reproduce

Run from the `Anneal/` project root (not from `paper_figures/`):

```
python3 tools/paper-figures/run_all.py
```
Regenerates every artifact above from `figures.config.json`. The homoclinic trace (~4–5 min) dominates; pass `--skip-curves` to reuse the cached `fig7_curves.json` when only the plots changed.
