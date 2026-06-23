# Pre-coding verification — Levy reuse surface, net-new layer, and clean-IP gate

The directive (GATE A) required auditing the Ouroboros reuse surface and the Minos house
template, and confirming Ouroboros has **no** CRLB/Fisher layer, before any code. Two
independent passes were run (a structured Explore audit and a direct grep). Verdict: **gate
PASSED**.

## The three required checks

### 1. What is reusable from Ouroboros, as-is?
**CONFIRMED reusable (read-only):** the Grünwald–Letnikov weight generator `gl_weights`
(`Ouroboros/ouroboros_fractional_sindy.py:7`), the time/space fractional-derivative operators
`gl_derivative_time/space` (`:14`, `:27`), the Gaussian `add_noise` routine
(`Ouroboros/ouroboros_identifiability.py:56`), and the noise-amplification law
`A(α) = dt^{-2α}·‖w(α)‖²` (`Ouroboros/ouroboros_noise_analysis.py:20-45`,
`compute_analytic_factors`). → Implemented: `levy/glreuse.py` carries a clean copy of
`gl_weights` and `noise_amplification` with provenance headers (read-only reuse). Levy adopts
the Minos explicit-Generator seeding discipline rather than Ouroboros' bare `np.random`, and
builds **Rician** (not Gaussian) noise, since the MRI forward model requires it.

### 2. Does Ouroboros contain any CRLB / Fisher / identifiability layer?
**CONFIRMED ABSENT.** A grep of all Ouroboros `*.py` for
`fisher|crlb|cram[ée]r|fim|information.matrix|profile.likelihood` returns **zero** code hits;
the only `bootstrap` matches are **Ensemble-SINDy** resampling over time-indices for model
*discovery* (`ouroboros_cp3_fairgrid.py`), not statistical CIs on a Fisher/CRLB estimate. The
word "identifiability" in Ouroboros is purely **empirical** (clean-limit SINDy R² specificity +
SNR sweeps), never Fisher-based. → Therefore the Fisher/CRLB/profile-likelihood layer
(`levy/fisher.py`, `levy/identifiability.py`, `levy/wall.py`) is **net-new and is Levy's
contribution.** This check is re-run live inside the CP0 gate (`verify_cp0.py` check 1/5).

### 3. Does any MRI signal-decay forward model already exist to reuse?
**CONFIRMED ABSENT.** Ouroboros has only (i) an integer-order 3-field reaction–diffusion PDE
and (ii) a linear *fractional-derivative ODE* test system. There is **no** b-indexed signal
attenuation `S(b; …)`, **no** Rician noise, and **no** stretched-exponential / CTRW /
Bloch–Torrey model (grep for `bloch.?torrey|stretched.?exp|ctrw|signal.?decay|b.?value|
attenuation` → zero hits). → The forward model `levy/forward.py` is built net-new.

## Findings that shaped the design
- **The wall is a 2-D object over (SNR, b-design), not a single SNR threshold.** The number of
  b-values is the dominant driver: α walls out within the realistic SNR band for sparse clinical
  acquisition (n_b ≲ 6) and recedes for dense research acquisition (n_b ≥ 8). The design maps the
  full surface (`levy/wall.wall_surface`) rather than reporting one cell.
- **Analytic CRLB wall ≠ empirical wall.** The CRLB is an information *lower bound* on the wall;
  the finite-sample MLE walls out at slightly **higher** SNR. Both are reported, separately
  labelled, so the CI is never implied to bracket the wrong point.
- **Extending b_max with few b-values trades against α–D collinearity** (ρ_αD → −0.87 at
  n_b=4, b_max=3000): a structural finding, reported as part of the deliverable.

## Clean-IP gate
**PASSED.** Fully synthetic; no medical / real data. See `ASSUMPTIONS.md` §3. Ouroboros is
reused read-only and is itself fully synthetic (PDE-generated, AGPLv3); no `pancData3`, no
pancreas, no TCIA, no clinical data is touched.

## Sources / provenance
- Ouroboros reuse surface: `Ouroboros/ouroboros_fractional_sindy.py`,
  `ouroboros_identifiability.py`, `ouroboros_noise_analysis.py` (read-only).
- House template: `Minos/` (pyproject flat layout, `pytest.ini`, `reproduce.sh`, `verify_cp1.py`,
  `_paths.py`, `VERIFICATION.md`, `POSITIONING.md`).
- Prior-art neighbours to distinguish: see `levy-core/POSITIONING.md` (Coeurjolly–Istas 2001;
  Spilling & Barrick 2022, PMID 36054778; Poot 2010 / Chuhutin 2017 — scoped out).
