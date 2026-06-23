# ASSUMPTIONS — projLevy

Pinned dependencies, regime scoping, and the clean-IP gate. Mirrors `Minos/future/ASSUMPTIONS.md`.

## 1. Environment / pinned versions
- Python: the `proteus` conda env (`/opt/homebrew/Caskroom/miniforge/base/envs/proteus/bin/python`).
- Verified: numpy 2.3.5, scipy 1.17.1 (`reproduce.sh` sets `PROT`; override with `PROT=...`).
- Core deps: numpy ≥ 1.23, scipy ≥ 1.9, matplotlib ≥ 3.6 (`levy-core/pyproject.toml`).

## 2. Regime scoping (every claim is scoped to these)
- **Forward model:** stretched-exponential `S(b; S₀, D, α) = S₀·exp(−(bD)^α)` (lead lane).
  Joint CTRW / fractional Bloch–Torrey (α, β) is Phase 3 only (built only if Gate C stands).
- **Noise:** Rician magnitude, `σ = S₀/SNR` (b=0 magnitude-SNR convention).
- **Truth ranges:** D = 1.5×10⁻³ mm²/s (tissue), α ∈ [0.6, 1.0] (stretched-exp tissue range).
- **Acquisition:** b ∈ {0 … b_max}, b_max ∈ {1000, 2000, 3000} s/mm²; n_b ∈ {4 … 16}.
  Clinical anomalous-diffusion DWI = **few b-values** (n_b 4–6); n_b ≥ 8–12 = dedicated research.
- **Realistic SNR band:** [20, 60] (clinical DWI at b=0; research up to ~100).
- **Pre-registered wall threshold:** `cv_α = √(CRLB_α)/α = 0.20` (fixed before results; not retuned).

## 3. DATA SOURCE — clean (synthetic only); the IP gate

**Resolved at CP0; confirmed clean. No `pancData3`, no MSK, no clinical/medical data is touched.**

| source | kind | license / provenance | used by |
|---|---|---|---|
| Levy forward model (`levy/forward.py`) | **synthetic** | seeded, reproducible, closed-form | CP0 (all) |
| Rician noise (`levy/noise.py`) | **synthetic** | seeded explicit Generators | CP0 (bootstrap, profile) |
| Ouroboros GL operators + A(α) law | **reused read-only** | Ouroboros (AGPLv3, fully synthetic) | `levy/glreuse.py` cross-check only |

**Default: fully synthetic. No in-vivo or open imaging data is required or used for CP0.** The
clean-IP gate passes: nothing private, nothing medical. (If a future stage validates against an
open DWI cohort, it would go behind an explicit, documented fetch — not required for the
identifiability result, which is forward-model theory + simulation.)

## 4. Re-validation contract
`reproduce.sh` re-runs CP0 (and later stages once built). FAST (default) is a smoke run with
modest bootstrap; `FULL=1` runs full-N bootstrap CIs. The pre-registered REFUTE is wired into
the gate: if a future change makes α recoverable across the realistic clinical band (no wall),
`verify_cp0.py` reports `refuted=True` honestly rather than failing.
