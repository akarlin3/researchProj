# VERIFICATION.md — Gnomon checkpoint ledger

Each checkpoint records what was *run* and what gate it had to clear. Nothing is
claimed here that was not executed.

## CP1 — scaffold + targets manifest  ✅ (this commit)

Gate: **mirrors the sibling subrepos, clean, self-consistent.**

- [x] Embedded as a top-level subrepo with its **own clean root history**, merged
      into the monorepo via `git merge --allow-unrelated-histories` (mirrors
      Lattice/Datum/Lethe). `git log -- Gnomon/` shows Gnomon's own commits.
- [x] Registered in the root `README.md` (Contents, Projects-at-a-glance, Project
      details, How-they-fit-together, Provenance).
- [x] **Reproduction-targets manifest frozen** (`gnomon/manifest.py`): 5 targets
      (T1, T3a, T3b, T3c, T4) with claimed values, frozen tolerances, and
      Fashion-**prose** provenance. `manifest.validate()` passes.
- [x] **Clean-room boundary enforced:** `_paths.py` allows `lattice` only; `caliper`
      in `FORBIDDEN`; static AST test confirms no `gnomon/` module imports Caliper.
- [x] **Clean IP:** no proprietary/clinical data, no data-like files in tree.
- [x] CP1 tests pass: `python -m pytest Gnomon/tests -q` → (recorded at commit).

## CP2 — clean-room rebuild  ✅

Gate: **runs, self-consistent.**

- [x] forward (`forward.py`), NLLS + railing (`nlls.py`), Laplace + per-voxel MCMC
      (`bayes.py`), MAF/NPE flow (`flow.py`), independent ruler (`metrics.py`),
      bootstrap (`bootstrap.py`), OSIPI loader (`osipi.py`), CP3 driver
      (`reproduce.py`) all implemented from spec — no Caliper import.
- [x] **Self-consistency gates pass** (`tests/test_cp2_selfconsistency.py`, 7 cases):
      analytic Jacobian matches finite differences; clean-signal NLLS round-trip
      recovers truth (no railing on clean data); continuity (f=0 → mono-exponential)
      holds; cohort draws ground truth from Lattice; ruler coverage exact; the
      end-to-end driver runs and emits a verdict + JSON.
- [x] Full suite green: `python -m pytest Gnomon/tests -q` → **16/16**.
- [x] Every design choice documented inline + in `docs/METHODS.md` (b-schemes, NLLS
      box/init/railing, MCMC sampler spec, MAF training spec, OSIPI ROI selection).

## CP3 — the reproduction gate  ✅  → verdict: **PARTIAL** (HARD HALT)

Ran `python -m gnomon.reproduce` (seed 20260620) → [`results/reproduction.json`](results/reproduction.json).
Full verdict + divergence report: [`VERDICT.md`](VERDICT.md).

- [x] **REPRODUCES (4/6):** T1 railing **54.2%** [52.0, 56.4] on real OSIPI (vs claimed
      54.7%, sha256-verified); T3c quantile coverage (D\* 0.90, D 0.95, f 0.96); T4 flow
      beats railed NLLS (ECE/sharpness/coverage gaps, CIs exclude 0).
- [x] **DIVERGES (2/6):** T3a Laplace-SD D\* 0.80 (claimed 0.30) and T3b MCMC-SD D\* 0.90
      (claimed 0.67) — the *severe marginal* Gaussian under-coverage.
- [x] **Cause pinned, run not asserted** (`scripts/divergence_diagnostic.py`): (1) cohort
      regime — under-coverage concentrates in the high-D\* tercile (0.63), diluted by a
      prior-spanning cohort; (2) railed-voxel uncertainty convention — the honest CRLB
      over-covers unidentified D\*, while Fashion's "overconfident" floored SD gives
      pooled 0.68 (≈ the 0.67 claim). Both are under-documented in Fashion's prose.

**Halt:** verdict is PARTIAL, so Gnomon does not auto-proceed to CP4. The divergence
changes how the numbers should be presented in the retool (re-frame marginal → conditional
coverage; document the uncertainty convention). Awaiting direction.

## CP4 — package for the retool  ⏸  (gated on the verdict)

The clean implementation + complete methods + the *reproduced* results (T1/T3c/T4) are
ready to hand to the retool; the *divergent* results (T3a/T3b) are documented with their
cause in [`VERDICT.md`](VERDICT.md). Proceed only on a decision about the re-framing.
