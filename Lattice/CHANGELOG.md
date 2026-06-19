# Changelog

All notable changes to Lattice are documented here.

## [0.1.0] — 2026-06-19

Initial digital reference object (DRO).

### Added
- Five clean-room IVIM forward generators (`lattice.generators`): bi-exponential
  plus four misspecification families (gamma & log-normal velocity dispersion,
  stretched-exponential, tri-exponential), each reducing to bi-exp at a
  continuity limit (three exact, one asymptotic).
- Cohort schema, physiological priors (`uniform`, `realistic` Gaussian copula),
  and seeded generation (`lattice.cohort`); family-invariant base draws.
- Standardized calibration-evaluation interface (`lattice.evaluate`): an
  estimator contract and a `to_scorer_inputs` adapter for any quantile scorer
  (canonically Caliper). Core imports nothing back (strict one-way dependency).
- NLLS round-trip self-validation (`lattice.selfcheck`, scipy optional).
- Optional OSIPI download-on-demand + provenance (`lattice.osipi`,
  `scripts/fetch_osipi.py`); synthetic-only in tree.
- Publication/citation gate (`lattice.publication`), OFF by default.
- Tests (continuity, round-trip, schema, interface, publication), worked
  examples, and docs (`DRO_SPEC`, `POSITIONING`, `CLEANROOM`).

### Verified
- Continuity residuals `0.000e+00` (exact families), `3.9e-09` (gamma, k=1e8).
- Clean bi-exp round-trip max relative error `7.99e-10` (gate `< 1e-2`).
