# Positioning: why Lattice stands alone

Lattice was vetted against the four nearest artifacts before being built, with an
explicit fallback to "fold into Caliper as a data module" if it proved too thin.
**Verdict: it stands alone**, conditional on the alternative-model generator
families being the shipped core (they are).

## Distinct contribution

A **UQ-calibration-focused reference object**: physiologically-grounded reference
*parameter distributions* + *alternative-model* generators (gamma/lognormal
velocity dispersion, tri-exponential, stretched-exponential, each with an exact
or asymptotic bi-exp continuity limit) + a *standardized calibration-evaluation
API*. That combination is not packaged anywhere else in the family.

## Versus the neighbours

- **OSIPI TF2.4** (upstream of `Fashion/`) is **accuracy-focused** — a fitting-
  method benchmark and a reference DRO for point-estimate correctness. Lattice
  targets *calibration* (coverage/sharpness), including under **model
  misspecification**, which an accuracy DRO does not provide.
- **Gauge's embedded cohort** is **internal-only** — seeded for the
  reproducibility of one paper's results, never released as a standalone reusable
  resource. Lattice is a packaged, versioned, documented DRO.
- **Caliper** is the **scorer** (the ruler: coverage/ECE/sharpness). Lattice is
  the **data** the scorer consumes. Complementary, not overlapping.
- **Datum** does not exist in the monorepo (no directory, no code, no planning
  note) — no collision.

## Why not fold into Caliper?

Caliper already ships a *bi-exponential* `synthetic_cohort()`. If Lattice were
only that, it would belong inside Caliper. It is not: its core is the
**alternative-model families** (absent from Caliper) and its role is a **shared
reference object** consumed by *many* downstream users — Caliper, the
Fashion/Gauge/Minos reproductions, and planned siblings (Vernier, Echo). A
reference object consumed by many, with its own versioning/DOI cadence, is
conventionally a separate citable resource (exactly as OSIPI's DRO is separate
from its fitting code). Burying it in one consumer would couple the resource to
the scorer and blur Caliper's clean "ruler" identity.

## Family map (IVIM program)

- **Fashion** — which UQ paradigms actually cover D\*.
- **Gauge** — distribution-free conformal coverage + the high-D\* identifiability wall.
- **Caliper** — the reusable calibration *scorer* (ruler).
- **Minos** — the *decision* value of a calibrated error bar + a label-free monitor.
- **Lattice** — the *reference object* (this project): the ground-truth data and
  generators all the above (and planned Vernier/Echo) can benchmark against.
