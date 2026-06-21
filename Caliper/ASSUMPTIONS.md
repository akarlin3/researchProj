# Caliper — assumptions manifest

> Caliper is the **ruler-as-code**: the numpy-only calibration ruler
> (`caliper.metrics`) plus the synthetic baselines and conformal layer that the
> downstream repos (Datum, Vernier, Minos, Lethe) score against. This manifest
> records what Caliper assumes about the **retooled Fashion** ruler it
> reproduces, and what is SOLID vs PROVISIONAL.

## Status split

| Component | Depends on Fashion (in review)? | Status |
|---|---|---|
| Calibration ruler `caliper.metrics` (coverage / ECE / sharpness / pinball / interval-score / tercile-conditional coverage) | No — numpy-only, definitions from the literature | **SOLID** |
| Conformal layer `caliper.conformal`, benchmark `caliper.benchmark` | No — synthetic, self-contained | **SOLID** |
| NLLS baseline `caliper.baselines` (railing rate; honest-CRLB SD convention) | No — synthetic phantoms only | **SOLID** |
| Fashion reproduction framing + publication gate `caliper.publication` | **Yes** — pins the retooled Fashion manuscript | **PROVISIONAL** |

## The retooled-Fashion pin (PROVISIONAL)

- **Manuscript:** *"Boundary-railing of conventional NLLS fits as an
  assumption-free pseudo-diffusion identifiability diagnostic in IVIM MRI …"*
  (retooled, **boundary-railing-first**).
- **Venue:** **NMR in Biomedicine** — *in review* (moved from MRM; submission ID
  pending). Machine-readable pin: `caliper.publication.PUBLICATION["fashion"]`
  (`venue="NMR in Biomedicine"`, `status="in_review"`, `paper_doi=None`).
- **Ruler scope:** the calibration ruler is a **scoped, ground-truth-only
  secondary** — it scores intervals against truth and lives on synthetic / DRO
  data; it cannot be applied to a real scan. The **assumption-free primary** is
  boundary railing on real data, which needs no ruler.
- **SD convention:** **honest CRLB is the default and the recommended
  convention** (`caliper.baselines.NLLSIVIMEstimator.sd_convention="honest"` —
  wide where D\* is unidentified). The `"floored"` convention (narrow railed-D\*
  SD, `railed_sd_floor=3.0` ≈ Gnomon's 0.003 mm²/s) is shipped only as a labelled
  illustration of how the now-dropped marginal 0.30 / 0.67 D\* severity arose; it
  is never the default. See Gnomon `docs/METHODS.md` §5b.
- **Load-bearing readout:** *conditional* coverage by D\* tercile (the high-D\*
  identifiability wall), not a marginal severity.

## The assumption

Caliper assumes the retooled Fashion's findings survive review: (a) boundary
railing is a real, assumption-free identifiability signature on real data;
(b) under the honest CRLB, symmetric Gaussian D\* intervals under-cover
*conditionally* in the high-D\* tercile; (c) the skew-aware quantile interval
restores marginal coverage with a residual high-D\* gap; and (d) an amortized
(flow) posterior is better-calibrated and sharper than the railed NLLS. Caliper
reproduces (a)–(d) *qualitatively on synthetic data only*; the real-data /
clinical numbers live in the paper.

## Clean IP

Synthetic-only by construction: every cohort is generated in-repo via
`caliper.forward.synthetic_cohort` (fixed seeds). No clinical / in-vivo data, no
`pancData3` / MSK, in tree or history. Caliper does **not** import any sibling's
ruler; downstream repos import Caliper's ruler read-only.

## Re-validation (one command)

```bash
pytest -q                              # ruler + baseline + publication-gate tests
python examples/fashion_repro.py       # synthetic railing + scoped-ruler scorecard
```

When NMR in Biomedicine clears, fill the real `paper_doi` in
`caliper.publication.PUBLICATION["fashion"]`; that single edit flips the gate,
the bibtex (`@unpublished` → `@article`), and the provenance language. Until then
the Fashion-dependent framing stays PROVISIONAL by construction.
