# Citing Caliper

Caliper is open-source tooling. If you use it, please cite the **software**. The
toolkit also ships *synthetic, qualitative* reproductions of two associated IVIM
manuscripts — but **both are pre-publication**, so this page renders them as
*submitted / in review*, never as *published*.

> **Single source of truth.** Everything below is generated from
> [`caliper.publication`](../caliper/publication.py). A paper counts as
> *published* **iff** it has a real `paper_doi` there; until then
> `caliper.publication.publication_enabled()` is `False` and the feature is OFF.
> The Zenodo DOIs shown are **software code-archive** DOIs for the sibling
> repositories — they are **not** publication DOIs and do not flip the gate.

```python
from caliper.publication import PUBLICATION, publication_enabled, bibtex
publication_enabled()            # False in the default (shipped) state
PUBLICATION["gauge"].published   # False -- no real paper DOI yet
print(bibtex("gauge"))           # @unpublished while pre-publication
```

## Software

```bibtex
@software{Caliper,
  author  = {Karlin, Avery},
  title   = {Caliper: an IVIM uncertainty-quantification calibration toolkit},
  version = {0.1.0},
  year    = {2026},
  url     = {https://github.com/akarlin3/ResearchProj/tree/main/Caliper},
  license = {MIT}
}
```

## Associated manuscripts (pre-publication)

These are the papers Caliper's reproduction modules map to. **Neither is
published**; there is no publication DOI for either. Cite them as below
(`@unpublished`), and **do not** present them as published or accepted. The
placeholder `10.XXXX/XXXXX` is to be replaced **only** once a paper actually
publishes.

### Gauge — conformal coverage & the D\* identifiability wall

- **Status:** submitted to *Magnetic Resonance in Medicine* (2026); pre-publication.
- **Software archive (NOT a publication DOI):** `doi:10.5281/zenodo.20686273`.
- **Reproduced (synthetic, qualitative):** `caliper.repro_gauge` —
  `python examples/gauge_repro.py`. This is *not* a published or independently
  validated result.

```bibtex
@unpublished{Karlin_Gauge_IVIM_Conformal,
  author = {Karlin, Avery},
  title  = {Distribution-Free Conformal Coverage for IVIM Parameter Maps, and the Identifiability Wall in the Pseudo-Diffusion Compartment},
  year   = {2026},
  note   = {Submitted to Magnetic Resonance in Medicine (2026); pre-publication -- no publication DOI yet (placeholder 10.XXXX/XXXXX). Software archive: doi:10.5281/zenodo.20686273.}
}
```

### Fashion — calibration & efficiency of IVIM uncertainty estimates

- **Status:** in review at *Magnetic Resonance in Medicine* (2026, Manuscript ID
  MRM-26-27109); pre-publication.
- **Software archive (NOT a publication DOI):** `doi:10.5281/zenodo.20649669`.
- **Reproduced (synthetic, qualitative):** `caliper.baselines` —
  `python examples/fashion_repro.py`; see
  [`fashion_reproduction.md`](fashion_reproduction.md). This is *not* a published
  or independently validated result.

```bibtex
@unpublished{Karlin_Fashion_IVIM_Calibration,
  author = {Karlin, Avery},
  title  = {Calibration and Efficiency of Uncertainty Estimates in Intravoxel Incoherent Motion Imaging: Quantile Intervals, Cross-Paradigm Comparison, and a Cramer-Rao Audit of Amortized Posteriors},
  year   = {2026},
  note   = {In review at Magnetic Resonance in Medicine (2026), Manuscript ID MRM-26-27109; pre-publication -- no publication DOI yet (placeholder 10.XXXX/XXXXX). Software archive: doi:10.5281/zenodo.20649669.}
}
```

## How the feature activates

The reproduction & citation feature is **dormant by default** and turns on one
paper at a time, only when that paper actually publishes:

1. Put the real publication DOI in the paper's `paper_doi` field in
   `caliper.publication.PUBLICATION` (replacing `None`).
2. `PaperRef.published` and `publication_enabled()` flip to `True` automatically;
   the bibtex entry becomes `@article` with the DOI, and the reproduction's
   provenance note flips to *"validated against the published result."*
3. Update this page's status line and `CITATION.cff` (`status:` and add the real
   `doi:`) to match.

Until then, everything here is honest by construction: **pre-publication, not
published.**
