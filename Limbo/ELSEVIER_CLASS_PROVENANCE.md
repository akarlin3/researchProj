# Vendored Elsevier class files — provenance

`limbo_phiro.tex` targets **Physics and Imaging in Radiation Oncology** (phiRO; Elsevier, open
access) as a **review article**, and is typeset with Elsevier's official `elsarticle` class and the
numbered (Vancouver) bibliography style phiRO specifies. The files are vendored here so the
manuscript builds self-contained and the exact class/style version is **pinned** for
reproducibility — mirroring how the IOP class was vendored for the prior PMB target (see
`IOP_CLASS_PROVENANCE.md`).

| file | role | license |
|---|---|---|
| `elsarticle.cls` | Elsevier journal article class | LaTeX Project Public License (LPPL) 1.3 or later — redistribution/modification permitted |
| `elsarticle-num.bst` | Elsevier numbered (Vancouver) bibliography style — references numbered in order of appearance, in square brackets, as phiRO's Guide for Authors requires | LaTeX Project Public License (LPPL) |

- **Upstream source:** CTAN, `https://ctan.org/pkg/elsarticle` (package "elsarticle", maintainer
  C. V. Radhakrishnan, Elsevier Ltd). The current CTAN release at audit time is **v3.5
  (2026-01-09)**.
- **Vendored version (what actually builds this manuscript):** the files here are the versions
  shipped by **tectonic's default bundle** (a frozen TeX Live snapshot), materialized into the
  build by compiling once with `tectonic` and copying the resolved files out:
  - `elsarticle.cls` — **v3.3, dated 2020/11/20** (`\ProvidesClass` header; Copyright 2007–2020
    Elsevier Ltd).
  - `elsarticle-num.bst` — **Version 2.1**, `$Id: elsarticle-num.bst 194 2020-11-23 ...$`.
  The numbered-Vancouver behaviour (front matter, `[number]` citation scheme, order-of-appearance
  reference list) is identical between v3.3 and the current v3.5; pinning v3.3 keeps the build
  byte-reproducible against the toolchain present in this environment.
- **Retrieved:** 2026-06-23.
- **Integrity (sha256):**
  - `elsarticle.cls`: `23a341acd0e35c20148c836fc6e9251549aa45b059fc1750356ca2f3f11f3cb0`
  - `elsarticle-num.bst`: `aaf9bce5453a0c191aa5ff4d200dc3a9b4170a6067ca723ec4f1edcd7ca14edf`
- **Why vendored:** unlike the IOP `iopjournal`/`iopart` class (which is **not** on CTAN and
  **cannot** be fetched by tectonic, forcing the vendoring of `iopjournal.cls` + `orcid.pdf`),
  `elsarticle` **is** on CTAN and tectonic fetches it automatically. Vendoring is therefore not
  strictly required for `elsarticle`; it is done here only to (a) pin the exact class + bibliography
  version so the rendered manuscript and reference list are reproducible, and (b) mirror the
  self-contained build contract of the IOP version. The class header carries the LPPL, which
  permits redistribution.

The class option `[review,number]` selects single-column review formatting (1.5× line spacing for
reviewers) with natbib's numbered citation scheme; `\biboptions{sort&compress}` adds ascending,
range-compressed multi-citations to match the prior rendering. `\bibliographystyle{elsarticle-num}`
selects the vendored numbered-Vancouver style. phiRO sets the real running head / journal branding
at production; `\journal{Physics and Imaging in Radiation Oncology}` only labels the author-prepared
preprint and does not affect content.
