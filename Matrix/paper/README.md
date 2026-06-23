# Matrix manuscript (`paper/`)

The Matrix paper: a closed-loop adaptive-RT harness on a synthetic digital twin,
grounded on real anatomy + dose geometry (Ferry). **Methods demonstration — no clinical
claim.** Target venue **PMB** (Physics in Medicine & Biology); formatted with IOP's
**`iopart` class**.

## IOP class files (vendored)

`iopart.cls`, `iopart10.clo`, `iopart12.clo`, `iopams.sty`, and `setstack.sty` are the
official IOP Publishing class files (© IOP Publishing, **freely distributable for preparing
IOP submissions**). They are **not on CTAN / in tectonic's bundle**, so they are vendored
beside `matrix.tex`; tectonic resolves them from the local directory. Do **not** load
`amsmath` with `iopart` (it clashes with iopart's `\equation*`); the manuscript uses core
math + `\boldsymbol` from `iopams` only. Author affiliation is a visible `[TODO-AUTHORS]`
placeholder until supplied at submission.

## Build

```bash
PROT=/opt/homebrew/Caskroom/miniforge/base/envs/proteus/bin/python
bash paper/build.sh        # consistency gate -> tectonic matrix.tex -> matrix.pdf
```

`build.sh` runs `consistency.py` first (it regenerates `numbers.tex` from the seeded
result JSONs and fails the build if any `\num*` macro is undefined or any load-bearing
assert breaks), then compiles with `tectonic` (fallback `pdflatex ×2`).

## Files

- `matrix.tex` — the manuscript. Every load-bearing number is a `\num*` macro from
  `numbers.tex`; assumption-dependent (Fashion/Minos) numbers are marked `\PROV`.
- `consistency.py` — the traceability gate. Reads `../results/RESULTS_CP2.json`,
  `RESULTS_CP4.json` (synthetic, regenerated offline by `../reproduce.sh`) and
  `RESULTS_FERRY_CP2.json` (grounded, regenerated from the public TCIA dataset and
  committed), writes `numbers.tex`, and asserts the headline claims (AUROC-by-construction,
  the F1 honest negative CI-excludes-0, the `loop.py` byte-identity, every CI).
- `numbers.tex` — **auto-generated; do not edit by hand.**
- `build.sh` — gate + compile.
- `matrix.pdf` — built artifact.

## Placeholders pending finalization

- `% TODO-AUTHORS` — institutional affiliation (sole author A. Karlin confirmed at GATE G;
  affiliation supplied at submission, not fabricated). Venue conversion to `iopart` is **done**.
- `% PROVISIONAL-{fashion,minos}-n` — forward-cite DOIs on publication (release gate).
- `% FUTURE-FORGE` — Forge dose engine, not built (2027); drop-in future work.
- `% TODO-CITE` / `% TODO-CITE-VERIFY` — foundational references for the author to supply/verify.
