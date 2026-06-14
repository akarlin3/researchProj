# Gauge — manuscript (`gauge/paper/`)

The Gauge 05 deliverable: the assembled paper for the **reframed** contribution
established across Gauge 01–04.

## Contents
- `gauge.tex` — the manuscript (LaTeX, `article` class).
- `refs.bib` — bibliography.
- `figures/` — the vector-PDF figures, copied from `gauge/figures/` (regenerable
  from seed `20260613` via the Gauge 02–04 modules).
- `consistency.py` — **GATE 3** check: re-reads the committed, seeded checkpoint
  printouts (`results/*.txt`, `POSITIONING.md`, `gauge/results*.md`) and asserts
  every headline number in `gauge.tex` appears verbatim in its source. Prints a
  one-paragraph pass/fail summary; exit 0 iff all trace.
- `build.sh` — runs the consistency check, then compiles with `tectonic`.

## Build
```bash
bash gauge/paper/build.sh        # consistency check + tectonic -> gauge.pdf
# or individually:
python gauge/paper/consistency.py
( cd gauge/paper && tectonic gauge.tex )
```
Requires `tectonic` (self-contained LaTeX; fetches its package bundle and BibTeX
style on first run). Compiles to an 11-page PDF; only cosmetic overfull-hbox
warnings, no errors.

## Headline (honest, reframed — matches the gated 01–04 verdicts)
1. Model-based IVIM UQ is **broadly overconfident** (not specifically
   perfusion); conformal restores marginal coverage (|gap| ≤ 0.024). *(Gauge 02)*
2. **Conformalize-the-MDN is the sharpest valid recipe** (0.65–0.79× pure-CQR
   width at equal coverage). *(Gauge 02)*
3. The **high-D\*** compartment is an **irreducible identifiability limit**
   (CRLB(D\*)/tercile-width reaches 1.12) — the paper says *characterize*, not
   *solve*. *(Gauge 03)*
4. Robustness (weighted conformal + deployment monitor), acquisition-robustness of
   the wall, and a **qualitative** (no-coverage-claim) in-vivo demonstration.
   *(Gauge 04)*

## Out of scope (human)
The venue call (the conditional-coverage result is methodological — candidates:
*Magnetic Resonance in Medicine*, or a methods venue such as *Medical Image
Analysis* / *MELBA*) and a human read of the abstract's open-problem framing.
