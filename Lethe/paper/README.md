# Lethe manuscript — the Echo portion

`lethe.tex` (`ebgaramond` + `microtype`) is the Echo portion of the Lethe paper: the
honest-limitation / constrained-validation write-up of the test–retest *scale* check, built
on the CP3 → **LETHE** verdict. Its traceability gate is `consistency.py` →
auto-generated `numbers.tex`.

Build:
```bash
bash build.sh        # tectonic preferred; falls back to pdflatex x2
```

Discipline (mirrors `Minos/future/paper/`):
- every `\num*` macro in `lethe.tex` is defined in `numbers.tex`, which `consistency.py`
  regenerates from the seeded `results/RESULTS_VALIDATION.json` (PROVISIONAL real-data
  verdict) and `results/RESULTS_HARNESS.json` (SOLID method self-test);
- every Fashion/Gauge/Minos-dependent number carries the PROVISIONAL marker (`\PROV`);
- `consistency.py` is the gate: zero undefined macros + load-bearing internal-consistency
  asserts (verdict==LETHE, coverage(D) under-scaled, R<1, self-test all-pass, distinctness).
