# Minos — the decision value of a calibrated error bar

*Paper:* **"Minos: the decision value of a calibrated uncertainty — A
decision–calibration gap and a label-free validity floor for quantitative MRI"**
(theory complete; applied half provisional; target *MRM*).

Minos prices *the error bar itself* — not the point estimate, not a population
parameter — on a treat/spare/escalate clinical decision. It defines the **Value of
Calibration** (utility lost when the error bar is mis-scaled) and the **Value of the
Trust-Gate** (utility recovered by detecting when uncertainty is untrustworthy under
shift). Its v2 result is a **decision–calibration gap** `G = tau* - tau_stat`: the
scale that achieves nominal coverage and the scale that maximizes expected utility
*diverge* under skew and cost asymmetry. Its v3 result is a label-free
deployment-validity monitor, honest about what it can detect (observable-driven
shift, AUC ≫ 0.5) versus cannot (hidden truth-shift, AUC ≈ 0.5). The theory core is
100% synthetic, deterministic, and gate-checked; the applied half is speculative by
construction (it assumes Fashion and Gauge survive to publication as submitted).

## Layout

- [`minos-core/`](minos-core/) — the validated, data-independent theory core: `VoC`/`VoTG`,
  the v2 decision–calibration gap, the v3 label-free validity monitor
  (`utility.py`, `decision.py`, `voi.py`, `calibration.py`, `correction.py`,
  `monitor.py`), 33 gate-as-assertion tests. See
  [`minos-core/README.md`](minos-core/README.md) and
  [`minos-core/POSITIONING.md`](minos-core/POSITIONING.md).
- [`theory/`](theory/) — **Plumbline**: analytic hardening of the minos-core findings
  — the gap-scaling law (Theorem 1), the detectability bound and hidden-channel
  impossibility proof (Theorem 2), the second-order value-of-information law
  (**Delphi**, Proposition 3) — all machine-checked against the v2/v3 simulation.
  See [`theory/README.md`](theory/README.md).
- [`future/`](future/) — the complete, applied Minos paper, built now under the
  assumption that **Fashion** and **Gauge** publish as submitted: wires the
  validated theory (read-only) to a decision/monitor layer that consumes
  Fashion's calibrated IVIM posteriors, contextualized by Gauge's coverage and
  high-D\* identifiability wall. CP1–CP4 all done, manuscript compiles
  (`future/paper/minos.tex`), every number traces to a seeded result — but
  **speculative by construction** and flagged PROVISIONAL throughout. See
  [`future/README.md`](future/README.md) and
  [`future/ASSUMPTIONS.md`](future/ASSUMPTIONS.md) for the SOLID-vs-PROVISIONAL
  split.
- [`paper/`](paper/) — standalone theory-only manuscript sections
  (`theory_standalone.tex`) independent of the Fashion/Gauge assumption.
- [`sibyl/`](sibyl/) — **Sibyl**: a related but distinct subproject reusing this
  house's tooling — out-of-distribution detection and trustworthiness flagging
  for quantitative diffusion-MRI microstructure estimation (IVIM via µGUIDE),
  validated against public ACRIN-6698 repeat-acquisition data (Tier 1 synthetic,
  Tier 2 in-vivo breast DWI done; Tier 3 glioma/liver is a stub). See
  [`sibyl/README.md`](sibyl/README.md).

## Reproduce

```bash
cd minos-core && python -m pytest              # theory-core gate suite (33 tests)
bash future/reproduce.sh                       # one-command CP1-CP4 re-validation (applied half)
```

`theory/`'s per-theorem gate scripts (`gap_scaling.py`, `detectability.py`,
`impossibility_check.py`, `confirm.py`, `voi_value.py`) are run individually — see
[`theory/README.md`](theory/README.md) for the exact commands.
