# sentinel-core

**Regret-targeted decision-stopping vs coverage-targeted stopping — a pre-registered
separation gate that returned 🔴 RED.**

projSentinel asked one question: does a *decision-value* stopping rule (Minos's
regret-targeted monitor, accumulated over a fractionated RT course) halt at a meaningfully
different time than the coverage-targeted baselines — ACI/conformal-PID (recalibrate
forever) and a WATCH-style conformal-martingale alarm? **On the mandated Matrix substrate,
no.** When the rules are placed on equal footing the separation is not there. See
`RESULTS.md` for the verdict and `POSITIONING.md` for the novelty record and the kill
mechanism.

## Design discipline
- **Built on the real baselines**, implemented from scratch and sanity-checked against
  their published behavior: ACI (arXiv:2106.00170), conformal PID (arXiv:2307.16895),
  WATCH (arXiv:2505.04608).
- **Matrix is HELD and byte-unchanged.** projSentinel imports the Matrix twin
  **read-only**; `matrix_bridge` asserts `matrix/loop.py`'s sha256 against the Gate-A
  anchor at import time and refuses to run against a modified Matrix.
- **Fair fight.** The regret-stop and WATCH are both *sequential* and calibrated to the
  *same* anytime false-alarm budget `δ`, so only the targeted statistic differs.
- Clean-IP, synthetic twin only, **no clinical claim**.

## Layout
```
sentinel-core/
  sentinel/
    config.py        frozen SentinelConfig (course, decision, drift, gate thresholds)
    seeding.py       GLOBAL_SEED, make_rng
    matrix_bridge.py read-only Matrix import, loop.py byte-identity anchor
    course.py        the fractionated-session axis ENABLER (accumulating drift wrapper)
    monitor.py       the regret-targeted monitor M (ported from Minos)
    baselines.py     ACI/conformal-PID + WATCH conformal test martingale
    stopping.py      the sequential regret-targeted CUSUM stop (anytime-δ calibrated)
    separation.py    CP0 separation test + voxel-bootstrap CI + pre-registered refute
  experiments/run_cp0.py   the verdict experiment (instrument control + Matrix sweep)
  tests/                   gate tests (enabler / baselines / separation)
```

## Install & run
```bash
# proteus env (numpy/scipy). Matrix path defaults to the matrix-subrepo worktree.
export SENTINEL_MATRIX_PATH=/path/to/Matrix
bash reproduce.sh                 # anchor check -> pytest -> CP0 verdict
# or:
PYTHONPATH=. python experiments/run_cp0.py
PYTHONPATH=. python -m pytest -q
```

## Status
🔴 **RED — no robust separation; no paper.** The enabler, the faithful baselines, and the
separation harness are kept as a reproducible negative result and a reusable instrument.
