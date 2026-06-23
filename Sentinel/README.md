# Sentinel

A pre-registered separation gate: does **regret-targeted decision-stopping** halt at a
different time than **coverage-targeted stopping** (ACI/conformal-PID; a WATCH-style
coverage-changepoint alarm) over a fractionated RT course?

**Verdict: 🔴 RED.** On the mandated Matrix substrate, with both rules made sequential and
calibrated to the same anytime false-alarm budget, the regret-stop does **not** robustly
halt earlier than WATCH. The separation was the paper; it is not there.

What the run leaves behind (all green, all honest):
- **`sentinel-core/`** — the fractionated-session **enabler** (imports the Matrix twin
  read-only, byte-identity enforced), faithful **ACI/PID + WATCH** baselines, and a
  **separation harness** with voxel-bootstrap CI and a pre-registered refute.
- **`sentinel-core/RESULTS.md`** — the verdict evidence (multi-patient + sweep).
- **`sentinel-core/POSITIONING.md`** — novelty record + the kill mechanism, including a
  disclosed earlier false positive (a per-session-vs-sequential confound) that was fixed,
  not rescued.

Reproduce: `cd sentinel-core && bash reproduce.sh`.

Clean-IP, synthetic twin only, no clinical claim. Matrix (PR #64) is HELD and remains
byte-unchanged — proven at import time and in `tests/test_enabler.py`.
