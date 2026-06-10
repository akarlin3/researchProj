"""Verification harness for a finite-N kinetic-theory correction to the
collective flow (anneal v6, exploratory branch).

Integrates the deterministic reduced flow (tools/reduced-ode/reduced_core.rhs_3d)
PLUS a pluggable correction

    f_corr(rho1, rho2, psi, N) -> (drift_vec (3,), diff_matrix B (3,3))

via Euler-Maruyama:  X' = X + dt*(rhs_3d(X) + drift) + sqrt(dt) * (B @ xi),
with rho1, rho2 kept in [0,1] by reflection (same operator as Appendix B,
tools/noise-test/em_core._reflect) and psi unbounded. The stepping, sampling
cadence (sampleStride), early stop and capture/breath detection replicate
em_core.em_run EXACTLY, so with B = sigma*I and zero drift the harness is
bit-identical to the Appendix B driver (asserted in null_test.py), and with
f_corr == 0 it reproduces the deterministic DOP853 capture to < 0.3%
(zero_corr_gate, asserted before any scoring run).

Ensemble: the section-6.4 seed-mapped collective ICs (reduced_results/
reduced_runs.jsonl, sorted by seed, first 200 per N), N in {8, 16, 32, 64},
200 realizations per cell, EM dt = 0.01, t_max = 2000 s, seed base 9_000_000
(disjoint from Appendix B's 7M and the mech probe's 8M).

Deterministic per-N references are the COMMITTED per-IC DOP853 capture times
from reduced_runs.jsonl; anchor_gate() asserts their medians match the
committed Appendix B anchors exactly and spot-recomputes rows with DOP853.

No theory lives in this file: the correction term is whatever the human
supplies in f_corr.py at CP2. Run: see run_probe.py (scoring) and
null_test.py (self-test that the harness rejects a known-wrong mechanism).
"""
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "tools/reduced-ode"))
sys.path.insert(0, str(ROOT / "tools/beta010-slice"))
sys.path.insert(0, str(ROOT / "tools/noise-test"))
import reduced_core as rc  # noqa: E402
from em_core import _detect, _reflect  # noqa: E402
from _reuse import load_funcs  # noqa: E402

CFG = rc.load_config()
rc.set_config(CFG)
A = 0.5
BETA = 0.05
NS = [8, 16, 32, 64]
N_REAL = 200
DT = 0.01
T_MAX = 2000.0
SEED_BASE = 9_000_000
OUT = ROOT / "kinetic_results"

# canonical PR#41 breath machinery (verbatim functions via the AST loader),
# wired exactly as tools/noise-test/mech_probe.py does
ABS = load_funcs(ROOT / "tools/absorption-recampaign/analysis.py")
BR = {"smoothWindowSec": 2.0, "minProminenceFrac": 0.1, "minPeakSepFrac": 0.5,
      "minCyclesForPhase": 2, "autoPeriodFloorSec": 8.0}
ABS["BR"] = BR
MIN_CYCLES = BR["minCyclesForPhase"]
detect_peaks = ABS["detect_peaks"]

# section-6.4 seed-mapped ICs + committed per-IC deterministic (DOP853) capture
_rr = [json.loads(s) for s in (ROOT / "reduced_results/reduced_runs.jsonl")
       .read_text().splitlines() if s.strip()]
ICS, DET = {}, {}
for _N in NS:
    _rows = sorted([r for r in _rr if r["N"] == _N], key=lambda r: r["seed"])[:N_REAL]
    assert len(_rows) == N_REAL, f"expected {N_REAL} committed rows at N={_N}"
    ICS[_N] = [(r["Rsync0"], r["Rincoh0"], r["dphi0"]) for r in _rows]
    DET[_N] = np.array([r["t_capture"] for r in _rows], float)

DET_MED = {N: float(np.median(DET[N])) for N in NS}
DET_POOLED_MED = float(np.median(np.concatenate([DET[N] for N in NS])))


def anchor_gate():
    """Committed deterministic references must match the Appendix B anchors
    exactly, and a 3-row-per-N DOP853 recompute must reproduce the committed
    per-IC values. Returns the gate dict; raises AssertionError on failure."""
    anchors = json.load(open(ROOT / "noise_results/noise_results.json"))["deterministic_ref"]
    med_ok = {str(N): abs(DET_MED[N] - anchors[str(N)]) < 1e-9 for N in NS}
    assert all(med_ok.values()), f"deterministic medians off anchors: {DET_MED} vs {anchors}"
    rows_checked, worst = [], 0.0
    for N in NS:
        for j in range(3):
            res = rc.reduced_run_3d(list(ICS[N][j]), rc.Params(A=A, beta=BETA), T_MAX, CFG)
            dev = abs((res["t_capture"] if res["captured"] else np.nan) - DET[N][j])
            worst = max(worst, dev)
            rows_checked.append({"N": N, "j": j, "committed": DET[N][j],
                                 "recomputed": res["t_capture"], "dev": dev})
    assert worst < 1e-9, f"DOP853 row recompute deviates: {worst}"
    return {"medians_match_appendixB_anchors": med_ok,
            "row_recompute_worst_dev_s": worst, "rows_checked": rows_checked}


def em_run_corr(state0, fc, N, seed, dt=DT, t_max=T_MAX):
    """One Euler-Maruyama realization of rhs_3d + f_corr. Mirrors
    em_core.em_run step-for-step (same rng usage: one standard_normal(3) per
    step when the diffusion is non-zero; same reflection, sampling, early
    stop, detection), generalized to an arbitrary (drift, B) correction.
    Returns the em_core detection dict plus the canonical PR#41 breath
    period / capture phase (Tb_canon, capture_phase_canon)."""
    p = rc.Params(A=A, beta=BETA)
    b = CFG["boundary"]
    theta, recThr, recWin = b["theta"], b["recoveryThreshold"], b["recoveryWindowSec"]
    dt_sample = CFG["timescale"]["sampleStride"]
    sample_every = max(1, int(round(dt_sample / dt)))
    sqdt = np.sqrt(dt)
    rng = np.random.default_rng(seed)
    n_steps = int(round(t_max / dt))

    x = np.array(state0, float)
    m_samples = [min(x[0], x[1])]
    m_lo = m_hi = m_samples[0]
    stop_level = 0.97
    for step in range(1, n_steps + 1):
        f = rc.rhs_3d(x, p)
        drift, B = fc(x[0], x[1], x[2], N)
        B = np.asarray(B)
        if B.any():
            xi = rng.standard_normal(3)
            # scale B first so B = sigma*I reproduces em_core's
            # (sigma*sqdt)*xi bit-for-bit (see null_test.py equivalence gate)
            x = x + dt * (np.asarray(f) + drift) + (B * sqdt) @ xi
        else:
            x = x + dt * (np.asarray(f) + drift)
        x[0] = _reflect(x[0])
        x[1] = _reflect(x[1])
        if step % sample_every == 0:
            mm = min(x[0], x[1])
            m_samples.append(mm)
            if mm < m_lo:
                m_lo = mm
            if mm > m_hi:
                m_hi = mm
            if mm >= stop_level and len(m_samples) > int(round(recWin / dt_sample)):
                break
    m = np.asarray(m_samples)
    det = _detect(m, dt_sample, theta, recThr, recWin)
    det["m_lo"], det["m_hi"] = float(m_lo), float(m_hi)

    # canonical PR#41 breath period + capture phase (measured-data convention)
    cap_idx = (int(round(det["t_capture"] / dt_sample)) if det["captured"] else None)
    pre = m[:cap_idx] if cap_idx is not None else m
    Tb = phase = None
    n_peaks_c = 0
    if len(pre) >= int(round(6.0 / dt_sample)):
        peaks = detect_peaks(np.asarray(pre, float), dt_sample)
        n_peaks_c = int(len(peaks))
        if len(peaks) - 1 >= MIN_CYCLES:
            ptimes = peaks * dt_sample
            Tb = float(np.median(np.diff(ptimes)))
            if det["captured"]:
                frac = ((det["t_capture"] - ptimes[-1]) / Tb) % 1.0
                phase = float(2 * np.pi * frac)
    det["Tb_canon"] = Tb
    det["n_peaks_canon"] = n_peaks_c
    det["capture_phase_canon"] = phase
    det.pop("peak_vals", None)
    return det


def _run_chunk(args):
    """Worker: one (fc, cell_index, N, j0, j1) chunk. fc must be picklable
    (module-level function or functools.partial of one)."""
    fc, ci, N, j0, j1 = args
    ics = ICS[N]
    rows = []
    for j in range(j0, j1):
        ic = ics[j % len(ics)]
        seed = SEED_BASE + ci * 100_000 + N * 1000 + j
        r = em_run_corr(list(ic), fc, N, seed)
        rows.append(dict(
            ci=ci, N=N, j=j, seed=seed,
            captured=bool(r["captured"]), t_capture=r["t_capture"],
            Tb_canon=r["Tb_canon"], n_peaks_canon=r["n_peaks_canon"],
            capture_phase_canon=r["capture_phase_canon"],
            breath_period_raw=r["breath_period"],
            capture_phase_raw=r["capture_phase"],
            m_lo=r["m_lo"], m_hi=r["m_hi"],
        ))
    return rows


def run_sweep(fc, cell_index=0, n_real=N_REAL, ns=NS, workers=9, chunk=25,
              progress=True):
    """Full sweep of one correction over ns x n_real. Returns rows grouped by
    N: {N: [row, ...]}. Deterministic given (fc, cell_index): seeds are
    SEED_BASE + cell_index*100_000 + N*1000 + j."""
    t0 = time.time()
    tasks = [(fc, cell_index, N, j0, min(j0 + chunk, n_real))
             for N in ns for j0 in range(0, n_real, chunk)]
    by_N = {N: [] for N in ns}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for i, rows in enumerate(ex.map(_run_chunk, tasks)):
            for r in rows:
                by_N[r["N"]].append(r)
            if progress and (i + 1) % 8 == 0:
                print(f"    {i + 1}/{len(tasks)} chunks ({time.time() - t0:.0f}s)")
    for N in ns:
        by_N[N].sort(key=lambda r: r["j"])
    return by_N


def zero_corr_gate(workers=9, n_real=N_REAL):
    """With f_corr == 0 the harness must reproduce the committed deterministic
    per-IC capture times to < 0.3% median per N (the Appendix B c=0 gate).
    Asserts and returns the per-N gate numbers."""
    from f_corr_zero import f_corr_zero
    by_N = run_sweep(f_corr_zero, cell_index=99, n_real=n_real, progress=False,
                     workers=workers)
    gate = {}
    for N in NS:
        rels = []
        for r in by_N[N]:
            d = DET[N][r["j"] % len(DET[N])]
            if r["captured"] and np.isfinite(d) and d > 0:
                rels.append(abs(r["t_capture"] - d) / d)
        med = float(np.median(rels))
        gate[str(N)] = {"median_rel": med, "n": len(rels), "passes": bool(med < 0.003)}
        assert med < 0.003, f"zero-correction gate FAILED at N={N}: median rel {med:.4%}"
    return gate


def save_rows(by_N, path):
    path.parent.mkdir(exist_ok=True)
    with open(path, "w") as f:
        for N in sorted(by_N):
            for r in by_N[N]:
                f.write(json.dumps(r) + "\n")
