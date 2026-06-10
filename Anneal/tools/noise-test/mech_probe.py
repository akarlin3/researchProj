"""Mechanism probe 2a — multiplicative, breath-locked, N-INDEPENDENT noise on the
reduced flow (A=0.5, beta=0.05).

The Discussion's open problem names "structured or multiplicative fluctuations
tied to the breath cycle" as the surviving candidate for the N-independent 3.2x
prolongation that additive c/sqrt(N) noise (Appendix B) failed to reproduce.
This script tests that candidate directly.

Noise model (pre-stated): Euler-Maruyama on (rho1, rho2, psi) with
    dX = f(X) dt + sigma_eff(X) dW,   sigma_eff = c_m * e(t),
    e(t) = 1 - min(rho1, rho2)   (the breath envelope: noise strongest
    mid-breath, vanishing at capture),
applied to all three collective variables with the SAME reflection of rho into
[0,1] as Appendix B (em_core). The amplitude c_m carries NO N factor: the noise
is N-independent by construction; N enters only through the seed-mapped IC
ensemble, exactly as in the deterministic reduced flow.

PRE-STATED PASS CONDITION (written before the campaign ran; the mechanism
"reproduces the puzzle" iff at the physically estimated amplitude):
  (1) prolongation factor ~3.2: median across N in [2.7, 3.7], where the per-N
      factor = median capture / per-N deterministic (DOP853) reference;
  (2) ~N-independent: CV of the per-N factor across N in {8,16,32,64} < 0.15;
  (3) breath-phase locking survives: Rayleigh p < 0.05 of capture phases
      (canonical PR#41 breath detector convention, as in the measured data)
      at every N;
  (4) k_cyc > 1: censored-Weibull shape in breath cycles (n_c = t/T_b), with
      95% profile-likelihood CI excluding 1, at every N.
Anything else = partial/fail, reported per condition with numbers.

Physical amplitude (pre-stated estimator): from the SAME baseline source as
Appendix B's c~0.05 (absorption_results/phase_traces.jsonl, A=0.5,
N in {8,16,32,64}, 100 traces each): per run, the high-frequency residual of
R_incoh about its 2-s moving average (the canonical breath smoother) over the
pre-absorption window gives sigma_HF; the envelope is e(t) = 1 - smooth(R);
c_m_run = sigma_HF / mean(e). Appendix-B consistency check: sigma_HF*sqrt(N)
should reproduce c ~ 0.05. The headline physical amplitude is the POOLED MEDIAN
of c_m_run over the four N (per-N medians reported alongside; if they are
N-dependent that is itself evidence against an N-independent multiplicative
amplitude and is reported as such).

Sweep: c_m = mult * c_m_phys, mult in {0.5, 1, 2, 4, 8, 16} (6 log-spaced
values bracketing the physical estimate), plus the mult=0 gate column;
N in {8,16,32,64}; 200 realizations per cell from the section-5.4 seed-mapped
ICs (reduced_results/reduced_runs.jsonl, identical to noise_sweep's ensemble —
asserted). Seeds: 8_000_000 + ci*100_000 + N*1000 + j (ci = index in the mult
list including 0), disjoint from the Appendix-B campaign (7_000_000 base).
EM dt = 0.01, t_max = 2000 (censoring counted and reported).

Gates before trusting (pre-stated):
  (i)  c_m = 0 must reproduce the deterministic (DOP853) capture to < 0.3%
       median per-IC (same gate as Appendix B);
  (ii) recomputed deterministic per-N reference must match the committed
       Appendix-B anchors (noise_results/noise_results.json) exactly;
  (iii) trajectory sanity: m stays in [0,1] (reflection), capture detection
       fires; per-run min/max recorded.

Outputs: noise_results/mech_probe_runs.jsonl (every run) and
noise_results/mech_probe_results.json (gates, estimate, per-cell matrix,
pass evaluation). Run: python3 tools/noise-test/mech_probe.py
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
sys.path.insert(0, str(HERE))
import reduced_core as rc  # noqa: E402
import noise_sweep as ns  # noqa: E402  (reuses CFG, ICS, deterministic_ref)
from em_core import _detect, _reflect  # noqa: E402
from _reuse import load_funcs  # noqa: E402

CFG = ns.CFG
A = ns.A                    # 0.5
BETA = ns.BETA              # 0.05
NS = ns.NS                  # [8, 16, 32, 64]
N_REAL = ns.N_REAL          # 200
DT = 0.01
T_MAX = 2000.0
SEED_BASE = 8_000_000
MULTS = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0]   # mult=0 is the gate column
OUT = ROOT / "noise_results"
OUT.mkdir(exist_ok=True)

# canonical PR#41 breath machinery (verbatim functions via the AST loader)
ABS = load_funcs(ROOT / "tools/absorption-recampaign/analysis.py")
BR = {"smoothWindowSec": 2.0, "minProminenceFrac": 0.1, "minPeakSepFrac": 0.5,
      "minCyclesForPhase": 2, "autoPeriodFloorSec": 8.0}
ABS["BR"] = BR
MIN_CYCLES = BR["minCyclesForPhase"]
detect_peaks = ABS["detect_peaks"]
moving_average = ABS["moving_average"]
rayleigh = ABS["rayleigh"]

# primary Weibull fitter: anneal-hazard profile-likelihood CI
import importlib.util  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "ah_survival", ROOT / "anneal-hazard/src/survival.py")
ah_survival = importlib.util.module_from_spec(_spec)
sys.modules["ah_survival"] = ah_survival  # dataclasses need the module registered
_spec.loader.exec_module(ah_survival)
fit_weibull_ci = ah_survival.fit_weibull

# ICs: reduced_results/reduced_runs.jsonl rows, sorted by seed, first N_REAL —
# asserted identical to noise_sweep's ensemble (transient_results features).
_rr = [json.loads(s) for s in (ROOT / "reduced_results/reduced_runs.jsonl")
       .read_text().splitlines() if s.strip()]
ICS = {}
for N in NS:
    rows = sorted([r for r in _rr if r["N"] == N], key=lambda r: r["seed"])[:N_REAL]
    ICS[N] = [(r["Rsync0"], r["Rincoh0"], r["dphi0"]) for r in rows]
for N in NS:
    assert ICS[N] == ns.ICS[N], f"IC ensemble mismatch at N={N}"


def em_run_mech(state0, c_m, seed, dt=DT, t_max=T_MAX):
    """One EM realization with multiplicative breath-locked noise
    sigma_eff = c_m * (1 - min(rho1, rho2)). Same drift (rhs_3d), reflection,
    sampling cadence, early stop and detection as em_core.em_run; only the
    noise amplitude differs (state-dependent, N-free)."""
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
    m_lo, m_hi = m_samples[0], m_samples[0]
    stop_level = 0.97
    for step in range(1, n_steps + 1):
        f = rc.rhs_3d(x, p)
        if c_m > 0:
            e = 1.0 - min(x[0], x[1])
            e = 0.0 if e < 0.0 else (1.0 if e > 1.0 else e)
            xi = rng.standard_normal(3)
            x = x + dt * np.asarray(f) + (c_m * e) * sqdt * xi
        else:
            x = x + dt * np.asarray(f)
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

    # canonical PR#41 breath period + capture phase (measured-data convention:
    # peaks on the smoothed pre-capture window, Tb = median peak spacing,
    # phase = 2*pi * frac((t_cap - last_peak)/Tb))
    cap_idx = (int(round(det["t_capture"] / dt_sample))
               if det["captured"] else None)
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


def run_chunk(args):
    """Worker: one (ci, N, j0, j1) chunk of realizations. Returns per-run rows."""
    ci, c_m, N, j0, j1 = args
    ics = ICS[N]
    rows = []
    for j in range(j0, j1):
        ic = ics[j % len(ics)]
        seed = SEED_BASE + ci * 100_000 + N * 1000 + j
        r = em_run_mech(list(ic), c_m, seed)
        rows.append(dict(
            ci=ci, mult=MULTS[ci], c_m=c_m, N=N, j=j, seed=seed,
            captured=bool(r["captured"]),
            t_capture=r["t_capture"],
            Tb_canon=r["Tb_canon"], n_peaks_canon=r["n_peaks_canon"],
            capture_phase_canon=r["capture_phase_canon"],
            breath_period_raw=r["breath_period"],
            capture_phase_raw=r["capture_phase"],
            m_lo=r["m_lo"], m_hi=r["m_hi"],
        ))
    return rows


def summarize_cell(rows, det_med):
    """Per-cell summary with censoring-aware fits (event=captured, censored at
    T_MAX). k (time) and k_cyc (cycles) via the anneal-hazard profile-CI fitter."""
    n = len(rows)
    cap = np.array([r["captured"] for r in rows], bool)
    t = np.array([r["t_capture"] if r["captured"] else T_MAX for r in rows], float)
    ev = cap.astype(int)
    n_cens = int(n - cap.sum())
    cap_frac = float(cap.mean())
    med_all = float(np.median(t))           # censoring-aware iff cap_frac > 0.5
    med_censored = bool(cap_frac <= 0.5)
    prolong = med_all / det_med if det_med and np.isfinite(det_med) else float("nan")

    phis = np.array([r["capture_phase_canon"] for r in rows
                     if r["capture_phase_canon"] is not None], float)
    n_ray, mphase, Rbar, z, p_ray = rayleigh(phis) if len(phis) else (
        0, np.nan, np.nan, np.nan, np.nan)

    wt = fit_weibull_ci(np.maximum(t, 0.05), ev) if ev.sum() >= 5 else None

    tb = np.array([r["Tb_canon"] if r["Tb_canon"] else np.nan for r in rows], float)
    ok = np.isfinite(tb) & (tb > 0)
    nc = t[ok] / tb[ok]
    evc = ev[ok]
    n_no_tb = int(n - ok.sum())
    wc = fit_weibull_ci(np.maximum(nc, 1e-3), evc) if evc.sum() >= 5 else None

    med_cyc = float(np.median(nc[evc == 1])) if (evc == 1).sum() else float("nan")
    return dict(
        n=n, n_captured=int(cap.sum()), n_censored=n_cens, capture_frac=cap_frac,
        median_capture=med_all, median_is_censored=med_censored,
        det_ref=det_med, prolongation=prolong,
        rayleigh_n=int(n_ray), rayleigh_Rbar=float(Rbar), rayleigh_p=float(p_ray),
        weibull_k=(float(wt.k) if wt else float("nan")),
        weibull_k_ci=([float(wt.k_ci[0]), float(wt.k_ci[1])] if wt else None),
        k_cyc=(float(wc.k) if wc else float("nan")),
        k_cyc_ci=([float(wc.k_ci[0]), float(wc.k_ci[1])] if wc else None),
        n_no_Tb=n_no_tb, median_cycles=med_cyc,
        median_Tb=float(np.nanmedian(tb)) if ok.any() else float("nan"),
        m_lo_min=float(min(r["m_lo"] for r in rows)),
        m_hi_max=float(max(r["m_hi"] for r in rows)),
    )


def physical_amplitude():
    """c_m estimate from the Appendix-B baseline source (phase traces)."""
    per_run = {N: [] for N in NS}
    with open(ROOT / "absorption_results/phase_traces.jsonl") as f:
        for line in f:
            d = json.loads(line)
            if d["A"] != A or d["N"] not in NS or d["abs_censored"]:
                continue
            if d.get("absIndex", -1) < 0:
                continue
            dt = d["sampleDt"]
            r = np.asarray(d["R_incoh"], float)[:d["absIndex"]]
            if len(r) < 100:
                continue
            sm = np.asarray(moving_average(r, round(BR["smoothWindowSec"] / dt)), float)
            k = 25  # trim smoother edge effects
            resid = (r - sm)[k:-k]
            env = (1.0 - sm)[k:-k]
            sig = float(np.std(resid))
            emean = float(np.mean(env))
            per_run[d["N"]].append((sig, emean, sig / emean))
    table = {}
    for N in NS:
        a = np.array(per_run[N])
        table[N] = dict(
            n_traces=int(len(a)),
            sigma_HF_median=float(np.median(a[:, 0])),
            sigma_HF_x_sqrtN=float(np.median(a[:, 0]) * np.sqrt(N)),
            env_mean_median=float(np.median(a[:, 1])),
            c_m_median=float(np.median(a[:, 2])),
        )
    pooled = np.concatenate([np.array(per_run[N])[:, 2] for N in NS])
    c_m_phys = float(np.median(pooled))
    cms = [table[N]["c_m_median"] for N in NS]
    return c_m_phys, table, float(np.std(cms) / np.mean(cms))


def evaluate(cells_by_N):
    """Apply the four pre-stated conditions to one amplitude's per-N cells."""
    facs = [cells_by_N[N]["prolongation"] for N in NS]
    med_fac = float(np.median(facs))
    cv = float(np.std(facs) / np.mean(facs)) if np.mean(facs) else float("nan")
    c1 = bool(2.7 <= med_fac <= 3.7)
    c2 = bool(cv < 0.15)
    c3 = bool(all(np.isfinite(cells_by_N[N]["rayleigh_p"])
                  and cells_by_N[N]["rayleigh_p"] < 0.05 for N in NS))
    kls = [(cells_by_N[N]["k_cyc_ci"][0] if cells_by_N[N]["k_cyc_ci"] else
            float("nan")) for N in NS]
    c4 = bool(all(np.isfinite(k) and k > 1.0 for k in kls))
    return dict(per_N_factor={str(N): facs[i] for i, N in enumerate(NS)},
                median_factor=med_fac, cv_factor=cv,
                cond1_prolong_2p7_3p7=c1, cond2_cv_lt_0p15=c2,
                cond3_rayleigh_all_N=c3, cond4_kcyc_gt1_all_N=c4,
                all_pass=bool(c1 and c2 and c3 and c4))


def main():
    t0 = time.time()
    print("[1/5] physical amplitude estimate from phase traces ...")
    c_m_phys, est_table, cm_cv = physical_amplitude()
    for N in NS:
        e = est_table[N]
        print(f"  N={N:2d}: sigma_HF={e['sigma_HF_median']:.4f} "
              f"(x sqrtN = {e['sigma_HF_x_sqrtN']:.4f}; Appendix-B c~0.05) "
              f"<e>={e['env_mean_median']:.3f}  c_m={e['c_m_median']:.4f}")
    print(f"  pooled-median c_m_phys = {c_m_phys:.4f}  "
          f"(per-N c_m CV = {cm_cv:.3f} -> "
          f"{'~N-independent' if cm_cv < 0.15 else 'N-DEPENDENT'})")

    print("[2/5] deterministic per-N reference (DOP853) ...")
    det_ref, det_per_ic = ns.deterministic_ref()
    anchors = json.load(open(OUT / "noise_results.json"))["deterministic_ref"]
    anchor_ok = all(abs(det_ref[N] - anchors[str(N)]) < 1e-9 for N in NS)
    for N in NS:
        print(f"  N={N:2d}: {det_ref[N]:.2f}s (committed anchor {anchors[str(N)]}; "
              f"{'OK' if abs(det_ref[N]-anchors[str(N)])<1e-9 else 'MISMATCH'})")
    if not anchor_ok:
        print("  !! deterministic anchors do not match Appendix B — aborting")
        sys.exit(1)

    cms = [m * c_m_phys for m in MULTS]
    tasks = []
    chunk = 50
    for ci, c_m in enumerate(cms):
        for N in NS:
            for j0 in range(0, N_REAL, chunk):
                tasks.append((ci, c_m, N, j0, min(j0 + chunk, N_REAL)))
    print(f"[3/5] EM campaign: {len(MULTS)}x{len(NS)} cells x {N_REAL} runs "
          f"(dt={DT}, t_max={T_MAX}, {len(tasks)} chunks) ...")
    all_rows = []
    with ProcessPoolExecutor(max_workers=9) as ex:
        for i, rows in enumerate(ex.map(run_chunk, tasks)):
            all_rows.extend(rows)
            if (i + 1) % 14 == 0:
                print(f"  {i+1}/{len(tasks)} chunks ({time.time()-t0:.0f}s)")

    with open(OUT / "mech_probe_runs.jsonl", "w") as f:
        for r in all_rows:
            f.write(json.dumps(r) + "\n")

    print("[4/5] gates + per-cell summaries ...")
    # gate (i): c_m = 0 vs deterministic per-IC
    gate_c0 = {}
    c0_ok = True
    for N in NS:
        det = det_per_ic[N]
        rels = []
        for r in all_rows:
            if r["ci"] != 0 or r["N"] != N:
                continue
            d = det[r["j"] % len(det)]
            if r["captured"] and np.isfinite(d) and d > 0:
                rels.append(abs(r["t_capture"] - d) / d)
        rels = np.array(rels)
        med = float(np.median(rels))
        ok = med < 0.003
        c0_ok = c0_ok and ok
        gate_c0[str(N)] = dict(median_rel=med, n=len(rels), passes=bool(ok))
        print(f"  [gate c_m=0] N={N:2d}: median per-IC rel = {med*100:.3f}% "
              f"{'OK' if ok else 'FAIL'}")

    matrix = {}
    for ci, c_m in enumerate(cms):
        for N in NS:
            rows = [r for r in all_rows if r["ci"] == ci and r["N"] == N]
            matrix[(ci, N)] = summarize_cell(rows, det_ref[N])

    hdr = (f"{'mult':>5} {'c_m':>7} {'N':>3} {'capfrac':>7} {'med':>7} "
           f"{'prolong':>7} {'ray_p':>9} {'k_t':>5} {'k_cyc':>6} "
           f"{'kcyc_CI':>13} {'cens':>4}")
    print(hdr)
    for ci, c_m in enumerate(cms):
        for N in NS:
            m = matrix[(ci, N)]
            ci_s = (f"[{m['k_cyc_ci'][0]:.2f},{m['k_cyc_ci'][1]:.2f}]"
                    if m["k_cyc_ci"] else "--")
            print(f"{MULTS[ci]:>5.1f} {c_m:>7.4f} {N:>3d} {m['capture_frac']:>7.2f} "
                  f"{m['median_capture']:>7.1f} {m['prolongation']:>7.2f} "
                  f"{m['rayleigh_p']:>9.1e} {m['weibull_k']:>5.2f} "
                  f"{m['k_cyc']:>6.2f} {ci_s:>13} {m['n_censored']:>4d}")
        print("-" * len(hdr))

    print("[5/5] pass evaluation against the pre-stated conditions ...")
    evals = {}
    for ci in range(1, len(MULTS)):
        cells = {N: matrix[(ci, N)] for N in NS}
        ev = evaluate(cells)
        evals[str(MULTS[ci])] = ev
        print(f"  mult={MULTS[ci]:>4.1f} (c_m={cms[ci]:.4f}): "
              f"med_factor={ev['median_factor']:.2f} CV={ev['cv_factor']:.3f} "
              f"c1={ev['cond1_prolong_2p7_3p7']} c2={ev['cond2_cv_lt_0p15']} "
              f"c3={ev['cond3_rayleigh_all_N']} c4={ev['cond4_kcyc_gt1_all_N']} "
              f"-> {'PASS' if ev['all_pass'] else 'no'}")
    phys_ev = evals["1.0"]
    any_pass = [m for m, e in evals.items() if e["all_pass"]]

    results = dict(
        pre_stated_pass_condition=(
            "At the physically estimated amplitude: (1) prolongation factor "
            "(median capture / per-N deterministic reference), median across "
            "N in {8,16,32,64}, within [2.7,3.7]; (2) CV of the per-N factor "
            "across N < 0.15; (3) Rayleigh p < 0.05 of capture phases "
            "(canonical breath detector) at every N; (4) censored-Weibull "
            "k_cyc > 1 with 95% profile CI excluding 1 at every N."),
        model=("sigma_eff = c_m * e, e = 1 - min(rho1,rho2), applied to "
               "(rho1,rho2,psi), Euler-Maruyama, reflection of rho into [0,1]; "
               "amplitude N-independent by construction"),
        A=A, beta=BETA, dt=DT, t_max=T_MAX, n_realizations=N_REAL, Ns=NS,
        seed_scheme=f"seed = {SEED_BASE} + ci*100000 + N*1000 + j",
        mults=MULTS, c_m_values={str(MULTS[i]): cms[i] for i in range(len(MULTS))},
        c_m_phys=c_m_phys,
        physical_estimate=dict(
            source="absorption_results/phase_traces.jsonl (A=0.5, pre-absorption)",
            method=("sigma_HF = std(R_incoh - 2s moving average); "
                    "e = 1 - smooth; c_m = sigma_HF/mean(e) per run; "
                    "headline = pooled median over N"),
            per_N={str(N): est_table[N] for N in NS},
            per_N_cm_cv=cm_cv,
            appendixB_consistency=(
                "sigma_HF*sqrt(N) = "
                + ", ".join(f"{est_table[N]['sigma_HF_x_sqrtN']:.3f}" for N in NS)
                + " vs Appendix B c ~ 0.05")),
        gates=dict(cm0_vs_deterministic=gate_c0, cm0_pass=bool(c0_ok),
                   det_ref_anchor_match=bool(anchor_ok)),
        deterministic_ref={str(N): det_ref[N] for N in NS},
        matrix={f"mult{MULTS[ci]}_N{N}": matrix[(ci, N)]
                for ci in range(len(MULTS)) for N in NS},
        pass_eval_by_mult=evals,
        physical_amplitude_eval=phys_ev,
        amplitudes_passing_all=any_pass,
        runtime_s=float(time.time() - t0),
        command="python3 tools/noise-test/mech_probe.py",
    )
    (OUT / "mech_probe_results.json").write_text(
        json.dumps(results, indent=2, default=float))
    print(f"\nWrote {OUT/'mech_probe_results.json'} and mech_probe_runs.jsonl "
          f"({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
