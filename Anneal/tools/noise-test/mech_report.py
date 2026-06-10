"""Assemble paper/revision-data-gated/results_mech.json from the two
re-injection mechanism probes:

  2a  noise_results/mech_probe_results.json (+ mech_probe_runs.jsonl)
      — multiplicative breath-locked N-independent noise on the reduced flow
        (tools/noise-test/mech_probe.py).
  2b  manifold_results/ws_decomposition_N{8,16,32,64}.jsonl
      — actual-constants vs constants-uniformized (Poisson-projected) partner
        ICs under the absorption campaign's exact integrator + labeler
        (tools/manifold-probe/ws_decomposition.mjs).

PRE-STATED PASS CONDITIONS (identical four-condition frame for both
mechanisms; written before the campaigns ran):
  (1) prolongation factor ~3.2 (median across N in [2.7, 3.7]);
  (2) ~N-independent (CV of the factor across N in {8,16,32,64} < 0.15);
  (3) breath-phase locking survives (Rayleigh p < 0.05 at every N);
  (4) k_cyc > 1 (censored-Weibull shape in breath cycles, 95% profile CI
      excluding 1, at every N).
For 2b the factor decomposition additionally requires (pre-set thresholds,
fixed before the full campaign): "partner collapses to reduced" iff
median(t_partner)/median(t_reduced) <= 1.3 at every N, and "actual stays
measured" iff the actual runs satisfy (1)-(4). If instead t_partner ~ t_actual,
the constants do not carry the prolongation (consistent with the prior D0
null) and 2b FAILS as a mechanism (reported straight).

Also re-verifies a few reduced-flow rows (t_capture) by direct DOP853
re-integration, and computes the measured per-N prolongation context.

Run: python3 tools/noise-test/mech_report.py
"""
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "tools/reduced-ode"))
sys.path.insert(0, str(ROOT / "tools/beta010-slice"))
import reduced_core as rc  # noqa: E402
from _reuse import load_funcs  # noqa: E402

NS = [8, 16, 32, 64]
ABS = load_funcs(ROOT / "tools/absorption-recampaign/analysis.py")
rayleigh = ABS["rayleigh"]

import importlib.util  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "ah_survival", ROOT / "anneal-hazard/src/survival.py")
ah_survival = importlib.util.module_from_spec(_spec)
sys.modules["ah_survival"] = ah_survival
_spec.loader.exec_module(ah_survival)
fit_weibull_ci = ah_survival.fit_weibull

T_MAX = 2000.0
BOOT_SEED = 20260610
PARTNER_COLLAPSE_MAX = 1.3   # pre-set: partner/reduced <= this at every N
RR = [json.loads(s) for s in (ROOT / "reduced_results/reduced_runs.jsonl")
      .read_text().splitlines() if s.strip()]


def boot_ratio_ci(num, den, n_boot=2000, seed=BOOT_SEED):
    rng = np.random.default_rng(seed)
    num, den = np.asarray(num, float), np.asarray(den, float)
    vals = []
    for _ in range(n_boot):
        a = rng.choice(num, len(num))
        b = rng.choice(den, len(den))
        vals.append(np.median(a) / np.median(b))
    return [float(np.quantile(vals, 0.025)), float(np.quantile(vals, 0.975))]


def reduced_verify(rows_per_N=3):
    """Re-integrate a few reduced rows; compare to committed t_capture."""
    cfg = rc.load_config()
    rc.set_config(cfg)
    out = []
    for N in NS:
        rows = sorted([r for r in RR if r["N"] == N], key=lambda r: r["seed"])
        for r in rows[:rows_per_N]:
            res = rc.reduced_run_3d([r["Rsync0"], r["Rincoh0"], r["dphi0"]],
                                    rc.Params(A=0.5, beta=0.05), T_MAX, cfg)
            out.append(dict(N=N, seed=r["seed"],
                            committed=r["t_capture"], recomputed=res["t_capture"],
                            match=bool(res["t_capture"] == r["t_capture"])))
    return out


def analyze_2b():
    rows = []
    for N in NS:
        p = ROOT / f"manifold_results/ws_decomposition_N{N}.jsonl"
        rows += [json.loads(l) for l in p.read_text().splitlines() if l.strip()]
    act = {(r["N"], r["seed"]): r for r in rows if r["kind"] == "actual"}
    par = {(r["N"], r["seed"]): r for r in rows if r["kind"] == "partner"}
    assert set(act) == set(par)

    # ---- gates ----
    ic_dev_max = max(r["ic_dev"] for r in act.values())
    tabs_match_frac = float(np.mean([r["t_abs_match"] for r in act.values()]))
    tabs_dev_max = max(abs(r["t_abs"] - r["t_abs_meas"]) for r in act.values())
    z1_err_max = max(max(r["z1_err_pop1"], r["z1_err_pop2"])
                     for r in par.values())
    d0_gate = {}
    for N in NS:
        a2 = [r["d0_pop2"] for r in act.values() if r["N"] == N]
        p2 = [r["d0_pop2"] for r in par.values() if r["N"] == N]
        a1 = [r["d0_pop1"] for r in act.values() if r["N"] == N]
        p1 = [r["d0_pop1"] for r in par.values() if r["N"] == N]
        d0_gate[str(N)] = dict(
            d0_incoh_actual_median=float(np.median(a2)),
            d0_incoh_partner_median=float(np.median(p2)),
            d0_sync_actual_median=float(np.median(a1)),
            d0_sync_partner_median=float(np.median(p1)))

    # ---- per-N outcome table ----
    table = {}
    for N in NS:
        keys = sorted(k for k in act if k[0] == N)
        ta = np.array([act[k]["t_abs"] for k in keys])
        ca = np.array([act[k]["abs_censored"] for k in keys], bool)
        tp = np.array([par[k]["t_abs"] for k in keys])
        cp = np.array([par[k]["abs_censored"] for k in keys], bool)
        tr = np.array([act[k]["t_capture_reduced"] for k in keys], float)
        ratio_pa = [par[k]["t_abs"] / act[k]["t_abs"] for k in keys
                    if act[k]["t_abs"] > 0]
        med_a, med_p, med_r = (float(np.median(ta)), float(np.median(tp)),
                               float(np.nanmedian(tr)))
        # conditions (3)/(4) on ACTUAL runs
        phis = np.array([act[k]["capture_phase"] for k in keys
                         if act[k]["capture_phase"] is not None], float)
        n_ray, mp_, Rbar, z, p_ray = rayleigh(phis)
        nc, evc = [], []
        for k in keys:
            r = act[k]
            if r["T_b"]:
                nc.append(max(r["t_abs"] / r["T_b"], 1e-3))
                evc.append(0 if r["abs_censored"] else 1)
        wc = fit_weibull_ci(np.array(nc), np.array(evc))
        # localization: where in the run does the prolongation accrue?
        tg_a = float(np.median([act[k]["t_graze"] for k in keys]))
        tg_p = float(np.median([par[k]["t_graze"] for k in keys]))
        ng_a = float(np.mean([act[k]["n_grazes_before_abs"] for k in keys]))
        ng_p = float(np.mean([par[k]["n_grazes_before_abs"] for k in keys]))
        table[str(N)] = dict(
            n_pairs=len(keys),
            median_t_actual=med_a, censored_actual=int(ca.sum()),
            median_t_partner=med_p, censored_partner=int(cp.sum()),
            median_t_reduced=med_r,
            median_t_graze_actual=tg_a, median_t_graze_partner=tg_p,
            graze_factor_actual=tg_a / med_r, graze_factor_partner=tg_p / med_r,
            mean_n_grazes_actual=ng_a, mean_n_grazes_partner=ng_p,
            factor_actual_vs_reduced=med_a / med_r,
            factor_partner_vs_reduced=med_p / med_r,
            ratio_actual_to_partner=med_a / med_p,
            ratio_actual_to_partner_ci=boot_ratio_ci(ta, tp),
            paired_median_partner_over_actual=float(np.median(ratio_pa)),
            rayleigh_n=int(n_ray), rayleigh_Rbar=float(Rbar),
            rayleigh_p=float(p_ray),
            k_cyc=float(wc.k), k_cyc_ci=[float(wc.k_ci[0]), float(wc.k_ci[1])],
            n_kcyc=len(nc),
        )

    fa = [table[str(N)]["factor_actual_vs_reduced"] for N in NS]
    fp = [table[str(N)]["factor_partner_vs_reduced"] for N in NS]
    med_fa, med_fp = float(np.median(fa)), float(np.median(fp))
    cv_fa = float(np.std(fa) / np.mean(fa))
    c1 = bool(2.7 <= med_fa <= 3.7)
    c2 = bool(cv_fa < 0.15)
    c3 = bool(all(table[str(N)]["rayleigh_p"] < 0.05 for N in NS))
    c4 = bool(all(table[str(N)]["k_cyc_ci"][0] > 1.0 for N in NS))
    partner_collapses = bool(all(
        table[str(N)]["factor_partner_vs_reduced"] <= PARTNER_COLLAPSE_MAX
        for N in NS))
    mech_pass = bool(partner_collapses and c1 and c2 and c3 and c4)

    return dict(
        gates=dict(
            ic_reproduction_max_abs_dev=float(ic_dev_max),
            t_abs_rerun_match_fraction=tabs_match_frac,
            t_abs_rerun_max_abs_dev=float(tabs_dev_max),
            partner_z1_match_max_err=float(z1_err_max),
            d0_actual_vs_partner=d0_gate,
        ),
        per_N=table,
        factor_actual_median=med_fa, factor_actual_cv=cv_fa,
        factor_partner_median=med_fp,
        partner_collapse_threshold=PARTNER_COLLAPSE_MAX,
        partner_collapses_to_reduced=partner_collapses,
        cond1_prolong_2p7_3p7=c1, cond2_cv_lt_0p15=c2,
        cond3_rayleigh_all_N=c3, cond4_kcyc_gt1_all_N=c4,
        all_pass=mech_pass,
    )


def measured_context():
    """Measured per-N prolongation (t_abs_meas / reduced t_capture medians)."""
    out = {}
    facs = []
    for N in NS:
        rows = [r for r in RR if r["N"] == N]
        mm = float(np.median([r["t_abs_meas"] for r in rows]))
        mr = float(np.nanmedian([r["t_capture"] if r["captured"] else np.nan
                                 for r in rows]))
        out[str(N)] = dict(measured_median=mm, reduced_median=mr, factor=mm / mr)
        facs.append(mm / mr)
    out["median_factor"] = float(np.median(facs))
    out["cv_factor"] = float(np.std(facs) / np.mean(facs))
    return out


def main():
    t0 = time.time()
    res2a = json.load(open(ROOT / "noise_results/mech_probe_results.json"))
    res2b = analyze_2b()
    rv = reduced_verify()
    ctx = measured_context()

    git_hash = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True,
                              cwd=ROOT).stdout.strip()

    out = dict(
        title=("Second re-injection mechanism probe: (2a) multiplicative "
               "breath-locked N-independent noise on the reduced flow; "
               "(2b) off-manifold Watanabe-Strogatz constants "
               "(actual vs constants-uniformized partner ICs)"),
        pre_stated_pass_condition=res2a["pre_stated_pass_condition"],
        pre_stated_2b_decomposition=(
            "partner collapses to reduced iff median(t_partner)/"
            f"median(t_reduced) <= {PARTNER_COLLAPSE_MAX} at every N; "
            "2b passes iff partner collapses AND actual runs satisfy (1)-(4)"),
        measured_context=ctx,
        reduced_row_verification=rv,
        experiment_2a=res2a,
        experiment_2b=res2b,
        contradictions_with_manuscript=[
            ("Sec. 5.4 sentence 'the off-manifold constants ... are precisely "
             "where the residual 3.2x prolongation must originate' is REFUTED "
             "by experiment 2b: constants-uniformized partners (Z1 matched to "
             "9e-16, incoherent D0 at the splay floor) are prolonged by the "
             "same factor as the actual seeds (actual/partner medians "
             "0.98-1.08, bootstrap CIs all spanning 1). Minimal fix proposed "
             "in snippet_mech.tex block (3)."),
            ("Discussion sentence 'the same factor of 3.2 across the full 16x "
             "range of N' holds only against the pooled reduced median: under "
             "per-N referencing the measured factor is 3.44/3.16/3.03/1.98 "
             "(CV 0.19; flat at CV 0.057 over N=8-32) because the reduced "
             "ensemble's median capture crosses one extra breath cycle "
             "(43->64 s) at N>=48 while the measured plateau stays flat "
             "(CV 0.046). Qualifier proposed in snippet_mech.tex block (4)."),
        ],
        notes=[
            ("Pre-stated condition (2) (CV<0.15) and condition (3) (per-N "
             "Rayleigh p<0.05) are stricter than the measured system itself "
             "satisfies: the measured per-N factor has CV 0.19, and the "
             "published absorption analysis finds per-N Rayleigh significant "
             "only at N=8,16 (N=32: p=0.43, N=64: p=0.095; pooled p=2.1e-6) — "
             "experiment 2b's actual runs reproduce this (N=32: p=0.165, "
             "N=64: p=0.099)."),
        ],
        commands=[
            "python3 tools/noise-test/mech_probe.py",
            "node tools/manifold-probe/ws_decomposition.mjs --N {8,16,32,64} --nseeds 200",
            "python3 tools/noise-test/mech_report.py",
            "python3 paper_figures/fig_mech_probe.py",
        ],
        git_hash=git_hash,
        generated_unix=time.time(),
        report_runtime_s=float(time.time() - t0),
    )
    dest = ROOT / "paper/revision-data-gated"
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "results_mech.json").write_text(json.dumps(out, indent=2,
                                                       default=float))
    print(f"Wrote {dest/'results_mech.json'}")

    print("\n=== 2a (multiplicative breath-locked) at physical amplitude ===")
    ev = res2a["physical_amplitude_eval"]
    print(json.dumps(ev, indent=1))
    print("amplitudes passing all four:", res2a["amplitudes_passing_all"])
    print("\n=== 2b (WS constants decomposition) ===")
    for N in NS:
        t = res2b["per_N"][str(N)]
        print(f"  N={N:2d}: actual={t['median_t_actual']:.1f} "
              f"partner={t['median_t_partner']:.1f} reduced={t['median_t_reduced']:.1f} "
              f"fa={t['factor_actual_vs_reduced']:.2f} fp={t['factor_partner_vs_reduced']:.2f} "
              f"ray_p={t['rayleigh_p']:.1e} k_cyc={t['k_cyc']:.2f} "
              f"CI=[{t['k_cyc_ci'][0]:.2f},{t['k_cyc_ci'][1]:.2f}]")
    print("partner collapses to reduced:", res2b["partner_collapses_to_reduced"])
    print("2b all_pass:", res2b["all_pass"])
    print("gates:", json.dumps(res2b["gates"], indent=1)[:400])


if __name__ == "__main__":
    main()
