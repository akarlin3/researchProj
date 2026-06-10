"""Null self-test: feed the harness the ALREADY-EXCLUDED additive 1/sqrt(N)
collective noise (Appendix B) as f_corr and confirm the scorer rejects it,
reproducing the committed Appendix B result. This is the evidence that the
harness is not rigged to pass: a known-wrong mechanism must come out FAIL.

Gates run first:
  (a) anchor gate     — committed deterministic refs match Appendix B anchors,
                        DOP853 row recompute bit-identical;
  (b) zero-corr gate  — f_corr == 0 reproduces deterministic capture < 0.3%;
  (c) equivalence gate— with B = (c/sqrt(N)) I and zero drift, em_run_corr is
                        BIT-IDENTICAL to the Appendix B driver em_core.em_run
                        at the same seed (same rng stream, same arithmetic).

Then the additive null at c = 0.05 (the physical Appendix B amplitude) and
c = 0.2 (the strongest committed amplitude) is swept (200 runs/cell,
N in {8,16,32,64}) and scored against the four pre-committed conditions.

Hard assertions (the "known Appendix B result"):
  - all_pass is False at BOTH amplitudes;
  - at c = 0.2 the prolongation is N-DEPENDENT: cond2 fails (CV >= 0.15) and
    factor(N=8)/factor(N=64) > 1.5 (committed: 2.363/0.682 = 3.46);
  - at c = 0.05 there is no ~3x prolongation: cond1 fails.
Per-N factors are also reported next to the committed Appendix B values
(noise_results/noise_results.json; different seed base, so agreement is
statistical, not bit-level).

Run: python3 tools/kinetic-probe/null_test.py
Outputs: kinetic_results/score_null_c0.05.json, score_null_c0.2.json,
         null_selftest.json, runs_null_*.jsonl
"""
import json
import sys
import time
from functools import partial
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import harness as H  # noqa: E402
import scorer as S  # noqa: E402
from f_corr_zero import _ZERO3  # noqa: E402
from em_core import em_run  # noqa: E402  (path set by harness import)

ROOT = H.ROOT
OUT = H.OUT
C_PHYS = 0.05   # Appendix B physical amplitude
C_STRONG = 0.2  # strongest committed Appendix B amplitude

_EYE = np.eye(3)
_B_CACHE = {}


def f_corr_additive(rho1, rho2, psi, N, c):
    """Appendix B additive collective noise as a correction term:
    zero drift, B = (c/sqrt(N)) * I."""
    key = (c, N)
    B = _B_CACHE.get(key)
    if B is None:
        B = _B_CACHE.setdefault(key, (c / np.sqrt(N)) * _EYE)
    return _ZERO3, B


def equivalence_gate():
    """em_run_corr with the additive correction must be bit-identical to the
    Appendix B driver em_core.em_run at the same seed."""
    checks = []
    for N in (8, 64):
        for seed in (12345, 54321, 99999):
            ic = list(H.ICS[N][0])
            a = em_run(ic, H.A, H.BETA, C_PHYS, N, H.CFG, seed=seed,
                       dt=H.DT, t_max=H.T_MAX)
            b = H.em_run_corr(ic, partial(f_corr_additive, c=C_PHYS), N, seed)
            same = (a["captured"] == b["captured"]
                    and a["t_capture"] == b["t_capture"]
                    and a["breath_period"] == b["breath_period"])
            checks.append({"N": N, "seed": seed, "em_core_t": a["t_capture"],
                           "harness_t": b["t_capture"], "bit_identical": same})
            assert same, f"equivalence gate FAILED at N={N} seed={seed}: {a['t_capture']} vs {b['t_capture']}"
    return checks


def main():
    t0 = time.time()
    OUT.mkdir(exist_ok=True)

    print("[1/5] anchor gate (committed det refs vs Appendix B anchors) ...")
    g_anchor = H.anchor_gate()
    print(f"  medians match: {g_anchor['medians_match_appendixB_anchors']}; "
          f"row recompute worst dev = {g_anchor['row_recompute_worst_dev_s']:.1e} s")

    print("[2/5] zero-correction determinism gate (<0.3% per N) ...")
    g_zero = H.zero_corr_gate()
    for N in H.NS:
        g = g_zero[str(N)]
        print(f"  N={N:2d}: median per-IC rel = {g['median_rel'] * 100:.3f}%  "
              f"{'OK' if g['passes'] else 'FAIL'}")

    print("[3/5] bit-equivalence gate vs Appendix B driver (em_core.em_run) ...")
    g_equiv = equivalence_gate()
    print(f"  {len(g_equiv)}/6 (N, seed) pairs bit-identical")

    committed = json.load(open(ROOT / "noise_results/noise_results.json"))["matrix"]
    scores = {}
    for ci, c in ((1, C_PHYS), (2, C_STRONG)):
        print(f"[4/5] additive null sweep c={c} (4 N x {H.N_REAL} runs) ...")
        by_N = H.run_sweep(partial(f_corr_additive, c=c), cell_index=ci)
        H.save_rows(by_N, OUT / f"runs_null_c{c}.jsonl")
        sc = S.score(by_N, label=f"NULL additive 1/sqrt(N) noise, c={c}",
                     extra={"mechanism": "Appendix B additive collective noise "
                                         "(already excluded)",
                            "c": c, "seed_cell_index": ci})
        # committed Appendix B comparison (statistical: different seed base)
        comp = {}
        for N in H.NS:
            mine = sc["per_N_factor"][str(N)]
            ref = committed[f"c{c}_N{N}"]["ratio_to_det"]
            comp[str(N)] = {"harness_factor": mine, "committed_appendixB": ref,
                            "rel_dev": abs(mine - ref) / ref}
        sc["committed_appendixB_comparison"] = comp
        scores[c] = sc
        S.write_score(sc, OUT / f"score_null_c{c}.json")
        print(f"  c={c}: per-N factor "
              + " ".join(f"{sc['per_N_factor'][str(N)]:.2f}" for N in H.NS)
              + f"  (committed: "
              + " ".join(f"{comp[str(N)]['committed_appendixB']:.2f}" for N in H.NS)
              + f")  median={sc['median_factor']:.2f} CV={sc['cv_factor']:.3f}")
        print(f"  conditions: c1={sc['cond1_prolong_2p9_3p5']} "
              f"c2={sc['cond2_cv_lt_0p15']} c3={sc['cond3_rayleigh_all_N']} "
              f"c4={sc['cond4_kcyc_gt1_all_N']} -> verdict={sc['verdict'].upper()}")
        sec = sc["secondary_pattern_match"]
        print(f"  secondary (measured-pattern): "
              f"S1={sec['S1_factor_pattern']['passes']} "
              f"S2={sec['S2_phase_pattern']['passes']} "
              f"S3={sec['S3_kcyc_pattern']['passes']} "
              f"-> matches={sec['matches_measured_pattern']}")

    print("[5/5] hard assertions (harness must reject the known-wrong mechanism)")
    sp, ss = scores[C_PHYS], scores[C_STRONG]
    assert not sp["all_pass"], "null at physical amplitude must not pass"
    assert not ss["all_pass"], "null at strong amplitude must not pass"
    assert not sp["cond1_prolong_2p9_3p5"], \
        "null at physical amplitude must show no ~3x prolongation (cond1 fail)"
    assert not ss["cond2_cv_lt_0p15"], \
        "null at strong amplitude must be N-dependent (cond2 fail, CV >= 0.15)"
    r8 = ss["per_N_factor"]["8"]
    r64 = ss["per_N_factor"]["64"]
    assert r8 / r64 > 1.5, f"expected N-dependent decay of the factor, got {r8:.2f}/{r64:.2f}"
    assert not sp["secondary_pattern_match"]["matches_measured_pattern"], \
        "null at physical amplitude must not match the measured pattern"
    assert not ss["secondary_pattern_match"]["matches_measured_pattern"], \
        "null at strong amplitude must not match the measured pattern"
    print("  all assertions hold: scorer FAILS the additive mechanism on the "
          "primary conditions AND the secondary measured-pattern criterion, "
          "matching the committed Appendix B exclusion")

    selftest = {
        "purpose": "prove the harness rejects a known-wrong mechanism "
                   "(additive 1/sqrt(N) noise, Appendix B)",
        "gates": {"anchor": g_anchor, "zero_corr": g_zero,
                  "bit_equivalence_vs_em_core": g_equiv},
        "null_verdicts": {str(c): {k: scores[c][k] for k in
                                   ("per_N_factor", "median_factor", "cv_factor",
                                    "cond1_prolong_2p9_3p5", "cond2_cv_lt_0p15",
                                    "cond3_rayleigh_all_N", "cond4_kcyc_gt1_all_N",
                                    "all_pass", "verdict",
                                    "secondary_pattern_match",
                                    "committed_appendixB_comparison")}
                          for c in (C_PHYS, C_STRONG)},
        "conclusion": "harness correctly FAILS the additive mechanism at both "
                      "the physical and 4x amplitudes; N-independence rejection "
                      "confirmed at c=0.2",
        "runtime_s": time.time() - t0,
        "command": "python3 tools/kinetic-probe/null_test.py",
    }
    (OUT / "null_selftest.json").write_text(json.dumps(selftest, indent=2,
                                                       default=float))
    print(f"\nWrote {OUT}/score_null_c*.json and null_selftest.json "
          f"({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
