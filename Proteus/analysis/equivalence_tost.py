"""Equivalence / non-superiority test for the central claim (revision Fix A.1).

The original claim — "the PET-branch exposed-cleft rate (96/294 = 32.65%) is
indistinguishable from the random floor (12/28 = 42.86%; Fisher p=0.30)" — is a
failure to reject on a tiny (n=28) baseline. Absence of evidence is not equivalence.
Here we (1) read the ENLARGED random triad-bearing floor (>=200 triad-bearers,
analysis/enlarge_floor.py) and (2) replace Fisher with a formal test of the actual
scientific claim: the PET-branch rate is NOT enriched above the random floor.

We report, on the difference d = PET_rate - baseline_rate:
  * TOST equivalence within a PRE-REGISTERED margin DELTA_PP (= 10 percentage
    points): two one-sided tests; equivalence iff the 90% CI of d lies inside
    (-DELTA, +DELTA).
  * NON-SUPERIORITY (the load-bearing direction): one-sided test of H0: d >= DELTA.
    Rejecting it means the PET branch is not enriched above the floor by DELTA.
  * The UPPER bound of the one-sided 95% CI on d — the largest plausible enrichment.

Honest reporting: if PET is *below* the floor by more than DELTA, strict two-sided
equivalence fails — we say so (PET is depleted, not merely equal) and report the
non-superiority result, which is the claim the paper actually makes.

Reproduce, from the repo root:
    PYTHONPATH=src python analysis/equivalence_tost.py
"""
from __future__ import annotations

import json
import math
import os

# ---- single config constant: the pre-registered equivalence margin ---------- #
DELTA_PP = 10.0                 # equivalence/non-superiority margin, percentage points
DELTA = DELTA_PP / 100.0
ALPHA = 0.05
Z_A = 1.6448536269514722        # one-sided z at alpha=0.05 (for TOST / 90% CI)
LINE = -1.1587

HERE = os.path.dirname(os.path.abspath(__file__))
WT = os.path.dirname(HERE)
WORKTREES = os.path.dirname(WT)
PROC = os.path.join(WT, "data", "processed")

ENLARGED_JSON = os.path.join(PROC, "floor_enlarged.json")
ORIG_FLOOR_JSON = os.path.join(WORKTREES, "floor-measurement", "data", "processed",
                               "floor.json")
PERQ_JSON = os.path.join(WORKTREES, "per-query-tiering", "data", "processed",
                         "per_query_tiering.json")
OUT_JSON = os.path.join(PROC, "equivalence.json")


def _phi(z: float) -> float:
    return 0.5 * math.erfc(-z / math.sqrt(2))


def wilson_ci(k: int, n: int, z: float = 1.959963984540054):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    z2 = z * z
    denom = 1 + z2 / n
    centre = (p + z2 / (2 * n)) / denom
    half = (z * math.sqrt(p * (1 - p) / n + z2 / (4 * n * n))) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def two_proportion_z(k1, n1, k2, n2):
    p1, p2 = k1 / n1, k2 / n2
    p = (k1 + k2) / (n1 + n2)
    se = math.sqrt(p * (1 - p) * (1 / n1 + 1 / n2))
    z = (p1 - p2) / se if se else float("nan")
    return z, math.erfc(abs(z) / math.sqrt(2))


def fisher_exact(k1, n1, k2, n2):
    # two-sided Fisher via hypergeometric enumeration (no scipy dependency)
    from math import comb
    a, b = k1, n1 - k1
    c, d = k2, n2 - k2
    r1, r2 = a + b, c + d
    c1 = a + c
    tot = r1 + r2

    def pmf(x):
        return (comb(r1, x) * comb(r2, c1 - x)) / comb(tot, c1)
    p_obs = pmf(a)
    lo = max(0, c1 - r2)
    hi = min(c1, r1)
    return sum(pmf(x) for x in range(lo, hi + 1) if pmf(x) <= p_obs * (1 + 1e-9))


def main() -> int:
    print(f"[tost] PINNED line = {LINE} | equivalence margin DELTA = "
          f"{DELTA_PP} pp | alpha = {ALPHA}")

    # ---- PET branch (fixed, from per_query_tiering.json) -------------------- #
    perq = json.load(open(PERQ_JSON))
    pet = perq["petase_rates"]["conditional_above_given_triad"]
    k1, n1 = int(pet["k"]), int(pet["n"])           # 96 / 294
    ache = perq["ache_rates"]["conditional_above_given_triad"]
    ka, na = int(ache["k"]), int(ache["n"])         # 26 / 297

    # ---- baseline: enlarged pool if present, else original n=28 ------------- #
    if os.path.exists(ENLARGED_JSON):
        d = json.load(open(ENLARGED_JSON))
        pf = d["pooled_floor"]
        k2, n2 = int(pf["above_line"]), int(pf["triad_positive_S4"])
        base_src = f"enlarged pool (complete={d.get('complete')})"
    else:
        fj = json.load(open(ORIG_FLOOR_JSON))
        fc = fj["floor_conditional_above_given_triad"]
        k2, n2 = int(fc["k"]), int(fc["n"])
        base_src = "ORIGINAL n=28 (enlarged pool not found)"

    p1, p2 = k1 / n1, k2 / n2
    diff = p1 - p2                                   # PET - baseline
    se = math.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)

    # TOST (90% CI of diff inside +/- DELTA)
    z_upper = (diff - DELTA) / se                    # H0: diff >= DELTA
    p_upper = _phi(z_upper)                           # reject if small -> not enriched above by DELTA
    z_lower = (diff + DELTA) / se                     # H0: diff <= -DELTA
    p_lower = 1 - _phi(z_lower)                        # reject if small -> not below by DELTA
    ci90 = (diff - Z_A * se, diff + Z_A * se)
    equivalent = (p_upper < ALPHA) and (p_lower < ALPHA)

    # one-sided 95% CI upper bound on enrichment (largest plausible PET - baseline)
    upper_one_sided_95 = diff + Z_A * se

    z2, z2p = two_proportion_z(k1, n1, k2, n2)
    fish = fisher_exact(k1, n1, k2, n2)

    # verdict
    if equivalent:
        verdict = (f"EQUIVALENT within {DELTA_PP} pp: the 90% CI of the PET-minus-"
                   f"baseline difference ({100*ci90[0]:.1f}, {100*ci90[1]:.1f}) pp "
                   f"lies inside (-{DELTA_PP}, +{DELTA_PP}).")
    elif p_upper < ALPHA and ci90[1] < DELTA:
        # non-superiority holds; check whether it is below by > DELTA (depleted)
        if ci90[1] < -DELTA:
            verdict = (f"NON-SUPERIORITY established and the PET branch is DEPLETED: "
                       f"its rate is significantly below the random floor by more "
                       f"than {DELTA_PP} pp (two-sided equivalence therefore fails "
                       f"on the lower side — PET is not merely 'not enriched', it is "
                       f"below the hydrolase floor). The largest plausible enrichment "
                       f"(one-sided 95% upper bound) is {100*upper_one_sided_95:+.1f} "
                       f"pp.")
        else:
            verdict = (f"NON-SUPERIORITY established: the PET branch is not enriched "
                       f"above the random floor by more than {DELTA_PP} pp "
                       f"(one-sided p={p_upper:.2g}); largest plausible enrichment "
                       f"(one-sided 95% upper bound) {100*upper_one_sided_95:+.1f} pp. "
                       f"Strict two-sided equivalence is not claimed.")
    else:
        verdict = (f"INCONCLUSIVE at delta={DELTA_PP} pp: non-superiority not "
                   f"established (one-sided p={p_upper:.2g}); report as 'no detectable "
                   f"enrichment, equivalence inconclusive'.")

    out = {
        "pinned_line": LINE,
        "delta_pp": DELTA_PP,
        "alpha": ALPHA,
        "baseline_source": base_src,
        "pet_branch": {"k": k1, "n": n1, "rate": p1, "wilson95": list(wilson_ci(k1, n1))},
        "baseline": {"k": k2, "n": n2, "rate": p2, "wilson95": list(wilson_ci(k2, n2))},
        "ache_branch": {"k": ka, "n": na, "rate": ka / na,
                        "wilson95": list(wilson_ci(ka, na))},
        "difference_pet_minus_baseline": diff,
        "se_unpooled": se,
        "tost": {
            "z_upper": z_upper, "p_upper_not_enriched_above_delta": p_upper,
            "z_lower": z_lower, "p_lower_not_below_delta": p_lower,
            "ci90_diff": list(ci90),
            "equivalent_within_delta": equivalent,
        },
        "one_sided_95_upper_bound_on_enrichment": upper_one_sided_95,
        "legacy_tests": {
            "two_proportion_z": z2, "two_proportion_p": z2p, "fisher_p": fish,
        },
        "verdict": verdict,
    }
    os.makedirs(PROC, exist_ok=True)
    with open(OUT_JSON, "w") as fh:
        json.dump(out, fh, indent=2)

    print(f"[tost] baseline: {base_src} -> {k2}/{n2} = {100*p2:.2f}% "
          f"(Wilson {100*out['baseline']['wilson95'][0]:.1f}-"
          f"{100*out['baseline']['wilson95'][1]:.1f}%)")
    print(f"[tost] PET: {k1}/{n1} = {100*p1:.2f}%  |  diff (PET-baseline) = "
          f"{100*diff:+.2f} pp  (90% CI {100*ci90[0]:+.1f},{100*ci90[1]:+.1f})")
    print(f"[tost] one-sided 95% upper bound on enrichment = {100*upper_one_sided_95:+.2f} pp")
    print(f"[tost] TOST p_upper={p_upper:.3g} p_lower={p_lower:.3g} "
          f"equivalent={equivalent} | two-prop p={z2p:.3g} Fisher p={fish:.3g}")
    print(f"[tost] VERDICT: {verdict}")
    print(f"[tost] wrote {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
