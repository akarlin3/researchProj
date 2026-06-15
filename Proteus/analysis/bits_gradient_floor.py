"""Recompute the PET-branch bits-stratified tail gradient against the ENLARGED
random floor (revision Fix A.1 follow-through).

The original per-query tiering compared each bits tier to the n=28 floor (42.86%).
The enlarged floor is 105/209 = 50.24%, so the vs-floor significances and the
crossing point must be recomputed. This script reads the committed bits_gradient
(per_query_tiering.json) and the enlarged floor (floor_enlarged.json) and re-tests
each tier; nothing is hand-typed.

Reproduce:  PYTHONPATH=src python analysis/bits_gradient_floor.py
"""
from __future__ import annotations

import json
import math
import os
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
WT = os.path.dirname(HERE)
PROC = os.path.join(WT, "data", "processed")
PERQ = os.path.join(PROC, "inputs_snapshot", "per_query_tiering.json")
ENL = os.path.join(PROC, "floor_enlarged.json")
OUT = os.path.join(PROC, "bits_gradient_enlarged.json")


def fisher_exact(k1, n1, k2, n2):
    a, b = k1, n1 - k1
    c, d = k2, n2 - k2
    r1, r2, c1, tot = a + b, c + d, a + c, n1 + n2

    def pmf(x):
        return comb(r1, x) * comb(r2, c1 - x) / comb(tot, c1)
    po = pmf(a)
    lo, hi = max(0, c1 - r2), min(c1, r1)
    return sum(pmf(x) for x in range(lo, hi + 1) if pmf(x) <= po * (1 + 1e-9))


def main() -> int:
    pq = json.load(open(PERQ))
    enl = json.load(open(ENL))["pooled_floor"]
    fk, fn = int(enl["above_line"]), int(enl["triad_positive_S4"])
    frate = fk / fn
    print(f"[bits] enlarged floor = {fk}/{fn} = {100*frate:.2f}%")
    tiers = []
    for b in pq["bits_gradient"]:
        k, n = int(b["k"]), int(b["n"])
        rate = k / n
        p = fisher_exact(k, n, fk, fn)
        tiers.append({"label": b["label"], "k": k, "n": n, "rate": rate,
                      "min_bits": b.get("min_bits"),
                      "rr_vs_enlarged_floor": rate / frate,
                      "fisher_p_vs_enlarged_floor": p,
                      "above_floor": rate > frate and p < 0.05,
                      "below_floor": rate < frate and p < 0.05})
        print(f"[bits] {b['label']:<12} {k:>3}/{n:<3} {100*rate:5.1f}%  "
              f"rr={rate/frate:.2f}x  Fisher p={p:.4g}")
    out = {"enlarged_floor": {"k": fk, "n": fn, "rate": frate}, "tiers": tiers}
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"[bits] wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
