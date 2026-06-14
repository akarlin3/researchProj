"""Per-query-tiering unit tests (proteus.per_query).

Pure-logic units only — no network, no fpocket. Covers: the statistics (Wilson,
two-proportion z, Katz rate ratio, Fisher) against hand-computable values; the
best-match partition of a synthetic result.m8 (including the global-best vs
within-anchor distinction that is the whole point of the module); the seeded
within-anchor selection (top-N + DETERMINISTIC tiebreak); the funnel; the pLDDT
parse; and the gradient-verdict classifier (above / flat / inverted).
"""
from __future__ import annotations

import os
import tempfile

import proteus.per_query as pq
from proteus.per_query import (
    branch_sizes,
    build_partition_records,
    build_verdict,
    fisher_p,
    funnel_of,
    katz_rate_ratio,
    partition_m8,
    select_branch,
    two_proportion_z,
    wilson_ci,
)

ANCHORS = {"PETASE": ["6EQE", "4EB0", "8B4U", "4WFI", "4CG1"],
           "ACHE": ["1EA5"],
           "OTHER_NEG": ["1TCA_A", "1CRL_A", "1EVQ"]}


# --------------------------------------------------------------------------- #
# Statistics
# --------------------------------------------------------------------------- #
def test_wilson_ci_known_values():
    lo, hi = wilson_ci(12, 28)                    # the floor conditional, 42.86%
    assert abs(lo - 0.2651) < 5e-4               # published Wilson lower bound
    assert abs(hi - 0.6093) < 5e-4               # published Wilson upper bound
    assert lo < 12 / 28 < hi


def test_wilson_ci_edges():
    assert wilson_ci(0, 0) == [float("nan"), float("nan")] or \
        wilson_ci(0, 0)[0] != wilson_ci(0, 0)[0]  # nan
    lo, hi = wilson_ci(0, 50)
    assert lo < 1e-9 and 0.0 < hi < 0.1
    lo, hi = wilson_ci(50, 50)
    assert hi == 1.0 and 0.9 < lo < 1.0


def test_two_proportion_z_direction_and_symmetry():
    z, p = two_proportion_z(96, 294, 12, 28)     # PETASE vs floor conditional
    assert z < 0                                  # PETASE point estimate below floor
    z2, p2 = two_proportion_z(12, 28, 96, 294)
    assert abs(z2 + z) < 1e-9                     # sign flips on swap
    assert abs(p2 - p) < 1e-12                    # p invariant


def test_two_proportion_z_strong_separation_tiny_p():
    z, p = two_proportion_z(96, 294, 26, 297)    # PETASE vs ACHE
    assert z > 5 and p < 1e-9


def test_katz_rate_ratio_basic_and_below_one():
    rr = katz_rate_ratio(96, 294, 12, 28)        # PETASE/floor
    assert abs(rr["rr"] - (96 / 294) / (12 / 28)) < 1e-9
    assert rr["ci"][0] < rr["rr"] < rr["ci"][1]
    assert rr["rr"] < 1.0                          # flat/below floor


def test_katz_rate_ratio_above_one():
    rr = katz_rate_ratio(96, 294, 26, 297)       # PETASE/ACHE
    assert rr["rr"] > 3 and rr["ci"][0] > 1       # CI excludes 1 -> separated


def test_fisher_p_matches_known():
    # PETASE vs ACHE conditional: highly significant
    p = fisher_p(96, 294 - 96, 26, 297 - 26)
    assert p < 1e-9
    # PETASE vs floor: not significant
    p2 = fisher_p(96, 294 - 96, 12, 28 - 12)
    assert p2 > 0.05


# --------------------------------------------------------------------------- #
# CP1 — partition (best-match), including the within-anchor vs global distinction
# --------------------------------------------------------------------------- #
def _m8_line(q, acc, bits):
    # 8-col custom foldseek: query, target(.pdb.gz), evalue, bits, fident, aln, ql, tl
    return f"{q}\t{acc}.pdb.gz\t1e-30\t{bits}\t0.5\t200\t200\t210\n"


def _write_m8(rows, gz=False):
    suffix = ".m8.gz" if gz else ".m8"
    fd, path = tempfile.mkstemp(suffix=suffix)
    os.close(fd)
    if gz:
        import gzip
        with gzip.open(path, "wt") as fh:
            fh.writelines(rows)
    else:
        with open(path, "w") as fh:
            fh.writelines(rows)
    return path


def test_partition_best_match_assigns_global_best_anchor():
    # Target T1 matches a PETASE query (900) and AChE (1000) — AChE wins (best-match),
    # even though T1 is a genuine PETase-neighbour. T2 best-matches PETASE.
    rows = [
        _m8_line("4CG1", "T1", 900),
        _m8_line("1EA5", "T1", 1000),
        _m8_line("6EQE", "T2", 800),
        _m8_line("1EA5", "T2", 700),
        _m8_line("1EVQ", "T3", 500),
    ]
    path = _write_m8(rows)
    try:
        part = partition_m8(path, ANCHORS)
    finally:
        os.unlink(path)
    recs = {r["accession"]: r for r in build_partition_records(part["per_acc"])}
    assert part["n_rows"] == 5
    assert recs["T1"]["anchor_class"] == "ACHE"      # global best wins
    assert recs["T1"]["best_bits"] == 1000.0
    # ...but the PETASE bits are retained for "bits-to-nearest-PETase-query"
    assert recs["T1"]["petase_bits"] == 900.0
    assert recs["T1"]["petase_query"] == "4CG1"
    assert recs["T2"]["anchor_class"] == "PETASE"
    assert recs["T3"]["anchor_class"] == "OTHER_NEG"
    assert branch_sizes(list(recs.values())) == {"ACHE": 1, "PETASE": 1, "OTHER_NEG": 1}


def test_partition_gzip_and_unknown_query_aborts():
    ok = _write_m8([_m8_line("6EQE", "T1", 500)], gz=True)
    try:
        part = partition_m8(ok, ANCHORS)
        assert part["n_rows"] == 1
    finally:
        os.unlink(ok)
    bad = _write_m8([_m8_line("ZZZZ", "T1", 500)])
    try:
        import pytest
        with pytest.raises(RuntimeError):
            partition_m8(bad, ANCHORS)
    finally:
        os.unlink(bad)


# --------------------------------------------------------------------------- #
# CP2 — within-anchor selection (seeded, deterministic, top-N)
# --------------------------------------------------------------------------- #
def _rec(acc, cls, bits):
    return {"accession": acc, "anchor_class": cls, "best_query": "Q",
            "best_bits": float(bits), "petase_bits": bits, "petase_query": "Q",
            "ache_bits": "", "other_bits": ""}


def test_select_branch_top_n_by_bits():
    recs = [_rec(f"P{i}", "PETASE", 1000 - i) for i in range(10)]
    recs += [_rec("A0", "ACHE", 9999)]            # other branch ignored
    sel = select_branch(recs, "PETASE", 3, seed=1729)
    assert [r["accession"] for r in sel] == ["P0", "P1", "P2"]   # highest bits first


def test_select_branch_all_when_smaller_than_n():
    recs = [_rec(f"P{i}", "PETASE", 500 + i) for i in range(4)]
    sel = select_branch(recs, "PETASE", 300, seed=1729)
    assert len(sel) == 4


def test_select_branch_tiebreak_is_deterministic_given_seed():
    # all equal bits -> ordering decided entirely by the seeded tiebreak
    recs = [_rec(f"P{i}", "PETASE", 500) for i in range(20)]
    a = [r["accession"] for r in select_branch(recs, "PETASE", 5, seed=1729)]
    b = [r["accession"] for r in select_branch(list(reversed(recs)), "PETASE", 5,
                                               seed=1729)]
    assert a == b                                 # input order does not matter
    c = [r["accession"] for r in select_branch(recs, "PETASE", 5, seed=2026)]
    assert a != c                                 # a different seed reorders ties


# --------------------------------------------------------------------------- #
# Funnel + pLDDT
# --------------------------------------------------------------------------- #
def test_funnel_of_counts_and_invariant():
    records = [
        {"accession": "A", "fetched": True, "triad_found": True, "pocket_ok": True,
         "above_threshold": True},
        {"accession": "B", "fetched": True, "triad_found": True, "pocket_ok": True,
         "above_threshold": False},
        {"accession": "C", "fetched": True, "triad_found": True, "pocket_ok": False,
         "above_threshold": None},
        {"accession": "D", "fetched": True, "triad_found": False, "pocket_ok": False,
         "above_threshold": None},
        {"accession": "E", "fetched": False, "triad_found": False, "pocket_ok": False,
         "above_threshold": None},
    ]
    f = funnel_of(records, -1.1587)
    assert f["attempted"] == 5 and f["fetched"] == 4
    assert f["triad_positive_S4"] == 3
    assert f["pocket_ok_S5"] == 2
    assert f["above_line"] == 1
    assert f["above_line_accessions"] == ["A"]


def test_mean_plddt_parses_ca_on_0_1_scale():
    # ESM Atlas stores pLDDT 0-1 in the B-factor column
    lines = (
        "ATOM      1  N   ALA A   1      11.100  11.100  11.100  1.00  0.80           N\n"
        "ATOM      2  CA  ALA A   1      12.100  12.100  12.100  1.00  0.90           C\n"
        "ATOM      3  CA  GLY A   2      13.100  13.100  13.100  1.00  0.70           C\n"
    )
    fd, path = tempfile.mkstemp(suffix=".pdb")
    os.close(fd)
    with open(path, "w") as fh:
        fh.write(lines)
    try:
        mp = pq._mean_plddt(path)
    finally:
        os.unlink(path)
    assert abs(mp - 0.80) < 1e-6                  # mean of the two CA B-factors


# --------------------------------------------------------------------------- #
# CP4 — gradient verdict classifier
# --------------------------------------------------------------------------- #
def _rates(k, n):
    return {"triad_rate": {"k": n, "n": n, "rate": 1.0, "wilson95": wilson_ci(n, n)},
            "conditional_above_given_triad": {"k": k, "n": n, "rate": k / n,
                                              "wilson95": wilson_ci(k, n)}}


def test_verdict_flat_when_indistinguishable_from_floor():
    fl = {"above_line": 12, "triad_positive_S4": 28, "screened": 1500}
    pr = _rates(96, 294)                          # 32.65%
    ar = _rates(26, 297)                          # 8.75%
    t_pf = pq.two_arm_test(96, 294, 12, 28, "PETASE", "floor")
    t_pa = pq.two_arm_test(96, 294, 26, 297, "PETASE", "ACHE")
    verdict, tldr = build_verdict(t_pf, pr, fl, t_pa, ar)
    assert "flat" in verdict.lower()
    assert "FLAT" in tldr
    assert "not blind" in verdict                 # PETASE-vs-ACHE nuance woven in


def test_verdict_present_when_substantially_above_floor():
    fl = {"above_line": 12, "triad_positive_S4": 28, "screened": 1500}
    pr = _rates(250, 294)                          # ~85%, well above floor
    ar = _rates(26, 297)
    t_pf = pq.two_arm_test(250, 294, 12, 28, "PETASE", "floor")
    t_pa = pq.two_arm_test(250, 294, 26, 297, "PETASE", "ACHE")
    verdict, tldr = build_verdict(t_pf, pr, fl, t_pa, ar)
    assert "gradient is present" in verdict.lower()
    assert "PRESENT" in tldr


def test_verdict_inverted_when_significantly_below_floor():
    # large floor n so a low PETASE rate is significantly below it
    fl = {"above_line": 600, "triad_positive_S4": 1000, "screened": 1500}  # 60%
    pr = _rates(30, 294)                            # ~10%, far below 60%
    ar = _rates(26, 297)
    t_pf = pq.two_arm_test(30, 294, 600, 1000, "PETASE", "floor")
    t_pa = pq.two_arm_test(30, 294, 26, 297, "PETASE", "ACHE")
    verdict, tldr = build_verdict(t_pf, pr, fl, t_pa, ar)
    assert "inverted" in verdict.lower()
    assert "INVERTED" in tldr
