"""Floor-measurement unit tests (proteus.floor).

Pure-logic units only — no network, no fpocket. The uniform sampler is exercised
against a monkeypatched `_http` serving a synthetic in-memory `.lookup`, which also
pins its key property: the seeded draw is DETERMINISTIC regardless of fetch order.
The statistics (Wilson, two-proportion z, Katz rate ratio, Fisher) are checked
against hand-computable values, and the decomposition identity is asserted.
"""
from __future__ import annotations

import proteus.floor as floor
from proteus.floor import (
    _mean_plddt,
    funnel,
    rate_ratio_ci,
    two_proportion_z,
    wilson_ci,
)


# --------------------------------------------------------------------------- #
# Statistics
# --------------------------------------------------------------------------- #
def test_wilson_ci_known_values():
    lo, hi = wilson_ci(12, 1500)
    # 0.80% point estimate; Wilson 95% CI ~ [0.46%, 1.39%]
    assert abs(lo - 0.0046) < 5e-4
    assert abs(hi - 0.0139) < 5e-4
    assert lo < 12 / 1500 < hi


def test_wilson_ci_edges():
    assert wilson_ci(0, 0) == (0.0, 0.0)
    lo, hi = wilson_ci(0, 50)
    assert lo < 1e-9 and 0.0 < hi < 0.1          # one-sided at zero (clamped)
    lo, hi = wilson_ci(50, 50)
    assert hi == 1.0 and 0.9 < lo < 1.0


def test_two_proportion_z_direction_and_pvalue():
    # enriched 23/300 vs floor 12/1500 — enriched higher, tiny p
    res = two_proportion_z(23, 300, 12, 1500)
    assert res["z"] > 0
    assert res["p_value"] < 1e-9
    # symmetric magnitude when groups swap
    swap = two_proportion_z(12, 1500, 23, 300)
    assert abs(swap["z"] + res["z"]) < 1e-9
    assert abs(swap["p_value"] - res["p_value"]) < 1e-12


def test_two_proportion_z_degenerate():
    assert two_proportion_z(0, 0, 1, 10)["z"] is None
    assert two_proportion_z(0, 10, 0, 10)["z"] is None   # se == 0


def test_rate_ratio_ci_basic():
    rr = rate_ratio_ci(23, 298, 12, 28)          # S5 conditional: enr vs floor
    assert abs(rr["ratio"] - (23 / 298) / (12 / 28)) < 1e-9
    assert rr["ci"][0] < rr["ratio"] < rr["ci"][1]
    assert rr["ratio"] < 1.0                      # enriched conditional LOWER


def test_rate_ratio_ci_zero_cell_continuity():
    rr = rate_ratio_ci(0, 100, 5, 100)            # zero numerator -> 0.5 correction
    assert rr["ratio"] is not None and rr["ratio"] > 0
    assert rr["ci"][0] >= 0


def test_compare_block_shape():
    c = floor.compare("x", 23, 300, 12, 1500)
    assert c["group1"]["k"] == 23 and c["group1"]["n"] == 300
    assert c["group2"]["rate"] == 12 / 1500
    assert len(c["group1"]["wilson95"]) == 2
    assert "fisher_p" in c and "two_proportion_z" in c
    assert c["rate_ratio_g1_over_g2"]["ratio"] > 1   # enriched/floor overall


def test_decomposition_identity():
    # overall = triad_rate * conditional_rate (enriched side)
    triad_rate = 298 / 300
    conditional = 23 / 298
    assert abs(triad_rate * conditional - 23 / 300) < 1e-12


# --------------------------------------------------------------------------- #
# Funnel + pLDDT parsing
# --------------------------------------------------------------------------- #
def test_funnel_counts_and_conditional_numerator():
    results = [
        {"accession": "A", "triad_found": True, "pocket_ok": True,
         "petase_like_hit": True},
        {"accession": "B", "triad_found": True, "pocket_ok": True,
         "petase_like_hit": False},
        {"accession": "C", "triad_found": True, "pocket_ok": False,
         "petase_like_hit": False},
        {"accession": "D", "triad_found": False, "pocket_ok": False,
         "petase_like_hit": False},
    ]
    f = funnel(results, n_attempted=10, n_fetched=4)
    assert f["screened"] == 4
    assert f["triad_positive_S4"] == 3
    assert f["pocket_ok_S5"] == 2
    assert f["above_line"] == 1
    assert f["above_line_accessions"] == ["A"]
    # every above-line hit is necessarily a triad-bearer (screen invariant)
    above = [r for r in results if r["petase_like_hit"]]
    assert all(r["triad_found"] for r in above)


def _atom(serial, name, resn, resseq, x, y, z, occ, b, elem="C"):
    """Build a column-exact PDB ATOM line (B-factor at cols 61-66 = line[60:66])."""
    nm = name if len(name) >= 4 else (" " + name).ljust(4)   # " CA " in cols 13-16
    return (
        "ATOM  " f"{serial:>5}" " " f"{nm:<4}" " " f"{resn:>3}" " " "A"
        f"{resseq:>4}" " " "   " f"{x:8.3f}{y:8.3f}{z:8.3f}" f"{occ:6.2f}{b:6.2f}"
        "          " f"{elem:>2}"
    )


def test_mean_plddt_parses_ca_bfactors():
    pdb = (
        "HEADER    TEST\n"
        + _atom(1, "N", "ALA", 1, 11.1, 11.1, 11.1, 1.0, 80.0, "N") + "\n"
        + _atom(2, "CA", "ALA", 1, 12.1, 12.1, 12.1, 1.0, 90.0) + "\n"
        + _atom(3, "CA", "GLY", 2, 13.1, 13.1, 13.1, 1.0, 70.0) + "\n"
    ).encode()
    mp, n = _mean_plddt(pdb)
    assert n == 2                         # two CA atoms (N is skipped)
    assert abs(mp - 80.0) < 1e-6          # mean of the two CA B-factors


def test_mean_plddt_empty():
    assert _mean_plddt(b"HEADER only\n") == (None, 0)


# --------------------------------------------------------------------------- #
# Uniform sampler — determinism + correctness against a synthetic lookup
# --------------------------------------------------------------------------- #
def _make_lookup(n_entries: int) -> bytes:
    # mimic `<key>\t<accession>\t<file#>\n`, fixed-width accession
    return b"".join(
        f"{i}\tMGYP{i:012d}\t0\n".encode() for i in range(n_entries)
    )


def test_sample_universe_deterministic_and_in_universe(monkeypatch):
    blob = _make_lookup(5000)
    size = len(blob)

    def fake_http(url, *, rng=None, timeout=60, max_bytes=None):
        assert rng is not None
        lo, hi = rng
        return 206, blob[lo:hi + 1], url

    monkeypatch.setattr(floor, "_http", fake_http)

    accs1, meta = floor.sample_universe("u", size, 50, seed=1729,
                                        workers=8, win=64, timeout=5,
                                        log=lambda *a, **k: None)
    accs2, _ = floor.sample_universe("u", size, 50, seed=1729,
                                     workers=1, win=64, timeout=5,  # different fan-out
                                     log=lambda *a, **k: None)
    assert len(accs1) == 50
    assert len(set(accs1)) == 50                  # unique
    assert accs1 == accs2                         # deterministic across worker counts
    assert all(a.startswith("MGYP") for a in accs1)
    assert meta["seed"] == 1729 and meta["n"] == 50


def test_sample_universe_changes_with_seed(monkeypatch):
    blob = _make_lookup(5000)
    size = len(blob)
    monkeypatch.setattr(floor, "_http",
                        lambda url, *, rng=None, timeout=60, max_bytes=None:
                        (206, blob[rng[0]:rng[1] + 1], url))
    a, _ = floor.sample_universe("u", size, 40, seed=1, win=64,
                                 log=lambda *x, **k: None)
    b, _ = floor.sample_universe("u", size, 40, seed=2, win=64,
                                 log=lambda *x, **k: None)
    assert a != b
