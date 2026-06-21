"""CP4 gate: the reframe pipeline runs and has the right structure / orderings.

Fast smoke (small cohort, short chain) -- the committed numbers come from the full
seeded run (`scripts/build_handoff.py`); here we only lock that the table is well-formed
and that the qualitative orderings the hand-off relies on hold.
"""
from __future__ import annotations

from gnomon import reframe


def test_reframe_table_structure_and_orderings(tmp_path):
    out = reframe.run_reframe(n_noise=30, out_dir=tmp_path, verbose=False)
    t = out["conditional_coverage"]

    # both affected estimators present, both SD conventions present
    assert set(t["Laplace_SD"]) == {"honest", "floored"}
    assert set(t["MCMC_SD"]) == {"honest", "floored"}
    assert "MCMC_quantile_recommended" in t

    for est in ("Laplace_SD", "MCMC_SD"):
        for conv in ("honest", "floored"):
            cells = t[est][conv]
            assert set(cells) >= {"low_Dstar", "mid_Dstar", "high_Dstar", "pooled"}
            for c in cells.values():
                assert "coverage" in c and "ci" in c and len(c["ci"]) == 2

        honest, floored = t[est]["honest"], t[est]["floored"]
        # the overconfident floor can only LOWER coverage (never raise it)
        assert floored["pooled"]["coverage"] <= honest["pooled"]["coverage"] + 1e-9
        # the failure concentrates in high D* (the wall): high <= low tercile
        assert honest["high_Dstar"]["coverage"] <= honest["low_Dstar"]["coverage"] + 1e-9

    assert (tmp_path / "conditional_coverage.json").exists()
    assert out["recommended_convention"] == "honest"
