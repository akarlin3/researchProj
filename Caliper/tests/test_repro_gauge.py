"""Tests for caliper.repro_gauge -- the Gauge reproduction module.

These pin the three qualitative Gauge findings the module exists to reproduce on
synthetic data, using only caliper.conformal + caliper.metrics:

1. raw (over-confident) intervals badly under-cover, and marginal CQR restores
   pooled / marginal coverage to near-nominal;
2. conditional coverage is NOT delivered -- the high-D* tercile stays
   under-covered relative to the low-D* tercile (the identifiability wall);
3. group-conditional (Mondrian) CQR equalizes per-tercile coverage only by
   inflating the high-D* interval width.

All data is in-repo synthetic; no clinical numbers; fixed seeds make it
deterministic.
"""
import numpy as np

from caliper import repro_gauge as RG


def test_raw_is_overconfident_and_marginal_cqr_restores():
    r = RG.reproduce()
    # raw reference estimator badly under-covers, worst on D*
    assert r.raw_marginal["Dstar"] < 0.5
    assert r.raw_marginal["f"] < 0.7
    # marginal CQR restores every parameter to within 0.03 of nominal
    assert r.marginal_restored
    for p in ("D", "f", "Dstar"):
        assert abs(r.cqr_marginal[p] - r.nominal) <= 0.03
    # pooled D* coverage under marginal CQR is near-nominal
    assert abs(r.dstar_pooled_cqr - r.nominal) <= 0.03


def test_high_dstar_tercile_stays_undercovered_under_marginal():
    r = RG.reproduce()
    cqr = r.per_tercile["marginal-CQR"]
    # the identifiability wall: high-D* under-covers, low-D* over-covers
    assert cqr[2].coverage < r.nominal              # high-D* below nominal
    assert cqr[0].coverage > cqr[2].coverage        # low-D* over- vs high under
    assert r.high_dstar_undercovered_marginal


def test_mondrian_restores_per_tercile_only_by_inflating_width():
    r = RG.reproduce()
    mond = r.per_tercile["Mondrian-CQR"]
    # per-tercile coverage equalized to near-nominal...
    for s in mond:
        assert abs(mond[s].coverage - r.nominal) <= 0.03
    # ...but only by inflating the high-D* interval width (constant under CQR)
    assert r.mondrian_width_ratio > 1.5
    cqr = r.per_tercile["marginal-CQR"]
    cqr_ratio = cqr[2].mean_width / cqr[0].mean_width
    assert np.isclose(cqr_ratio, 1.0, atol=1e-6)    # CQR holds one width
    assert r.mondrian_restores_by_inflation


def test_phenomenon_holds_and_is_deterministic():
    r1 = RG.reproduce()
    r2 = RG.reproduce()
    assert r1.phenomenon_holds
    # fixed seeds -> byte-identical key numbers across runs
    assert r1.dstar_pooled_cqr == r2.dstar_pooled_cqr
    assert r1.mondrian_width_ratio == r2.mondrian_width_ratio
    assert r1.cqr_marginal == r2.cqr_marginal


def test_format_is_honest_about_pre_publication():
    text = RG.reproduce().format().lower()
    # the human-readable summary must not claim publication
    assert "pre-publication" in text
    assert "published" not in text
