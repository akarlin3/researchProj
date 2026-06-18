"""Tests for the evaluation harness (caliper.benchmark).

These are the CP1 gate: the grid runs end-to-end, every CSV number traces to a
run, and fixed seeds reproduce the table. They use the numpy-only reference
estimator and small cohorts so the whole file runs in ~1s with no torch.
"""
import csv
import math

import pytest

from caliper import benchmark as B

# small, fast grid; SNR 80 high enough to expose the D* identifiability gap.
GRID = dict(
    estimators=["reference"],
    snrs=(20.0, 80.0),
    seeds=(0, 1),
    n_cal=1500,
    n_test=2500,
    verbose=False,
)


@pytest.fixture(scope="module")
def rows():
    return B.run_grid(**GRID)


def _avg(rows, **eq):
    sel = [r for r in rows if all(r[k] == v for k, v in eq.items())]
    vals = [r["coverage"] for r in sel if r["coverage"] is not None]
    return sum(vals) / len(vals)


# --------------------------------------------------------------------------- #
# Schema / end-to-end
# --------------------------------------------------------------------------- #
def test_grid_runs_and_has_expected_shape(rows):
    # 1 est x 4 calib x 2 snr x 2 seed x 3 param x 4 strata
    assert len(rows) == 1 * 4 * 2 * 2 * 3 * 4
    assert {r["calibration"] for r in rows} == set(B.CALIBRATIONS)
    assert {r["param"] for r in rows} == {"D", "f", "Dstar"}
    assert {r["stratum"] for r in rows} == {"all", "dstar_lo", "dstar_mid", "dstar_hi"}
    assert all(set(r) == set(B.CSV_COLUMNS) for r in rows)


def test_marginal_metrics_only_on_all_stratum(rows):
    for r in rows:
        if r["stratum"] == "all":
            assert r["ece"] is not None and not math.isnan(r["ece"])
            assert r["mean_pinball"] is not None
        else:
            # tercile rows carry coverage + width but not marginal-only metrics
            assert r["ece"] is None
            assert r["mean_interval_score"] is None
            assert r["coverage"] is not None and r["mean_width"] is not None


# --------------------------------------------------------------------------- #
# The headline patterns (qualitative, robust to small-n noise)
# --------------------------------------------------------------------------- #
def test_raw_is_overconfident_and_degrades_as_snr_drops(rows):
    hi = _avg(rows, calibration="raw", param="Dstar", stratum="all", snr=80.0)
    lo = _avg(rows, calibration="raw", param="Dstar", stratum="all", snr=20.0)
    assert lo < hi          # lower SNR -> worse coverage
    assert hi < 0.75        # raw D* under-covers even at SNR 80 (nominal 0.90)


def test_conformal_restores_marginal_coverage(rows):
    for calib in ("split", "CQR", "Mondrian"):
        for snr in (20.0, 80.0):
            cov = _avg(rows, calibration=calib, param="Dstar", stratum="all", snr=snr)
            assert 0.85 <= cov <= 0.95, f"{calib}@{snr}: {cov}"


def test_split_and_cqr_coincide_for_homoscedastic_reference(rows):
    # the reference reports constant widths, so CQR's adaptive term cancels and
    # CQR reduces exactly to split-conformal residual.
    s = {(r["param"], r["snr"], r["seed"], r["stratum"]): r["coverage"]
         for r in rows if r["calibration"] == "split"}
    for r in rows:
        if r["calibration"] != "CQR":
            continue
        k = (r["param"], r["snr"], r["seed"], r["stratum"])
        assert r["coverage"] == pytest.approx(s[k], abs=1e-9)


def test_marginal_cqr_leaves_high_dstar_conditional_gap_at_high_snr(rows):
    lo = _avg(rows, calibration="CQR", param="Dstar", stratum="dstar_lo", snr=80.0)
    hi = _avg(rows, calibration="CQR", param="Dstar", stratum="dstar_hi", snr=80.0)
    # low-D* over-covers, high-D* under-covers: marginal CQR cannot close it.
    assert lo - hi > 0.05, f"expected conditional gap, got lo={lo} hi={hi}"


def test_mondrian_equalizes_coverage_at_a_width_cost(rows):
    # per-tercile coverage roughly equalized...
    covs = [_avg(rows, calibration="Mondrian", param="Dstar", stratum=s, snr=80.0)
            for s in ("dstar_lo", "dstar_mid", "dstar_hi")]
    assert max(covs) - min(covs) < 0.06
    # ...bought with high/low D* width inflation.
    def width(stratum):
        sel = [r for r in rows if r["calibration"] == "Mondrian" and r["param"] == "Dstar"
               and r["stratum"] == stratum and r["snr"] == 80.0]
        return sum(r["mean_width"] for r in sel) / len(sel)
    assert width("dstar_hi") / width("dstar_lo") > 2.0


# --------------------------------------------------------------------------- #
# Reproducibility + CSV traceability (CP1 gate)
# --------------------------------------------------------------------------- #
def test_fixed_seeds_reproduce_the_table():
    assert B.check_reproducible(**GRID) is True


def test_csv_roundtrip_matches_rows(tmp_path, rows):
    path = B.write_csv(rows, str(tmp_path / "bench.csv"))
    back = list(csv.DictReader(open(path)))
    assert len(back) == len(rows)
    assert list(back[0].keys()) == B.CSV_COLUMNS
    for r, b in zip(rows, back):
        assert b["estimator"] == r["estimator"]
        assert b["calibration"] == r["calibration"]
        assert float(b["snr"]) == r["snr"]
        # numeric round-trip at the written precision
        if r["coverage"] is not None:
            assert float(b["coverage"]) == pytest.approx(r["coverage"], abs=1e-9)
        if r["ece"] in ("", None):
            assert r["ece"] is None
