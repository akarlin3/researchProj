"""CP2 gate: the benchmark runs, is PROVISIONAL-flagged, and tells the honest story.

These assertions are *relative and structural* (e.g. "conformal coverage gap is no
worse than raw"), never tuned targets -- the honest gate forbids fitting numbers.
"""
from __future__ import annotations

from dataclasses import replace

import pytest

from datum import run as R
from datum.baselines import BASELINES
from datum.task import CURRENT_TASK


@pytest.fixture(scope="module")
def quick_rows():
    task = replace(CURRENT_TASK, n_train=200, n_cal=400, n_test=600, n_bootstrap=100)
    include = [k for k, b in BASELINES.items() if b.estimator != "maf"]  # skip torch
    rows, meta = R.run_benchmark(task=task, include=include, seed=7, verbose=False)
    return rows, meta


def test_runs_and_emits_rows(quick_rows):
    rows, meta = quick_rows
    assert rows, "no rows produced"
    # marginal rows for every non-MAF baseline x 3 params
    marg = [r for r in rows if r["stratum"] == "all"]
    params = {r["param"] for r in marg}
    assert params == {"D", "f", "Dstar"}


def test_every_number_is_provisional(quick_rows):
    rows, _ = quick_rows
    assert all(r["provisional"] is True for r in rows)
    assert all("Fashion" in r["ruler"] for r in rows)


def test_bootstrap_ci_brackets_point(quick_rows):
    rows, _ = quick_rows
    for r in rows:
        if r["stratum"] == "all" and r.get("coverage_gap_lo") is not None:
            assert r["coverage_gap_lo"] <= r["coverage_gap"] + 1e-9
            assert r["coverage_gap"] <= r["coverage_gap_hi"] + 1e-9


def test_conformal_improves_dstar_coverage_gap(quick_rows):
    """Honest, relative check: CQR's |D* gap| is no worse than the raw segmented ref."""
    rows, _ = quick_rows
    def gap(key):
        return next(r["coverage_gap"] for r in rows
                    if r["baseline"] == key and r["param"] == "Dstar"
                    and r["stratum"] == "all")
    raw = abs(gap("reference_segmented"))
    cqr = abs(gap("reference_cqr"))
    assert cqr <= raw + 1e-6, f"CQR |gap| {cqr} worse than raw {raw}"
