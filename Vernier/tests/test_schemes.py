"""b-scheme registry: validation, scan-time, segmented-fit support."""
from __future__ import annotations

import numpy as np
import pytest

from vernier.schemes import (
    CANDIDATE_POOL,
    EFFICIENCY_FRONTIER,
    MATCHED_SCANTIME,
    BScheme,
)


def test_validation_rejects_unsorted():
    with pytest.raises(ValueError):
        BScheme("bad", (0.0, 50.0, 20.0))


def test_validation_rejects_negative():
    with pytest.raises(ValueError):
        BScheme("bad", (-1.0, 50.0, 200.0))


def test_averages_length_must_match():
    with pytest.raises(ValueError):
        BScheme("bad", (0.0, 50.0, 200.0, 400.0), averages=(1, 1))


def test_scan_minutes_scales_with_acquisitions():
    s1 = BScheme("a", (0.0, 50.0, 200.0, 400.0))
    s2 = BScheme("b", (0.0, 50.0, 200.0, 400.0), averages=(1, 1, 2, 2))
    assert s1.n_acquisitions == 4
    assert s2.n_acquisitions == 6
    assert s2.scan_minutes() > s1.scan_minutes()


def test_matched_scantime_is_actually_matched():
    counts = {s.n_acquisitions for s in MATCHED_SCANTIME}
    minutes = {round(s.scan_minutes(), 6) for s in MATCHED_SCANTIME}
    assert len(counts) == 1, "matched-scan-time schemes must share acquisition count"
    assert len(minutes) == 1
    assert all(s.n_b == 11 for s in MATCHED_SCANTIME)


def test_all_registry_schemes_support_segmented_fit():
    for s in MATCHED_SCANTIME + CANDIDATE_POOL + EFFICIENCY_FRONTIER:
        s.require_segmented_fit()  # raises if not


def test_segmented_fit_rejects_high_b_starved_scheme():
    # only one b >= 200 -> cannot run the two-point tissue fit
    s = BScheme("starved", (0.0, 5.0, 10.0, 20.0, 35.0, 800.0))
    assert not s.supports_segmented_fit()


def test_efficiency_frontier_spans_scantime():
    nbs = sorted(s.n_b for s in EFFICIENCY_FRONTIER)
    assert nbs == [7, 11, 15, 22]
    mins = [s.scan_minutes() for s in EFFICIENCY_FRONTIER]
    assert mins == sorted(mins), "more b-values must cost more scanner-minutes"


def test_b_array_roundtrip():
    s = MATCHED_SCANTIME[0]
    assert isinstance(s.b, np.ndarray)
    assert s.b.shape[0] == s.n_b
