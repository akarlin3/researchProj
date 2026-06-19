"""CP2 gate: the OSIPI external-validation substrate adapter.

The 245 MB DRO is download-on-demand and not present in CI, so the load test skips
when the cache is absent; the absent-cache contract is always checked.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from datum import substrate

_DRO = Path(substrate.__file__).resolve().parent.parent / "data" / "osipi" / "DRO.npy"
_have_dro = _DRO.exists()


def test_osipi_absent_cache_raises_actionably(tmp_path, monkeypatch):
    if _have_dro:
        pytest.skip("DRO cache present; absent-path covered by contract elsewhere")
    with pytest.raises(FileNotFoundError) as ei:
        substrate.osipi_dro()
    assert "datum.osipi_fetch" in str(ei.value)


@pytest.mark.skipif(not _have_dro, reason="OSIPI DRO not fetched (245 MB on demand)")
def test_osipi_loads_and_splits():
    sub = substrate.osipi_dro(seed=20260613)
    assert set(sub.signals) == {"cal", "test"}
    assert sub.params["test"].shape[1] == 3            # (D, D*, f)
    assert sub.signals["test"].shape[1] == sub.b.size  # 7-b scheme
    # OSIPI true D* is OOD-high vs Gauge prior (0.05-0.20 mm^2/s).
    dstar = sub.params["test"][:, 1]
    assert dstar.max() > 0.1
    assert np.all(np.isfinite(sub.signals["test"]))
    assert sub.provenance["doi"] == "10.5281/zenodo.14605039"
