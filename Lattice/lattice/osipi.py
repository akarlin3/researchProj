"""Optional OSIPI TF2.4 reference-DRO integration via download-on-demand.

Lattice ships only *synthetic* data in-tree. The OSIPI TF2.4 IVIM reference DRO
is external CC-BY-4.0 data; this module fetches it on demand (never committed)
and writes a provenance manifest, mirroring Gauge's posture. Importing this
module pulls in nothing external; the network is touched only when you call
:func:`fetch`.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional

import numpy as np

__all__ = ["OSIPI_SOURCE", "ExternalDRO", "provenance_record", "load_dro", "fetch"]

# External source coordinates (CC-BY-4.0 data; Apache-2.0 code upstream).
OSIPI_SOURCE = {
    "name": "OSIPI TF2.4 IVIM reference DRO",
    "zenodo_doi": "10.5281/zenodo.14605039",
    "zenodo_url": "https://zenodo.org/records/14605039",
    "data_license": "CC-BY-4.0",
    "code_license": "Apache-2.0",
    "citation": (
        "OSIPI Taskforce 2.4. OSIPI TF2.4 IVIM-MRI Code Collection. "
        "Zenodo; doi:10.5281/zenodo.14605039."
    ),
}


@dataclass
class ExternalDRO:
    """Ground truth + native signals loaded from an external reference DRO."""

    params: np.ndarray   # (N, 3) ordered (D, Dstar, f)
    signals: np.ndarray  # (N, n_b)
    bvalues: np.ndarray   # (n_b,)
    source: dict

    def __len__(self) -> int:
        return self.params.shape[0]


def provenance_record(local_path: str | Path, sha256: Optional[str] = None, fetched: Optional[str] = None) -> dict:
    """Build a provenance manifest dict for a fetched external artifact.

    ``fetched`` (an ISO timestamp) is passed in by the caller -- this module
    never reads the wall clock, so manifests stay reproducible/auditable.
    """
    return {
        "source": OSIPI_SOURCE,
        "local_path": str(local_path),
        "sha256": sha256,
        "fetched": fetched,
        "note": "External CC-BY-4.0 data; not redistributed in-tree.",
    }


def load_dro(path: str | Path) -> ExternalDRO:
    """Load an external DRO previously written by :func:`fetch`.

    Expects an ``.npz`` with ``params (N,3) [D, Dstar, f]``, ``signals (N,n_b)``,
    and ``bvalues (n_b,)``.
    """
    path = Path(path)
    with np.load(path) as z:
        return ExternalDRO(
            params=z["params"],
            signals=z["signals"],
            bvalues=z["bvalues"],
            source=dict(OSIPI_SOURCE),
        )


def fetch(dest_dir: str | Path, fetched: Optional[str] = None) -> Path:  # pragma: no cover - network
    """Download the external OSIPI DRO into ``dest_dir`` and write provenance.

    Network-dependent; not exercised in unit tests. Raises with a clear message
    if the optional ``requests`` dependency or network is unavailable. The
    downloaded artifact is *not* committed to the repository.
    """
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    try:
        import requests  # optional; declared in the 'external-data' extra
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError(
            "fetching the external OSIPI DRO needs the 'external-data' extra "
            "(pip install lattice[external-data])"
        ) from exc

    raise NotImplementedError(
        "Wire the concrete Zenodo asset URL from OSIPI_SOURCE['zenodo_url'] here; "
        "write provenance via provenance_record(...) alongside the artifact. "
        "Left unimplemented so no network call is hard-coded into the package."
    )
