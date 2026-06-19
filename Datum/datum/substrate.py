"""Data substrates for Datum -- read-only adapters over Gauge (and, later, Lattice).

Datum benchmarks methods *on data*. The natural substrate, **Lattice**, is not
built yet, so the primary substrate is **Gauge's synthetic cohort**; an OSIPI
digital reference object (DRO) is wired as an external-validation substrate via
Gauge's fetch/provenance scripts. Both are reused, not reinvented. A ``lattice()``
swap-in point is reserved so that when Lattice exists it replaces ``gauge_cohort``
as ``primary`` (the build-order dependency recorded in ``datum.manifest``).

Guardrail: synthetic-only. The OSIPI DRO is itself synthetic; it is downloaded on
demand and git-ignored (only its provenance manifest is tracked). No clinical /
in-vivo / MSK data is ever materialised in this tree.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from datum import _paths

_paths.ensure_deps()

from gauge.cohort import generate_cohort, DEFAULT_SEED, DEFAULT_SNR_GRID  # noqa: E402
from gauge.forward import DEFAULT_B_VALUES  # noqa: E402

from datum.manifest import SUBSTRATE  # noqa: E402


@dataclass
class Substrate:
    """A scored-on data substrate: signals + ground truth, keyed by split."""
    name: str
    b: np.ndarray                 # (n_b,) b-values
    signals: dict                 # split -> (N, n_b) noisy signals
    params: dict                  # split -> (N, 3) ground truth (D, D*, f)
    snr: dict                     # split -> (N,) per-sample SNR (may be None)
    provenance: dict              # how this substrate was produced/sourced

    @property
    def sizes(self) -> dict:
        return {k: v.shape[0] for k, v in self.signals.items()}


def gauge_cohort(n_train: int = 3000, n_cal: int = 2000, n_test: int = 3000,
                 snr_grid=DEFAULT_SNR_GRID, b=DEFAULT_B_VALUES,
                 seed: int = DEFAULT_SEED) -> Substrate:
    """Primary substrate: Gauge's seeded synthetic IVIM cohort (read-only).

    Wraps ``gauge.cohort.generate_cohort`` unchanged. Defaults reproduce Gauge's
    production cohort; pass smaller sizes for smoke tests.
    """
    c = generate_cohort(n_train, n_cal, n_test, snr_grid=snr_grid, b=b, seed=seed)
    return Substrate(
        name=SUBSTRATE["primary"]["name"],
        b=np.asarray(c.b, dtype=float),
        signals=c.signals, params=c.params, snr=c.snr,
        provenance={
            "kind": "synthetic (Gauge cohort)",
            "entrypoint": "gauge.cohort.generate_cohort",
            "seed": seed,
            "snr_grid": tuple(snr_grid),
            "n_b": int(np.asarray(c.b).shape[0]),
            "commit": SUBSTRATE["primary"]["commit"],
        },
    )


def osipi_dro(data_dir: str | Path | None = None) -> Substrate:
    """External-validation substrate: the OSIPI TF2.4 IVIM DRO (synthetic).

    The DRO is fetched on demand by ``Gauge/scripts/fetch_osipi.py`` and cached
    under a git-ignored ``data/`` directory; only the provenance manifest is
    tracked. This adapter does NOT download: it reads an already-fetched cache and
    its provenance, and raises an actionable error if the cache is absent so the
    benchmark never silently runs on missing data.
    """
    root = _paths.find_monorepo_root()
    prov_path = (root / SUBSTRATE["external_validation"]["provenance"]) if root else None
    if prov_path is None or not prov_path.exists():
        raise FileNotFoundError(
            "OSIPI DRO provenance not found. Fetch it first with Gauge's script:\n"
            f"    python {SUBSTRATE['external_validation']['fetch_script']}\n"
            f"(DOI {SUBSTRATE['external_validation']['doi']}). The DRO is synthetic "
            "and download-on-demand; only provenance is committed."
        )
    provenance = json.loads(prov_path.read_text())
    raise NotImplementedError(
        "OSIPI DRO loading is wired at CP1 (provenance located at "
        f"{prov_path}) and is populated by the CP2 benchmark build once the DRO "
        "cache is fetched. Provenance keys available: "
        f"{sorted(provenance)[:6]}..."
    )


def lattice(*args, **kwargs) -> Substrate:
    """Reserved swap-in point for the Lattice substrate (NOT BUILT yet)."""
    raise NotImplementedError(
        "Lattice substrate is not built yet. " + SUBSTRATE["planned"]["note"]
        + " Until then, use gauge_cohort() (primary) or osipi_dro() (validation)."
    )


# Registry the task spec references by name.
SUBSTRATES = {
    "gauge_cohort": gauge_cohort,
    "osipi_dro": osipi_dro,
    "lattice": lattice,
}
