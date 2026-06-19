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

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from datum import _paths

_paths.ensure_deps()

from gauge.cohort import DEFAULT_SEED, DEFAULT_SNR_GRID, generate_cohort  # noqa: E402
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


_DRO_REL = "data/osipi/DRO.npy"      # git-ignored cache, under the Datum package root
_DATUM_ROOT = Path(__file__).resolve().parent.parent


def osipi_dro(n_cal: int | None = None, seed: int = 20260613) -> Substrate:
    """External-validation substrate: the OSIPI TF2.4 IVIM DRO (synthetic).

    The DRO (``Utilities/DRO.npy`` from Zenodo record 14605039 -- 5000 voxels, a
    fixed sparse 7-b acquisition pre-noised at ~SNR 80, ground truth (D, D*, f) we
    did NOT generate) is downloaded on demand and cached under the git-ignored
    ``Datum/data/osipi/`` directory; only its provenance is committed. The fetch
    reuses Gauge's pinned DOI/URL/MD5 (see ``datum.osipi_fetch``) -- Gauge itself
    is never modified. This loader reads the cache and splits it into cal/test so
    the conformal baselines have a calibration set; it raises an actionable error
    if the cache is absent so the benchmark never silently runs on missing data.

    Ground truth is returned in Gauge convention ``(D, D*, f)`` (physical mm^2/s),
    matching :func:`gauge_cohort`, so the same downstream conversion applies. NB:
    OSIPI's true D* (~0.05-0.20) is shifted HIGHER than Gauge's prior -- a genuine
    out-of-distribution external test for the analytic, b-flexible baselines.
    """
    dro_path = _DATUM_ROOT / _DRO_REL
    if not dro_path.exists():
        raise FileNotFoundError(
            f"OSIPI DRO cache not found at {dro_path}. Fetch it on demand:\n"
            "    python -m datum.osipi_fetch\n"
            f"(DOI {SUBSTRATE['external_validation']['doi']}; synthetic, 245 MB "
            "download, md5-verified, cached git-ignored).")
    dro = np.load(dro_path, allow_pickle=True)
    D = np.array([float(e["D"]) for e in dro])
    Dstar = np.array([float(e["Dp"]) for e in dro])
    f = np.array([float(e["f"]) for e in dro])
    sig = np.array([np.asarray(e["signals"], dtype=float) for e in dro])
    b = np.asarray(dro[0]["bvals"], dtype=float)
    params = np.stack([D, Dstar, f], axis=1)        # Gauge convention (D, D*, f)

    n = sig.shape[0]
    n_cal = n // 2 if n_cal is None else int(n_cal)
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    cal_idx, test_idx = perm[:n_cal], perm[n_cal:]
    return Substrate(
        name=SUBSTRATE["external_validation"]["name"],
        b=b,
        signals={"cal": sig[cal_idx], "test": sig[test_idx]},
        params={"cal": params[cal_idx], "test": params[test_idx]},
        snr={"cal": None, "test": None},
        provenance={
            "kind": "synthetic external DRO (OSIPI TF2.4)",
            "doi": SUBSTRATE["external_validation"]["doi"],
            "n_voxels": int(n), "n_b": int(b.size), "split_seed": seed,
            "dstar_range_mm2_s": [float(Dstar.min()), float(Dstar.max())],
            "note": "OOD-high D* vs Gauge prior; pre-noised ~SNR 80; 7-b sparse.",
        },
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
