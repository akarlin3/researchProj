"""S5 — Cleft filter: fpocket metrics + two real discriminators, control-anchored.

fpocket on its own collapses to "is there a druggable pocket" — every serine
hydrolase (CalB, AChE, lipases) passes that. S5 adds the features that actually
separate a PET hydrolase from a non-PET serine hydrolase that shares the fold and
triad:

  1. EXPOSURE / lid-absence (primary). PET hydrolases present the catalytic serine on
     an open surface groove; lipases bury it under a mobile lid and esterases/AChE
     down a deep gorge. The robust proxy is the catalytic-Ser *peripherality* relative
     to the protein CA centroid (large = surface, small = core). This is preferred over
     raw OG SASA, which on crystal controls is dominated by the open-vs-closed lid state
     (an open-lid lipase crystal shows a deceptively exposed Ser). Raw OG SASA is still
     computed and reported as a diagnostic column.
     Peripherality is expressed via `peripherality_mode` (config): `percentile` (rank of
     the Ser among its own residues; PRODUCTION default), `rg_norm` (OG->centroid over
     Rg), or `absolute` (raw Angstrom; size-dependent, comparison only). The first two
     are size-invariant so the metric transfers across the size-heterogeneous dark
     proteome; `absolute` cannot. See _peripherality_modes.
  2. AROMATIC subsites. Count Trp/Tyr/Phe within a shell of the Ser OG (PET binds via
     aromatic stacking, e.g. IsPETase Trp185). Down-weighted: raw counts are noisy and
     aromatic-gorge esterases inflate them; on folded models the rotamers are noisy too.

For each triad-positive model S5 selects the catalytic pocket (fpocket pocket whose
centre is nearest the catalytic Ser OG), reads the A-E base metrics, computes the two
discriminators, then z-scores every metric against the POSITIVE controls and combines
them with config weights/orientation into one composite. Nothing is hardcoded; weights,
orientation, shells and the positive set all come from config (s5_cleft_filter).

Local usage, from the repo root:
    PYTHONPATH=src python -m proteus.s5_cleft_filter --pdb structures/6EQE.pdb --ser 160
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

from proteus.utils import (
    DEFAULT_CONFIG,
    REPO,
    euclidean,
    load_structure,
    per_atom_sasa,
    protein_atoms,
    residue_iter,
)

_BACKBONE = {"N", "CA", "C", "O", "OXT"}


# --------------------------------------------------------------------------- #
# fpocket
# --------------------------------------------------------------------------- #
def _parse_info(info_path: str) -> dict:
    """Parse fpocket <name>_info.txt into {pocket_id: {field: float}}."""
    pockets, cur = {}, None
    with open(info_path) as fh:
        for line in fh:
            m = re.match(r"^Pocket (\d+)\s*:", line)
            if m:
                cur = int(m.group(1))
                pockets[cur] = {}
                continue
            if cur is not None and ":" in line:
                key, _, val = line.strip().partition(":")
                try:
                    pockets[cur][key.strip()] = float(val.strip())
                except ValueError:
                    pass
    return pockets


def _pocket_center(vert_pqr: str):
    """Centroid of a pocket's alpha-sphere centres (Voronoi vertices)."""
    coords = []
    with open(vert_pqr) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return np.mean(coords, axis=0) if coords else None


def run_fpocket(pdb_path: str) -> dict:
    """Run fpocket on a copy of `pdb_path` in a temp dir; return
    {pocket_id: {**info_fields, 'center': np.ndarray}}. Empty dict if it produced none."""
    if shutil.which("fpocket") is None:
        raise RuntimeError("fpocket not on PATH")
    work = tempfile.mkdtemp(prefix="proteus_fpocket_")
    try:
        base = os.path.basename(pdb_path)
        stem = os.path.splitext(base)[0]
        shutil.copy(pdb_path, os.path.join(work, base))
        subprocess.run(["fpocket", "-f", base], cwd=work,
                       capture_output=True, text=True)
        out_dir = os.path.join(work, f"{stem}_out")
        info_path = os.path.join(out_dir, f"{stem}_info.txt")
        if not os.path.exists(info_path):
            return {}
        info = _parse_info(info_path)
        pockets = {}
        for pid, fields in info.items():
            vert = os.path.join(out_dir, "pockets", f"pocket{pid}_vert.pqr")
            center = _pocket_center(vert) if os.path.exists(vert) else None
            if center is None:
                continue
            pockets[pid] = {**fields, "center": center}
        return pockets
    finally:
        shutil.rmtree(work, ignore_errors=True)


# --------------------------------------------------------------------------- #
# Geometry helpers for the discriminators
# --------------------------------------------------------------------------- #
def _catalytic_og(protein, ser_res_id: int):
    mask = (protein.res_id == ser_res_id) & (protein.atom_name == "OG")
    if not mask.any():
        return None
    return np.asarray(protein.coord[mask][0], dtype=float)


PERIPHERALITY_MODES = ("percentile", "rg_norm", "absolute")


def _peripherality_modes(protein, og, ser_res_id: int) -> dict:
    """Catalytic-Ser peripherality (exposure proxy) in all three forms, against ONE
    shared CA centroid. Deterministic — pure geometry, no seed.

    - absolute:   ||OG - centroid||, Angstrom. Size-DEPENDENT (kept for comparison).
    - rg_norm:    ||OG - centroid|| / Rg, where Rg = sqrt(mean_i ||CA_i - centroid||^2).
                  Numerator and denominator both scale with size -> size-invariant.
    - percentile: fraction of CA atoms whose centroid-distance is < the catalytic-Ser
                  CA centroid-distance (0 = dead centre, 1 = most peripheral). Compares
                  the Ser to its OWN protein's residues -> size- and shape-invariant.

    Under a uniform coordinate scaling (all coords x k): absolute scales by k, while
    rg_norm and percentile are unchanged (ratio / rank order are scale-free)."""
    ca = protein.coord[protein.atom_name == "CA"]
    centroid = ca.mean(axis=0)
    d_og = euclidean(og, centroid)
    d_ca = np.linalg.norm(ca - centroid, axis=1)
    rg = float(np.sqrt(np.mean(d_ca ** 2)))

    ser_ca_mask = (protein.res_id == ser_res_id) & (protein.atom_name == "CA")
    if ser_ca_mask.any():
        d_ser_ca = euclidean(np.asarray(protein.coord[ser_ca_mask][0], dtype=float),
                             centroid)
        percentile = float(np.mean(d_ca < d_ser_ca))
    else:
        percentile = float("nan")

    return {
        "absolute": d_og,
        "rg_norm": (d_og / rg) if rg > 0 else float("nan"),
        "percentile": percentile,
    }


def _peripherality(protein, og, ser_res_id: int, mode: str) -> float:
    """Catalytic-Ser peripherality in the requested `mode` (see _peripherality_modes)."""
    modes = _peripherality_modes(protein, og, ser_res_id)
    if mode not in modes:
        raise ValueError(f"unknown peripherality_mode {mode!r}; "
                         f"expected one of {PERIPHERALITY_MODES}")
    return modes[mode]


def _og_sasa(protein, ser_res_id: int, probe: float) -> float:
    sasa = per_atom_sasa(protein, probe_radius=probe)
    mask = (protein.res_id == ser_res_id) & (protein.atom_name == "OG")
    return float(np.nansum(sasa[mask])) if mask.any() else float("nan")


def _aromatic_count(protein, og, resnames, shell_min, shell_max) -> int:
    arom = set(resnames)
    count = 0
    for _chain, _rid, rname, sub in residue_iter(protein):
        if rname not in arom:
            continue
        side = sub.coord[~np.isin(sub.atom_name, list(_BACKBONE))]
        if len(side) == 0:
            continue
        d = min(euclidean(og, x) for x in side)
        if shell_min <= d <= shell_max:
            count += 1
    return count


# --------------------------------------------------------------------------- #
# Per-model cleft analysis
# --------------------------------------------------------------------------- #
def analyze_cleft(pdb_path: str, ser_res_id: int, cfg: dict) -> dict:
    """Select the catalytic pocket nearest the Ser OG and compute the cleft metrics.

    Returns {pocket_id, dist_og_pocket, n_pockets, metrics:{...}, raw_og_sasa, ...}.
    pocket_id is None if no pocket lies within catalytic_pocket_max_dist of the OG."""
    s5 = cfg["s5_cleft_filter"]
    max_dist = float(s5["catalytic_pocket_max_dist"])
    protein = protein_atoms(load_structure(pdb_path))
    og = _catalytic_og(protein, ser_res_id)
    mode = s5.get("peripherality_mode", "percentile")
    result = {"model": os.path.basename(pdb_path), "ser_res_id": ser_res_id,
              "pocket_id": None, "dist_og_pocket": None, "n_pockets": 0,
              "metrics": {}, "raw_og_sasa": None,
              "peripherality_mode": mode, "peripherality": None}
    if og is None:
        result["error"] = f"no Ser OG at residue {ser_res_id}"
        return result

    pockets = run_fpocket(pdb_path)
    result["n_pockets"] = len(pockets)

    best_pid, best_d = None, None
    for pid, p in pockets.items():
        d = euclidean(og, p["center"])
        if best_d is None or d < best_d:
            best_pid, best_d = pid, d

    if s5.get("exposure_metric", "peripherality") == "peripherality":
        periph = _peripherality_modes(protein, og, ser_res_id)
        result["peripherality"] = {k: round(v, 5) for k, v in periph.items()}
        exposure = periph[mode]
    else:
        exposure = _og_sasa(protein, ser_res_id, float(s5["sasa_probe_radius"]))
    aromatics = _aromatic_count(protein, og, s5["aromatic_resnames"],
                                float(s5["aromatic_shell_min"]),
                                float(s5["aromatic_shell_max"]))
    result["raw_og_sasa"] = _og_sasa(protein, ser_res_id, float(s5["sasa_probe_radius"]))

    if best_pid is None or best_d > max_dist:
        result["dist_og_pocket"] = (round(best_d, 3) if best_d is not None else None)
        result["error"] = "no pocket within catalytic_pocket_max_dist"
        return result

    p = pockets[best_pid]
    solvent_access = p.get("Mean alp. sph. solvent access", float("nan"))
    result["pocket_id"] = best_pid
    result["dist_og_pocket"] = round(best_d, 3)
    result["metrics"] = {
        "exposure": round(float(exposure), 3),
        "aromatics": float(aromatics),
        "volume": round(p.get("Volume", float("nan")), 3),
        "hydrophobicity": round(p.get("Hydrophobicity score", float("nan")), 3),
        "polarity": round(p.get("Polarity score", float("nan")), 3),
        "druggability": round(p.get("Druggability Score", float("nan")), 3),
        # depth = pocket buriedness = 1 - mean alpha-sphere solvent access
        "depth": round(1.0 - solvent_access, 3),
    }
    return result


# --------------------------------------------------------------------------- #
# Control-anchored scoring
# --------------------------------------------------------------------------- #
def _robust_scale(pos_vals: np.ndarray, all_vals: np.ndarray) -> float:
    """Scale for z-scoring: positive-control std, floored by half the overall spread
    (positive std is unreliable with ~2 positives) and a tiny epsilon."""
    pos_std = float(np.std(pos_vals)) if len(pos_vals) > 1 else 0.0
    all_std = float(np.std(all_vals)) if len(all_vals) > 1 else 0.0
    return max(pos_std, 0.5 * all_std, 1e-6)


def composite_from_anchor(metrics: dict, anchor: dict, cfg: dict) -> dict:
    """Score ONE cleft's metrics against a PRECOMPUTED anchor (center/scale per
    metric), using the configured weights + orientation. Returns {z, composite}.

    This is the scoring kernel shared by control calibration (anchor built from the
    positives) and dark-candidate screening (same positive-control anchor applied
    to a folded model). `metrics.exposure` must already be in the anchor's
    peripherality mode."""
    s5 = cfg["s5_cleft_filter"]
    weights = s5["weights"]
    orient = s5["orientation"]
    zr, composite = {}, 0.0
    for metric, w in weights.items():
        z = (metrics[metric] - anchor[metric]["center"]) / anchor[metric]["scale"]
        sign = orient.get(metric, 0)
        contrib = (sign * z) if sign != 0 else (-abs(z))  # 0 => penalise deviation
        zr[metric] = round(z, 3)
        composite += w * contrib
    return {"z": zr, "composite": round(composite, 4)}


def score_controls(metrics_by_id: dict, positive_ids, cfg: dict) -> dict:
    """z-score each metric against the positive controls and combine into a composite.

    metrics_by_id: {control_id: {metric: value}}  (only models with a catalytic pocket)
    Returns {control_id: {z:{metric:z}, composite: float}} plus the per-metric
    center/scale used (under key '_anchor')."""
    s5 = cfg["s5_cleft_filter"]
    weights = s5["weights"]
    ids = list(metrics_by_id)
    pos_ids = [i for i in positive_ids if i in metrics_by_id]

    anchor = {}
    for metric in weights:
        all_vals = np.array([metrics_by_id[i][metric] for i in ids], dtype=float)
        pos_vals = np.array([metrics_by_id[i][metric] for i in pos_ids], dtype=float)
        center = float(np.mean(pos_vals)) if len(pos_vals) else float(np.mean(all_vals))
        scale = _robust_scale(pos_vals, all_vals)
        anchor[metric] = {"center": round(center, 4), "scale": round(scale, 4)}

    out = {"_anchor": anchor}
    for cid in ids:
        out[cid] = composite_from_anchor(metrics_by_id[cid], anchor, cfg)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb", required=True, help="PDB to analyze")
    ap.add_argument("--ser", required=True, type=int, help="catalytic Ser residue id")
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    args = ap.parse_args(argv)

    from proteus.utils import load_config  # noqa: PLC0415
    cfg = load_config(args.config)
    if not os.path.exists(args.pdb):
        print(f"PDB not found: {args.pdb}", file=sys.stderr)
        return 2
    r = analyze_cleft(args.pdb, args.ser, cfg)
    if r["pocket_id"] is None:
        print(f"[s5] {r['model']}: no catalytic pocket "
              f"({r.get('error','')}; nearest={r['dist_og_pocket']}A)")
        return 1
    m = r["metrics"]
    print(f"[s5] {r['model']} Ser{args.ser}: pocket {r['pocket_id']} "
          f"@{r['dist_og_pocket']}A  exposure={m['exposure']} "
          f"(mode={r['peripherality_mode']}) aromatics={m['aromatics']:.0f} "
          f"vol={m['volume']} drug={m['druggability']} depth={m['depth']} "
          f"(raw OG SASA={r['raw_og_sasa']:.2f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
