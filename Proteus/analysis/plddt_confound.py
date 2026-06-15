"""pLDDT confound test (revision Fix A.2): is the divergent-tail exposed-cleft
depletion an ESMFold prediction-error artifact, or biology?

The manuscript flags but never runs this: "prediction error at low identity cannot
be excluded entirely and merits an explicit pLDDT comparison between the divergent
band and the baseline." This script runs it.

Three groups, all from CACHED data (no re-fetch):
  (a) DIVERGENT band  — the sub-25%-identity screened tail (sensitivity-research,
        2,100 hits). Global + active-site-local pLDDT from the cached ESMFold V0
        PDBs (B-factor column, 0-1). <20% subset called out.
  (b) BASELINE        — random-floor triad-bearers (enlarged pool if available,
        else the original n=28). Global pLDDT from the floor CSV; active-site-local
        where the triad+ PDB was retained.
  (c) NEAR-HOMOLOG    — the PET-branch above-line hits (seqid-analysis), which sit
        predominantly at >=30% identity. Global pLDDT from the seqid CSV.

Compare distributions (Mann-Whitney U + KS + Cliff's delta). Decision:
  * divergent pLDDT materially lower -> pLDDT-MATCHED re-screen: subsample the
    divergent triad-bearers to the baseline pLDDT distribution and recompute the
    exposed-cleft rate; report whether the depletion PERSISTS after matching.
  * divergent pLDDT comparable -> prediction quality does not explain the depletion
    (strengthens the claim).

Reproduce, from the repo root:
    PYTHONPATH=src python analysis/plddt_confound.py
Prints the pinned line (-1.1587) and the matching RNG seed used for any subsampling.
"""
from __future__ import annotations

import csv
import glob
import json
import os

import numpy as np
from scipy.stats import ks_2samp, mannwhitneyu

LINE = -1.1587            # pinned decision line (for provenance; not re-screened here)
MATCH_SEED = 1729         # RNG seed for pLDDT-matched subsampling (config random_seed)

HERE = os.path.dirname(os.path.abspath(__file__))
WT = os.path.dirname(HERE)
WORKTREES = os.path.dirname(WT)
PROC = os.path.join(WT, "data", "processed")

SENS = os.path.join(WORKTREES, "sensitivity-research")
SENS_CSV = os.path.join(SENS, "data", "processed", "sensitivity_per_hit.csv")
SENS_CACHE = os.path.join(SENS, "data", "interim", "sensitivity_cache")
SEQID_CSV = os.path.join(WORKTREES, "seqid-analysis", "data", "processed",
                         "seqid_per_hit.csv")
PERQ_CACHE = os.path.join(WORKTREES, "per-query-tiering", "data", "interim",
                          "per_query_cache")
ORIG_FLOOR_CSV = os.path.join(WORKTREES, "floor-measurement", "data", "processed",
                              "floor.csv")
ENLARGED_CSV = os.path.join(PROC, "floor_enlarged_pooled.csv")
ENLARGED_HITS = os.path.join(WT, "structures", "floor_enlarged_hits")

OUT_JSON = os.path.join(PROC, "plddt_confound.json")


def _tb(v) -> bool:
    return str(v).strip().lower() in ("true", "1")


def ca_bfactors(pdb_path: str) -> dict[int, float]:
    """resnum -> CA pLDDT (B-factor, 0-1) for one ESMFold V0 PDB."""
    out: dict[int, float] = {}
    with open(pdb_path, "rb") as fh:
        for line in fh:
            if line.startswith(b"ATOM") and line[12:16].strip() == b"CA":
                try:
                    res = int(line[22:26])
                    out[res] = float(line[60:66])
                except ValueError:
                    pass
    return out


def local_plddt(bf: dict[int, float], triad: list[int], win: int = 1) -> float | None:
    """Mean pLDDT over triad residues +/- win (active-site-local proxy; fpocket
    pocket-residue lists are not persisted, so we use the catalytic triad and an
    immediate neighbourhood)."""
    res = set()
    for t in triad:
        if t is None:
            continue
        for d in range(-win, win + 1):
            res.add(t + d)
    vals = [bf[r] for r in res if r in bf]
    return float(np.mean(vals)) if vals else None


def find_pdb(cache_dirs: list[str], acc: str) -> str | None:
    for d in cache_dirs:
        p = os.path.join(d, f"{acc}.pdb")
        if os.path.exists(p):
            return p
    return None


def load_csv(path: str) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh))


def _intick(v):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return None


# --------------------------------------------------------------------------- #
# Group assembly
# --------------------------------------------------------------------------- #
def group_divergent():
    """(a) sub-25% divergent tail; pLDDT from cached PDBs."""
    rows = load_csv(SENS_CSV)
    recs = []
    for r in rows:
        acc = r["accession"]
        pdb = find_pdb([SENS_CACHE], acc)
        if not pdb:
            continue
        bf = ca_bfactors(pdb)
        if not bf:
            continue
        g = float(np.mean(list(bf.values())))
        triad = [_intick(r.get("catalytic_ser")), _intick(r.get("his")),
                 _intick(r.get("acid"))]
        loc = local_plddt(bf, triad)
        seqid = None
        try:
            seqid = float(r["seqid_nearest"]) if r.get("seqid_nearest") else None
        except ValueError:
            seqid = None
        recs.append({"acc": acc, "global": g, "local": loc, "seqid": seqid,
                     "triad": _tb(r.get("triad_found")),
                     "above": _tb(r.get("above_threshold")),
                     "low_cov": _tb(r.get("low_coverage"))})
    return recs


def group_baseline():
    """(b) random-floor triad-bearers; enlarged pool if present else original 28."""
    path = ENLARGED_CSV if os.path.exists(ENLARGED_CSV) else ORIG_FLOOR_CSV
    rows = load_csv(path)
    recs = []
    for r in rows:
        if not _tb(r["triad_found"]):
            continue
        try:
            g = float(r["mean_plddt"]) if r.get("mean_plddt") else None
        except ValueError:
            g = None
        if g is None:
            continue
        acc = r["accession"]
        pdb = find_pdb([os.path.join(ENLARGED_HITS, d) for d in
                        (os.listdir(ENLARGED_HITS) if os.path.isdir(ENLARGED_HITS)
                         else [])], acc)
        loc = None
        if pdb:
            bf = ca_bfactors(pdb)
            triad = [_intick(r.get("catalytic_ser")), _intick(r.get("his")),
                     _intick(r.get("acid"))]
            loc = local_plddt(bf, triad)
        recs.append({"acc": acc, "global": g, "local": loc,
                     "above": _tb(r.get("petase_like_hit"))})
    return recs, os.path.basename(path)


def group_nearhomolog():
    """(c) PET-branch above-line hits (near-homologs); pLDDT from seqid CSV."""
    rows = load_csv(SEQID_CSV)
    recs = []
    for r in rows:
        if not _tb(r.get("above_line")):
            continue
        try:
            g = float(r["mean_plddt"]) if r.get("mean_plddt") else None
        except ValueError:
            g = None
        if g is None:
            continue
        seqid = None
        try:
            seqid = float(r["seqid_nearest"]) if r.get("seqid_nearest") else None
        except ValueError:
            pass
        recs.append({"acc": r.get("mgyp"), "global": g, "seqid": seqid})
    return recs


# --------------------------------------------------------------------------- #
# Stats
# --------------------------------------------------------------------------- #
def cliffs_delta(a: list[float], b: list[float]) -> float:
    """Cliff's delta via the Mann-Whitney U rank-biserial identity:
    delta = 2U/(n*m) - 1, where U is for sample `a`. Positive => a > b."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    n, m = len(a), len(b)
    if n == 0 or m == 0:
        return float("nan")
    U, _ = mannwhitneyu(a, b, alternative="two-sided")
    return float(2.0 * U / (n * m) - 1.0)


def compare(name: str, a: list[float], b: list[float]) -> dict:
    a = [x for x in a if x is not None]
    b = [x for x in b if x is not None]
    mwu_p = ks_p = None
    if a and b:
        _, mwu_p = mannwhitneyu(a, b, alternative="two-sided")
        _, ks_p = ks_2samp(a, b)
    return {
        "comparison": name,
        "n_a": len(a), "n_b": len(b),
        "median_a": float(np.median(a)) if a else None,
        "median_b": float(np.median(b)) if b else None,
        "mean_a": float(np.mean(a)) if a else None,
        "mean_b": float(np.mean(b)) if b else None,
        "mannwhitney_p": float(mwu_p) if mwu_p is not None else None,
        "ks_p": float(ks_p) if ks_p is not None else None,
        "cliffs_delta_a_vs_b": cliffs_delta(a, b) if a and b else None,
    }


def matched_rescreen(divergent: list[dict], baseline_global: list[float]) -> dict:
    """pLDDT-matched re-screen: subsample the divergent triad-bearers so their
    global-pLDDT histogram matches the baseline's, then recompute the above-line
    (exposed-cleft) rate. Reports whether the depletion persists after matching."""
    rng = np.random.default_rng(MATCH_SEED)
    # restrict to credible-coverage triad-bearers (the manuscript's screened set)
    dv = [d for d in divergent if d["triad"] and not d["low_cov"]
          and d["global"] is not None]
    if not dv or not baseline_global:
        return {"note": "insufficient data for matched re-screen"}
    bins = np.linspace(0.7, 1.0, 7)  # HQ set is pLDDT>0.7
    base_hist, _ = np.histogram(baseline_global, bins=bins)
    base_frac = base_hist / base_hist.sum()
    dv_glob = np.array([d["global"] for d in dv])
    dv_bin = np.digitize(dv_glob, bins) - 1
    dv_bin = np.clip(dv_bin, 0, len(bins) - 2)
    # target total = as large as possible while honouring baseline bin fractions
    cap = []
    for bi in range(len(bins) - 1):
        avail = int(np.sum(dv_bin == bi))
        if base_frac[bi] > 0:
            cap.append(avail / base_frac[bi])
    total = int(min(cap)) if cap else 0
    picks = []
    for bi in range(len(bins) - 1):
        k = int(round(base_frac[bi] * total))
        idx = np.where(dv_bin == bi)[0]
        if k > 0 and len(idx) > 0:
            picks.extend(rng.choice(idx, size=min(k, len(idx)), replace=False))
    picks = np.array(picks, dtype=int)
    matched = [dv[i] for i in picks]
    raw_rate = float(np.mean([d["above"] for d in dv])) if dv else None
    matched_rate = float(np.mean([d["above"] for d in matched])) if matched else None
    return {
        "raw_divergent_triad_credible_n": len(dv),
        "raw_above_line_rate": raw_rate,
        "matched_n": len(matched),
        "matched_above_line_rate": matched_rate,
        "matched_global_plddt_mean": float(np.mean([d["global"] for d in matched]))
        if matched else None,
        "baseline_global_plddt_mean": float(np.mean(baseline_global)),
        "note": "above-line rate among credible-coverage divergent triad-bearers, "
                "raw vs after matching the baseline global-pLDDT distribution.",
    }


def main() -> int:
    print(f"[plddt] PINNED line = {LINE} | matched-subsample seed = {MATCH_SEED}")
    a = group_divergent()
    b, b_src = group_baseline()
    c = group_nearhomolog()
    print(f"[plddt] groups: divergent={len(a)} (pdb-backed), "
          f"baseline={len(b)} (from {b_src}), near-homolog={len(c)}")

    a_g = [r["global"] for r in a]
    a_sub20_g = [r["global"] for r in a if r.get("seqid") is not None and r["seqid"] < 20]
    b_g = [r["global"] for r in b]
    c_g = [r["global"] for r in c]
    a_loc = [r["local"] for r in a if r.get("local") is not None]
    b_loc = [r["local"] for r in b if r.get("local") is not None]

    summary = {
        "pinned_line": LINE,
        "match_seed": MATCH_SEED,
        "baseline_source": b_src,
        "n": {"divergent": len(a), "divergent_sub20": len(a_sub20_g),
              "baseline": len(b), "near_homolog": len(c)},
        "global_plddt": {
            "divergent_median": float(np.median(a_g)) if a_g else None,
            "divergent_sub20_median": float(np.median(a_sub20_g)) if a_sub20_g else None,
            "baseline_median": float(np.median(b_g)) if b_g else None,
            "near_homolog_median": float(np.median(c_g)) if c_g else None,
        },
        "tests_global": [
            compare("divergent vs baseline", a_g, b_g),
            compare("divergent vs near-homolog", a_g, c_g),
            compare("baseline vs near-homolog", b_g, c_g),
            compare("divergent_sub20 vs baseline", a_sub20_g, b_g),
        ],
        "tests_local": [
            compare("divergent vs baseline (active-site-local)", a_loc, b_loc),
        ] if b_loc else [],
        "matched_rescreen": matched_rescreen(a, b_g),
    }

    # decision
    div_med = summary["global_plddt"]["divergent_median"]
    base_med = summary["global_plddt"]["baseline_median"]
    dvb = summary["tests_global"][0]
    delta = dvb["cliffs_delta_a_vs_b"]
    materially_lower = (div_med is not None and base_med is not None
                        and (base_med - div_med) >= 0.02
                        and dvb["mannwhitney_p"] is not None
                        and dvb["mannwhitney_p"] < 0.05
                        and delta is not None and delta <= -0.15)
    summary["decision"] = {
        "divergent_materially_lower_than_baseline": bool(materially_lower),
        "interpretation": (
            "Divergent-band pLDDT is materially lower; ran a pLDDT-matched re-screen "
            "to test whether the exposed-cleft depletion persists after matching."
            if materially_lower else
            "Divergent-band pLDDT is comparable to the baseline (and near-homolog) "
            "groups; ESMFold prediction quality does not explain the exposed-cleft "
            "depletion in the divergent tail."),
    }

    # arrays for the figure
    summary["_arrays"] = {
        "divergent_global": a_g,
        "divergent_sub20_global": a_sub20_g,
        "baseline_global": b_g,
        "near_homolog_global": c_g,
    }
    os.makedirs(PROC, exist_ok=True)
    with open(OUT_JSON, "w") as fh:
        json.dump(summary, fh, indent=2)

    # console report
    print(f"[plddt] global pLDDT medians: divergent={div_med:.4f} "
          f"baseline={base_med:.4f} near-homolog="
          f"{summary['global_plddt']['near_homolog_median']:.4f}")
    print(f"[plddt] divergent vs baseline: MWU p={dvb['mannwhitney_p']:.3g} "
          f"KS p={dvb['ks_p']:.3g} Cliff's d={delta:+.3f}")
    print(f"[plddt] DECISION: divergent materially lower? {materially_lower}")
    mr = summary["matched_rescreen"]
    if "matched_above_line_rate" in mr and mr["matched_above_line_rate"] is not None:
        print(f"[plddt] matched re-screen: raw above-line "
              f"{mr['raw_above_line_rate']:.3f} -> matched "
              f"{mr['matched_above_line_rate']:.3f} (n={mr['matched_n']})")
    print(f"[plddt] wrote {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
