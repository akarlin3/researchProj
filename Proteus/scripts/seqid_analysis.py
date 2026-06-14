#!/usr/bin/env python
"""Figure 1 — sequence-identity reach analysis for the PETASE branch.

Recomputes, for every screened PETASE-branch triad-bearer, the % sequence
identity (with alignment coverage) to its structural anchor query and to the
nearest PETase query overall, using a local Smith-Waterman alignment
(BLOSUM62, gap -10/-1) in biotite. Then re-stratifies the above-line|triad
gradient by sequence identity (replacing Foldseek bits), correlates bits vs
identity, and triages the 96 above-line hits into near-homolog vs structure-first.

Deterministic; no network (hit + query sequences are read from local PDBs).
"""
import json, os, sys, glob, csv, math
import numpy as np
import biotite.structure.io.pdb as pdb
import biotite.structure as struc
import biotite.sequence as bseq
import biotite.sequence.align as balign
from scipy import stats

WT = "/Users/averykarlin/projProteus/.claude/worktrees/per-query-tiering"
STRUCT = "/Users/averykarlin/projProteus/structures"
OUT_DIR = sys.argv[1]               # worktree data/processed
TMP = sys.argv[2]                   # job tmp for the summary json

COV_MIN = 0.50                      # query-coverage floor for a credible homolog call
QUERIES = {"6EQE": "IsPETase", "4EB0": "LCC_WT", "8B4U": "PET46",
           "4WFI": "Cut190", "4CG1": "TfCut2"}
MATRIX = balign.SubstitutionMatrix.std_protein_matrix()   # BLOSUM62
GAP = (-10, -1)

def seq_from_pdb(path, chain_first=True):
    arr = pdb.PDBFile.read(path).get_structure(model=1)
    aa = arr[struc.filter_amino_acids(arr) & (arr.atom_name == "CA")]
    if chain_first and len(aa) and len(set(aa.chain_id)) > 1:
        aa = aa[aa.chain_id == sorted(set(aa.chain_id))[0]]
    out = []
    for r3 in aa.res_name:
        try:
            out.append(bseq.ProteinSequence.convert_letter_3to1(r3))
        except Exception:
            out.append("X")          # unknown -> placeholder (never counts as identity)
    return "".join(out)

def to_protseq(s):
    # biotite ProteinSequence rejects X? It accepts ambiguous incl 'X'. Replace anything odd.
    clean = "".join(c if c in bseq.ProteinSequence.alphabet else "X" for c in s)
    return bseq.ProteinSequence(clean)

def local_identity(q_seq, h_seq):
    """Smith-Waterman local align; return pident (n_ident/aln_len), cov_q, cov_h, aln_len, n_ident."""
    aln = balign.align_optimal(q_seq, h_seq, MATRIX, gap_penalty=GAP,
                               terminal_penalty=False, local=True)[0]
    tr = aln.trace                      # (L,2); -1 == gap
    L = tr.shape[0]
    if L == 0:
        return dict(pident=0.0, cov_q=0.0, cov_h=0.0, aln_len=0, n_ident=0)
    q_arr = np.frombuffer(bytes(str(q_seq), "ascii"), dtype=np.uint8)
    h_arr = np.frombuffer(bytes(str(h_seq), "ascii"), dtype=np.uint8)
    qi, hi = tr[:, 0], tr[:, 1]
    both = (qi != -1) & (hi != -1)
    n_ident = int(np.sum((qi != -1) & (hi != -1) &
                         (q_arr[np.where(qi != -1, qi, 0)] == h_arr[np.where(hi != -1, hi, 0)]) & both))
    n_q = int(np.sum(qi != -1))
    n_h = int(np.sum(hi != -1))
    return dict(pident=n_ident / L, cov_q=n_q / len(q_seq), cov_h=n_h / len(h_seq),
                aln_len=L, n_ident=n_ident)

# ---- load query sequences ----
qseqs = {acc: to_protseq(seq_from_pdb(f"{STRUCT}/{acc}.pdb")) for acc in QUERIES}
qlen = {acc: len(str(qseqs[acc])) for acc in QUERIES}
print("[queries]", {f"{a}({QUERIES[a]})": qlen[a] for a in QUERIES})

# ---- load Foldseek fident (structure-derived seq identity) for cross-check ----
import gzip
fident = {}
with gzip.open(f"{WT}/data/interim/result.m8.gz", "rt") as fh:
    for line in fh:                     # query, target, evalue, bits, fident, alnlen, qlen, tlen
        p = line.rstrip("\n").split("\t")
        fident[(p[0], p[1].split(".")[0])] = float(p[4])

# ---- load screened PETASE records ----
recs = json.load(open(f"{WT}/data/interim/per_query_screen_PETASE.json"))
cache = f"{WT}/data/interim/per_query_cache"

rows = []
for i, r in enumerate(recs):
    acc = r["accession"]
    hpath = f"{cache}/{acc}.pdb"
    hseq = to_protseq(seq_from_pdb(hpath))
    hlen = len(str(hseq))
    anchor = r["best_query"]
    # identity to every PETase query
    per_q = {}
    for acc_q in QUERIES:
        per_q[acc_q] = local_identity(qseqs[acc_q], hseq)
    # anchor identity (structural best-match query, by bits)
    a = per_q[anchor]
    # nearest PETase query OVERALL by identity, among queries meeting coverage floor
    credible = {k: v for k, v in per_q.items() if v["cov_q"] >= COV_MIN}
    if credible:
        nq = max(credible, key=lambda k: credible[k]["pident"])
        low_cov = False
    else:
        nq = max(per_q, key=lambda k: per_q[k]["pident"])   # best available, but flagged
        low_cov = True
    nn = per_q[nq]
    rows.append(dict(
        mgyp=acc, anchor=anchor, anchor_name=QUERIES[anchor],
        best_bits=r["best_bits"],
        composite=(round(r["composite"], 4) if r.get("composite") is not None else None),
        above_line=bool(r["above_threshold"]), triad=bool(r["triad_found"]),
        mean_plddt=r["mean_plddt"], hit_len=hlen,
        seqid_anchor=round(100 * a["pident"], 2), cov_anchor=round(a["cov_q"], 3),
        nearest_query=nq, nearest_name=QUERIES[nq],
        seqid_nearest=round(100 * nn["pident"], 2), cov_nearest=round(nn["cov_q"], 3),
        cov_nearest_hit=round(nn["cov_h"], 3), aln_len_nearest=nn["aln_len"],
        low_coverage=low_cov,
        # Foldseek's own structure-derived fident (independent cross-check)
        fident_anchor=(round(100 * fident[(anchor, acc)], 2) if (anchor, acc) in fident else None),
        fident_nearest=(round(100 * fident[(nq, acc)], 2) if (nq, acc) in fident else None),
        # full per-query identity for audit
        **{f"id_{q}": round(100 * per_q[q]["pident"], 2) for q in QUERIES},
        **{f"cov_{q}": round(per_q[q]["cov_q"], 3) for q in QUERIES},
    ))
    if (i + 1) % 50 == 0:
        print(f"  ... {i+1}/{len(recs)} aligned")

# ---- write per-hit table ----
os.makedirs(OUT_DIR, exist_ok=True)
csv_path = f"{OUT_DIR}/seqid_per_hit.csv"
cols = ["mgyp", "anchor", "anchor_name", "best_bits", "composite", "above_line",
        "triad", "mean_plddt", "hit_len", "seqid_anchor", "cov_anchor",
        "nearest_query", "nearest_name", "seqid_nearest", "cov_nearest",
        "cov_nearest_hit", "aln_len_nearest", "low_coverage",
        "fident_anchor", "fident_nearest"] + \
       [f"id_{q}" for q in QUERIES] + [f"cov_{q}" for q in QUERIES]
with open(csv_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    for row in rows:
        w.writerow({k: row[k] for k in cols})
print("[wrote]", csv_path, len(rows), "rows")

# ============ analysis ============
triad = [r for r in rows if r["triad"]]
above = [r for r in rows if r["above_line"]]
assert len(triad) == 294 and len(above) == 96, (len(triad), len(above))

FLOOR = 12 / 28   # 42.857%  (floor.json conditional above-line|triad)

def wilson(k, n, z=1.96):
    if n == 0: return (0.0, 0.0)
    p = k / n; d = 1 + z*z/n
    c = (p + z*z/(2*n)) / d
    h = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / d
    return (max(0, c-h), min(1, c+h))

def fisher_vs_floor(k, n):
    # 2x2: [[k, n-k],[12, 28-12]]
    _, p = stats.fisher_exact([[k, n-k], [12, 28-12]])
    return p

# ---- CP2: stratify above-line|triad by nearest-PETase seq identity ----
BINS = [("<20%", 0, 20), ("20-30%", 20, 30), ("30-40%", 30, 40),
        ("40-60%", 40, 60), (">60%", 60, 1e9)]
strat = []
for label, lo, hi in BINS:
    grp = [r for r in triad if lo <= r["seqid_nearest"] < hi]
    n = len(grp); k = sum(r["above_line"] for r in grp)
    lo_ci, hi_ci = wilson(k, n)
    strat.append(dict(bin=label, n_triad=n, k_above=k,
                      rate=(k/n if n else None),
                      wilson95=[round(lo_ci, 4), round(hi_ci, 4)],
                      fisher_p_vs_floor=(fisher_vs_floor(k, n) if n else None),
                      rr_vs_floor=((k/n)/FLOOR if n else None)))

# ---- bits <-> seqid correlation on the 294 triad-bearers ----
bits = np.array([r["best_bits"] for r in triad])
sid = np.array([r["seqid_nearest"] for r in triad])
pear = stats.pearsonr(bits, sid)
spear = stats.spearmanr(bits, sid)
# where does the bits crossover (~1090 bits, the per-query top-100 cut) map in seq-id?
near_1090 = [r["seqid_nearest"] for r in triad if 1080 <= r["best_bits"] <= 1100]
sid_at_1090 = float(np.median(near_1090)) if near_1090 else None

# ---- CP3: triage the 96 above-line ----
def tier(r):
    if r["low_coverage"]:           # no credible full-domain alignment to ANY PETase query
        return "structure_first"     # sequence search does not reach a full homolog
    if r["seqid_nearest"] >= 30: return "near_homolog"
    if r["seqid_nearest"] >= 25: return "twilight"
    return "structure_first"
for r in above:
    r["_tier"] = tier(r)
triage = {t: [r for r in above if r["_tier"] == t]
          for t in ("near_homolog", "twilight", "structure_first")}
sf = sorted(triage["structure_first"], key=lambda r: -r["composite"])

# ---- dual-metric robustness: triage the 96 above-line under Foldseek fident too ----
def tier_fident(r):
    fn = r["fident_nearest"]
    if fn is None: return "structure_first"
    if fn >= 30: return "near_homolog"
    if fn >= 25: return "twilight"
    return "structure_first"
sf_fident = sorted([r for r in above if tier_fident(r) == "structure_first"],
                   key=lambda r: -r["composite"])
triage_fident_counts = {t: sum(1 for r in above if tier_fident(r) == t)
                        for t in ("near_homolog", "twilight", "structure_first")}
# Foldseek-vs-biotite agreement on the matched anchor rows
bio = np.array([r["seqid_anchor"] for r in rows])
fid = np.array([r["fident_anchor"] for r in rows if r["fident_anchor"] is not None])
bio_m = np.array([r["seqid_anchor"] for r in rows if r["fident_anchor"] is not None])
fs_corr = stats.pearsonr(bio_m, fid)

summary = dict(
    n_screened=len(rows), n_triad=len(triad), n_above=len(above),
    floor_rate=FLOOR, cov_min=COV_MIN,
    foldseek_crosscheck=dict(
        pearson_r=round(fs_corr.statistic, 4), n=len(fid),
        mean_biotite=round(float(bio_m.mean()), 2), mean_foldseek=round(float(fid.mean()), 2),
        median_abs_diff=round(float(np.median(np.abs(bio_m - fid))), 2)),
    triage_counts_foldseek_fident=triage_fident_counts,
    structure_first_foldseek_fident=[dict(
        mgyp=r["mgyp"], anchor=r["anchor"], fident_nearest=r["fident_nearest"],
        seqid_nearest=r["seqid_nearest"], cov_nearest=r["cov_nearest"],
        composite=r["composite"], mean_plddt=r["mean_plddt"], best_bits=r["best_bits"])
        for r in sf_fident],
    seqid_strat=strat,
    seqid_distribution=dict(
        triad=dict(min=float(sid.min()), p25=float(np.percentile(sid, 25)),
                   median=float(np.median(sid)), p75=float(np.percentile(sid, 75)),
                   max=float(sid.max())),
        above_line=dict(
            min=float(np.min([r["seqid_nearest"] for r in above])),
            median=float(np.median([r["seqid_nearest"] for r in above])),
            max=float(np.max([r["seqid_nearest"] for r in above]))),
    ),
    correlation=dict(pearson_r=round(pear.statistic, 4), pearson_p=pear.pvalue,
                     spearman_rho=round(spear.statistic, 4), spearman_p=spear.pvalue,
                     median_seqid_at_1090bits=sid_at_1090),
    triage_counts={t: len(v) for t, v in triage.items()},
    triage_above_seqid={t: sorted([r["seqid_nearest"] for r in v]) for t, v in triage.items()},
    structure_first=[dict(mgyp=r["mgyp"], anchor=r["anchor"], nearest_query=r["nearest_query"],
                          seqid_nearest=r["seqid_nearest"], cov_nearest=r["cov_nearest"],
                          seqid_anchor=r["seqid_anchor"], cov_anchor=r["cov_anchor"],
                          composite=r["composite"], mean_plddt=r["mean_plddt"],
                          best_bits=r["best_bits"], low_coverage=r["low_coverage"]) for r in sf],
    low_coverage_count=sum(r["low_coverage"] for r in rows),
)
json.dump(summary, open(f"{TMP}/seqid_summary.json", "w"), indent=2)

# ---- console report ----
print("\n==== CP2: above-line|triad by NEAREST-PETase sequence identity ====")
print(f"{'bin':>8} {'n_triad':>7} {'k_above':>7} {'rate':>7} {'vs floor 42.9%':>22}")
for s in strat:
    rate = f"{100*s['rate']:.1f}%" if s["rate"] is not None else "—"
    vs = (f"p={s['fisher_p_vs_floor']:.3f} RR={s['rr_vs_floor']:.2f}x"
          if s["n_triad"] else "—")
    print(f"{s['bin']:>8} {s['n_triad']:>7} {s['k_above']:>7} {rate:>7} {vs:>22}")

print("\n==== bits <-> seq-id correlation (294 triad-bearers) ====")
print(f"Pearson r={pear.statistic:.3f} (p={pear.pvalue:.2e}); "
      f"Spearman rho={spear.statistic:.3f} (p={spear.pvalue:.2e})")
print(f"median nearest-PETase seq-id at ~1090 bits (top-100 cut): {sid_at_1090}")

print("\n==== seq-id distribution ====")
print("  all 294 triad nearest-id: min %.1f  p25 %.1f  median %.1f  p75 %.1f  max %.1f"
      % (sid.min(), np.percentile(sid,25), np.median(sid), np.percentile(sid,75), sid.max()))
av = [r["seqid_nearest"] for r in above]
print("  96 above-line nearest-id: min %.1f  median %.1f  max %.1f"
      % (min(av), np.median(av), max(av)))

print("\n==== Foldseek fident cross-check (independent structure-derived id) ====")
fc = summary["foldseek_crosscheck"]
print(f"  biotite vs Foldseek anchor identity: Pearson r={fc['pearson_r']} "
      f"(mean {fc['mean_biotite']}% vs {fc['mean_foldseek']}%, median|diff| {fc['median_abs_diff']}pts)")

print("\n==== CP3: triage of the 96 above-line ====")
print("  [primary metric = biotite Smith-Waterman, BLOSUM62]")
for t in ("near_homolog", "twilight", "structure_first"):
    print(f"  {t:>16}: {len(triage[t])}")
print(f"  low-coverage flagged (any of 300): {summary['low_coverage_count']}")
print("  [robustness = Foldseek structural fident]")
for t in ("near_homolog", "twilight", "structure_first"):
    print(f"  {t:>16}: {triage_fident_counts[t]}")
if sf_fident:
    print("  structure-first under Foldseek fident (<25%):")
    for r in sf_fident:
        print(f"    {r['mgyp']} anchor={r['anchor']} fident={r['fident_nearest']}% "
              f"(biotite {r['seqid_nearest']}%) comp={r['composite']} pLDDT={r['mean_plddt']}")
print("\n  structure-first finds (<25% id OR no credible homolog), above-line:")
if sf:
    print(f"  {'MGYP':>18} {'anchor':>6} {'nearest':>7} {'id%':>5} {'covQ':>5} {'comp':>7} {'pLDDT':>6} {'bits':>6} {'lowcov'}")
    for r in sf:
        print(f"  {r['mgyp']:>18} {r['anchor']:>6} {r['nearest_query']:>7} "
              f"{r['seqid_nearest']:>5.1f} {r['cov_nearest']:>5.2f} {r['composite']:>7.2f} "
              f"{r['mean_plddt']:>6.3f} {r['best_bits']:>6.0f} {r['low_coverage']}")
else:
    print("  (none — every above-line hit is >=25% identity to a PETase query with credible coverage)")
print("\n[done] summary ->", f"{TMP}/seqid_summary.json")
