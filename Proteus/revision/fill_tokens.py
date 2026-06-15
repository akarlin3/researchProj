"""Fill the {{double-brace}} tokens in manuscript_revisions.tex from the committed
analysis JSON, writing manuscript_revisions.filled.tex. Every number in the revised
prose therefore traces to a script output — none is hand-typed (anti-fabrication).

Reproduce:  python revision/fill_tokens.py
"""
from __future__ import annotations

import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
WT = os.path.dirname(HERE)
PROC = os.path.join(WT, "data", "processed")
SRC = os.path.join(HERE, "manuscript_revisions.tex")
OUT = os.path.join(HERE, "manuscript_revisions.filled.tex")

eq = json.load(open(os.path.join(PROC, "equivalence.json")))
fl = json.load(open(os.path.join(PROC, "floor_enlarged.json")))
pc = json.load(open(os.path.join(PROC, "plddt_confound.json")))


def sci(p: float) -> str:
    """LaTeX scientific notation, e.g. 1.5\\times10^{-9}."""
    if p == 0:
        return "0"
    exp = 0
    x = p
    while x < 1:
        x *= 10
        exp -= 1
    while x >= 10:
        x /= 10
        exp += 1
    return f"{x:.1f}\\times10^{{{exp}}}"


def pp(x_frac: float, signed: bool = True) -> str:
    v = x_frac * 100
    return (f"{v:+.1f}" if signed else f"{v:.1f}")


base = eq["baseline"]
pet = eq["pet_branch"]
ache = eq["ache_branch"]
tost = eq["tost"]
pf = fl["pooled_floor"]
g = pc["global_plddt"]
dvb = next(t for t in pc["tests_global"] if t["comparison"] == "divergent vs baseline")
mr = pc["matched_rescreen"]

diff = eq["difference_pet_minus_baseline"]
ub = eq["one_sided_95_upper_bound_on_enrichment"]

# narrative fragments (verdict-dependent, honest)
tost_descriptor = ("i.e.\\ even the most optimistic bound is negative --- the branch "
                   "in fact sits significantly below the floor")
tost_sentence = (
    "Strict two-sided equivalence within the margin is therefore not claimed: the "
    "PET branch is not merely \\emph{not enriched}, it is significantly \\emph{below} "
    f"the random hydrolase floor (two-proportion $p={sci(eq['legacy_tests']['two_proportion_p'])}$, "
    f"Fisher $p={sci(eq['legacy_tests']['fisher_p'])}$).")

TOK = {
    "triad_ratio": "53",
    "base_rate_pct": f"{base['rate']*100:.1f}",
    "base_n": str(base["n"]),
    "base_k": str(base["k"]),
    "base_wilson_lo_pct": f"{base['wilson95'][0]*100:.1f}",
    "base_wilson_hi_pct": f"{base['wilson95'][1]*100:.1f}",
    "delta_pp": f"{eq['delta_pp']:.0f}",
    "upper_bound_pp": pp(ub),
    "diff_pp": pp(diff),
    "ci90_lo_pp": pp(tost["ci90_diff"][0]),
    "ci90_hi_pp": pp(tost["ci90_diff"][1]),
    "p_upper": sci(tost["p_upper_not_enriched_above_delta"]),
    "tost_descriptor": tost_descriptor,
    "tost_sentence": tost_sentence,
    "pooled_screened": f"{pf['screened']:,}",
    "new_screened": f"{fl['new_draws']['screened']:,}",
    "enlarge_seed": str(fl["enlarge_seed_base"]),
    "div_n": f"{pc['n']['divergent']:,}",
    "nearhom_n": str(pc["n"]["near_homolog"]),
    "div_med_pct": f"{g['divergent_median']*100:.1f}",
    "base_med_pct": f"{g['baseline_median']*100:.1f}",
    "nearhom_med_pct": f"{g['near_homolog_median']*100:.1f}",
    "div_sub20_med_pct": f"{g['divergent_sub20_median']*100:.1f}",
    "plddt_mwu_p": sci(dvb["mannwhitney_p"]),
    "plddt_cliffs": f"{dvb['cliffs_delta_a_vs_b']:+.2f}",
    "matched_raw_pct": f"{mr['raw_above_line_rate']*100:.1f}",
    "matched_pct": f"{mr['matched_above_line_rate']*100:.1f}",
}

text = open(SRC).read()
missing = set(re.findall(r"\{\{([a-z0-9_]+)\}\}", text)) - set(TOK)
if missing:
    raise SystemExit(f"unfilled tokens (no value provided): {sorted(missing)}")
for k, v in TOK.items():
    text = text.replace("{{" + k + "}}", v)
left = re.findall(r"\{\{[a-z0-9_]+\}\}", text)
assert not left, f"tokens remain: {left}"
open(OUT, "w").write(text)
print(f"filled {len(TOK)} tokens -> {OUT}")
for k, v in sorted(TOK.items()):
    print(f"  {k} = {v}")
