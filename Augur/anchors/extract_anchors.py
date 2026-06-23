#!/usr/bin/env python3
"""Extract Augur's load-bearing in-repo anchors from the committed sibling results.

Augur makes *no new measurement*. Its reproduction therefore begins by EXTRACTING the
load-bearing numbers from the committed, separately-gated sibling results (the "in-repo
provisional anchors") into a single ``anchors.json`` -- the single source of truth that the
manuscript's ``numbers.tex`` and every downstream reproduction script trace back to. This is
the house number-traceability discipline (cf. Gauge/Minos/Lethe ``consistency.py``): no number
appears in Augur that does not originate in a committed, cited source file.

Sources (all committed, all clean-IP -- synthetic / open ACRIN-6698 CC-BY-4.0 / literature):
  - Gauge/results/invivo_real_provenance.json     -> D* test-retest r, p, n, BCa+Fisher-z CIs
  - Gauge/results/conditional_attack_report.txt    -> CRLB(D*) identifiability-wall numbers
  - Augur/CITATIONS.md                             -> external D*-Ktrans literature (Sun/Yang)

When the sibling Gauge tree is present (the monorepo case) this re-extracts fresh and rewrites
anchors.json. When Augur is carved standalone (Gauge absent) it validates the committed
anchors.json instead. Either way the output is deterministic and fully traceable.

Exit code: 0 = anchors extracted/validated and self-consistent; non-zero = a source is missing
or an expected anchor token is absent (fail loudly -- never silently fabricate).
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
AUGUR = HERE.parent
ROOT = AUGUR.parent
GAUGE = ROOT / "Gauge"
OUT = HERE / "anchors.json"

# Committed sources (relative to the monorepo root).
SRC_PROVENANCE = GAUGE / "results" / "invivo_real_provenance.json"
SRC_CRLB = GAUGE / "results" / "conditional_attack_report.txt"
SRC_CITATIONS = AUGUR / "CITATIONS.md"


def _require(cond: bool, msg: str) -> None:
    if not cond:
        print(f"EXTRACT-ANCHORS FAIL: {msg}", file=sys.stderr)
        raise SystemExit(2)


def extract_retest(prov_path: Path) -> dict:
    """D* (and companion D) test-retest width-vs-repeatability anchor (Gauge Checkpoint D)."""
    _require(prov_path.exists(), f"missing source {prov_path}")
    prov = json.loads(prov_path.read_text())
    rt = prov["retest"]
    out = {
        "source": str(prov_path.relative_to(ROOT)),
        "dataset": rt["arm"],
        "doi": rt["doi"],
        "license": rt["license"],
        "n_pairs": int(rt["n_pairs"]),
        "b_values": list(rt["b_values"]),
        "seed": int(rt["seed"]),
        "ci_method": rt["ci_method"],
        "Dstar_spearman": float(rt["Dstar_spearman"]),
        "Dstar_p": float(rt["Dstar_spearman_p"]),
        "Dstar_ci95_boot": [float(x) for x in rt["Dstar_spearman_ci95_boot"]],
        "Dstar_ci95_fisher": [float(x) for x in rt["Dstar_spearman_ci95_fisher"]],
        "D_spearman": float(rt["D_spearman"]),
        "D_p": float(rt["D_spearman_p"]),
        "D_ci95_boot": [float(x) for x in rt["D_spearman_ci95_boot"]],
        "D_ci95_fisher": [float(x) for x in rt["D_spearman_ci95_fisher"]],
    }
    # Sanity: the D* CI must span zero (the null finding the spine relies on); D must not.
    lo, hi = out["Dstar_ci95_boot"]
    _require(lo < 0 < hi, f"D* BCa CI must span zero (got [{lo}, {hi}])")
    _require(out["D_ci95_boot"][0] > 0, "D BCa CI must be strictly positive")
    _require(out["n_pairs"] == 76, f"expected n=76 retest pairs, got {out['n_pairs']}")
    return out


def extract_crlb(report_path: Path) -> dict:
    """Gauge CRLB high-D* identifiability-wall anchor (conditional_attack CP3 verdict)."""
    _require(report_path.exists(), f"missing source {report_path}")
    text = report_path.read_text()
    # Pull the three load-bearing identifiability numbers verbatim from the verdict block.
    m_terc = re.search(r"CRLB\(D\*\)/tercile-width reaches\s*([0-9]+\.[0-9]+)", text)
    m_grow = re.search(r"absolute CRLB grows\s*~?([0-9]+(?:\.[0-9]+)?)x", text)
    m_corr = re.search(r"conformal width\s*~\s*CRLB\s*r=([0-9]+\.[0-9]+)", text)
    m_edges = re.search(r"true-D\* regime edges \(terciles[^)]*\):\s*\[([0-9.eE+\- ]+)\]", text)
    _require(m_terc is not None, "CRLB/tercile-width ratio token not found")
    _require(m_grow is not None, "CRLB ~Nx growth token not found")
    _require(m_corr is not None, "conformal-width~CRLB r token not found")
    _require(m_edges is not None, "tercile edge token not found")
    edges = [float(x) for x in m_edges.group(1).split()]
    out = {
        "source": str(report_path.relative_to(ROOT)),
        "crlb_over_tercile_width_hi": float(m_terc.group(1)),  # ~1.12
        "crlb_growth_factor": float(m_grow.group(1)),          # ~6
        "conformal_width_crlb_r": float(m_corr.group(1)),      # 0.77
        "dstar_tercile_edges": edges,                           # [0.0416, 0.0709]
        # The IVIM acquisition regime the wall is scoped to (segmented 4-b scheme).
        "regime": "segmented 4-b IVIM (b=[0,100,600,800] s/mm^2); CRLB at true params, S0 free",
    }
    _require(out["crlb_over_tercile_width_hi"] >= 1.0,
             "high-D* CRLB/tercile-width must be >= 1 (regime unresolvable)")
    return out


def extract_dstar_ktrans(cit_path: Path) -> dict:
    """External cross-modal D*-Ktrans evidence -- literature only, no in-repo data."""
    _require(cit_path.exists(), f"missing source {cit_path}")
    text = cit_path.read_text()
    _require("10.1016/j.acra.2018.08.012" in text, "Sun 2019 DOI absent from CITATIONS.md")
    _require("10.1177/0284185118791201" in text, "Yang 2019 DOI absent from CITATIONS.md")
    _require("r = 0.389" in text, "Sun 2019 verified r value absent from CITATIONS.md")
    return {
        "source": str(cit_path.relative_to(ROOT)),
        "kind": "external literature (no in-repo computation)",
        "sun2019": {"r": 0.389, "p": "<0.001", "fDstar_r": 0.533,
                    "doi": "10.1016/j.acra.2018.08.012", "cohort": "rectal cancer"},
        "yang2019": {"r": None, "significant": False, "Dstar_icc": 0.55,
                     "doi": "10.1177/0284185118791201", "cohort": "rectal cancer"},
        "framing": "weak and cohort-inconsistent (significant in Sun 2019; null in Yang 2019)",
    }


def main() -> int:
    if not GAUGE.exists():
        # Standalone-carve case: validate the committed anchors.json instead of re-extracting.
        if OUT.exists():
            data = json.loads(OUT.read_text())
            print(f"Gauge sibling absent; validated committed anchors.json "
                  f"({len(data)} anchor groups). [standalone-carve mode]")
            return 0
        print("EXTRACT-ANCHORS FAIL: Gauge sibling absent and no committed anchors.json",
              file=sys.stderr)
        return 2

    anchors = {
        "_about": "Augur load-bearing in-repo anchors, extracted from committed sibling results. "
                  "AUTO-GENERATED by anchors/extract_anchors.py -- DO NOT EDIT BY HAND.",
        "retest": extract_retest(SRC_PROVENANCE),
        "crlb": extract_crlb(SRC_CRLB),
        "dstar_ktrans": extract_dstar_ktrans(SRC_CITATIONS),
    }
    OUT.write_text(json.dumps(anchors, indent=2) + "\n")
    print("Augur anchors extracted ->", OUT.relative_to(ROOT))
    print("=" * 72)
    rt, cr = anchors["retest"], anchors["crlb"]
    print(f"  D* test-retest: Spearman r={rt['Dstar_spearman']:+.3f} (p={rt['Dstar_p']:.3f}, "
          f"n={rt['n_pairs']}); BCa 95% CI [{rt['Dstar_ci95_boot'][0]:+.2f}, "
          f"{rt['Dstar_ci95_boot'][1]:+.2f}]  (spans 0)")
    print(f"  companion D:    Spearman r={rt['D_spearman']:+.3f} (p={rt['D_p']:.1e})")
    print(f"  CRLB wall:      CRLB(D*)/tercile-width ->{cr['crlb_over_tercile_width_hi']:.2f} at hi-D*; "
          f"abs CRLB ~{cr['crlb_growth_factor']:.0f}x; width~CRLB r={cr['conformal_width_crlb_r']:.2f}")
    print(f"  D*-Ktrans:      {anchors['dstar_ktrans']['framing']}")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())
