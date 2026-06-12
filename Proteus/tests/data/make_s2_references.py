#!/usr/bin/env python3
"""Generate the deterministic known-answer S2 fold-class reference fixtures.

S2 triages the S1 query DB against the alpha/beta-hydrolase FOLD CLASS using BOTH
references (the decision recorded in envlog/env-failures.md):

  * CURATED reference (precision anchors) — a small representative set of real
    alpha/beta-hydrolase structures. Stands in for the curated ESTHER /
    representative set in production.
  * BROAD reference (recall + post-filter) — the same hydrolases PLUS non-fold
    distractor structures, mimicking a broad AF-DB/PDB search that returns
    everything. A post-filter allow-list keeps only the alpha/beta-hydrolase-fold
    targets, so spurious hits to distractors do NOT shortlist a query.

NONE of these reference proteins is among the 6 query controls in
mini_corpus.fasta, so a query passes by FOLD similarity, never by being identical
to a reference (true fold-class triage, not a seeded-template/homology gate).

CURATED set (6 real alpha/beta-hydrolase chains, from RCSB):
    1CEX  Fusarium solani cutinase
    1EDE  Xanthobacter haloalkane dehalogenase  (archetypal alpha/beta-hydrolase)
    1I6W  Bacillus subtilis lipase A            (minimal alpha/beta-hydrolase)
    1AUO  Pseudomonas esterase
    3TGL  Rhizomucor miehei lipase
    1QLW  Streptomyces lipase

BROAD set = CURATED + 2 non-alpha/beta-hydrolase-fold DISTRACTORS:
    1MBN  sperm-whale myoglobin   (all-alpha)
    1UBQ  ubiquitin               (beta-grasp)
The post-filter allow-list (s2_ref_broad.allow.txt) names ONLY the 6 hydrolases.

EXPECTED S2 OUTCOME on the mini corpus (known answer the test asserts)
---------------------------------------------------------------------
S0 yields 8 representatives; S1 tokenizes them to 3Di; S2 searches both refs:
    SHORTLISTED (alpha/beta-hydrolase fold): IsPETase, LCC_WT, CalB, AChE, CRL, Est2
    DROPPED     (no fold-class hit):          decoy_allalpha, decoy_random
    => 8 representatives -> 6 shortlisted, 2 dropped.

Note the 4 NEGATIVE controls (CalB/AChE/CRL/Est2) are NOT PETases yet ARE
shortlisted: they share the alpha/beta-hydrolase fold and are only separated
later at S4 (catalytic geometry) / S5 (cleft). That is the point of fold-class
triage — keep the whole fold class, including the divergent dark tail.
"""
from __future__ import annotations

import os
import textwrap

HERE = os.path.dirname(os.path.abspath(__file__))
CURATED = os.path.join(HERE, "s2_ref_curated.fasta")
BROAD = os.path.join(HERE, "s2_ref_broad.fasta")
ALLOW = os.path.join(HERE, "s2_ref_broad.allow.txt")

# Real alpha/beta-hydrolase chain-A sequences (RCSB), distinct from the controls.
HYDROLASES = {
    "1CEX": (
        "LPTSNPAQELEARQLGRTTRDDLINGNSASCADVIFIYARGSTETGNLGTLGPSIASNLESAFGKDGVW"
        "IQGVGGAYRATLGDNALPRGTSSAAIREMLGLFQQANTKCPDATLIAGGYSQGAALAAASIEDLDSAIR"
        "DKIAGTVLFGYTKNLQNRGRIPNYPADRTKVFCNTGDLVCTGSLIVAAPHLAYGPDARGPAPEFLIEKV"
        "RAVRGSA"
    ),
    "1EDE": (
        "MINAIRTPDQRFSNLDQYPFSPNYLDDLPGYPGLRAHYLDEGNSDAEDVFLCLHGEPTWSYLYRKMIPV"
        "FAESGARVIAPDFFGFGKSDKPVDEEDYTFEFHRNFLLALIERLDLRNITLVVQDWGGFLGLTLPMADPS"
        "RFKRLIIMNACLMTDPVTQPAFSAFVTQPADGFTAWKYDLVTPSDLRLDQFMKRWAPTLTEAEASAYAAP"
        "FPDTSYQAGVRKFPKMVAQRDQACIDISTEAISFWQNDWNGQTFMAIGMKDKLLGPDVMYPMKALINGCP"
        "EPLEIADAGHFVQEFGEQVAREALKHFAETE"
    ),
    "1I6W": (
        "AEHNPVVMVHGIGGASFNFAGIKSYLVSQGWSRDKLYAVDFWDKTGTNYNNGPVLSRFVQKVLDETGAKK"
        "VDIVAHSMGGANTLYYIKNLDGGNKVANVVTLGGANRLTTGKALPGTDPNQKILYTSIYSSADMIVMNYL"
        "SRLDGARNVQIHGVGHIGLLYSSQVNSLIKEGLNGGGQNTN"
    ),
    "1AUO": (
        "MTEPLILQPAKPADACVIWLHGLGADRYDFMPVAEALQESLLTTRFVLPQAPTRPVTINGGYEMPSWYDI"
        "KAMSPARSISLEELEVSAKMVTDLIEAQKRTGIDASRIFLAGFSQGGAVVFHTAFINWQGPLGGVIALST"
        "YAPTFGDELELSASQQRIPALCLHGQYDDVVQNAMGRSAFEHLKSRGVTVTWQEYPMGHEVLPQEIHDIG"
        "AWLAARLG"
    ),
    "3TGL": (
        "SINGGIRAATSQEINELTYYTTLSANSYCRTVIPGATWDCIHCDATEDLKIIKTWSTLIYDTNAMVARGD"
        "SEKTIYIVFRGSSSIRNWIADLTFVPVSYPPVSGTKVHKGFLDSYGEVQNELVATVLDQFKQYPSYKVAV"
        "TGHSLGGATVLLCALDLYQREEGLSSSNLFLYTQGQPRVGDPAFANYVVSTGIPYRRTVNERDIVPHLPP"
        "AAFGFLHAGEEYWITDNSPETVQVCTSDLETSDCSNSIVPFTSVLDHLSYFGINTGLCT"
    ),
    "1QLW": (
        "APPPVPKTPAGPLTLSGQGSFFVGGRDVTSETLSLSPKYDAHGTVTVDQMYVRYQIPQRAKRYPITLIHG"
        "CCLTGMTWETTPDGRMGWDEYFLRKGYSTYVIDQSGRGRSATDISAINAVKLGKAPASSLPDLFAAGHEA"
        "AWAIFRFGPRYPDAFKDTQFPVQAQAELWQQMVPDWLGSMPTPNPTVANLSKLAIKLDGTVLLSHSQSGI"
        "YPFQTAAMNPKGITAIVSVEPGECPKPEDVKPLTSIPVLVVFGDHIEEFPRWAPRLKACHAFIDALNAAG"
        "GKGQLMSLPALGVHGNSHMMMQDRNNLQVADLILDWIGRNTAKPAHGR"
    ),
}

# Non-alpha/beta-hydrolase-fold distractors — present ONLY in the broad reference;
# the post-filter allow-list excludes them.
DISTRACTORS = {
    "1MBN": (  # sperm-whale myoglobin (all-alpha)
        "VLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKHGVTVLT"
        "ALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRK"
        "DIAAKYKELGYQG"
    ),
    "1UBQ": (  # ubiquitin (beta-grasp)
        "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLV"
        "LRLRGG"
    ),
}


def _write_fasta(path: str, records: dict[str, str]) -> None:
    with open(path, "w") as fh:
        for rid, seq in records.items():
            fh.write(f">{rid}\n")
            for line in textwrap.wrap(seq, 60):
                fh.write(line + "\n")


if __name__ == "__main__":
    _write_fasta(CURATED, HYDROLASES)
    broad = {**HYDROLASES, **DISTRACTORS}
    _write_fasta(BROAD, broad)
    with open(ALLOW, "w") as fh:
        fh.write("# alpha/beta-hydrolase-fold targets in s2_ref_broad.fasta "
                 "(post-filter allow-list)\n")
        for rid in HYDROLASES:
            fh.write(rid + "\n")
    print(f"curated: {len(HYDROLASES)} hydrolases -> {os.path.relpath(CURATED)}")
    print(f"broad:   {len(broad)} ({len(DISTRACTORS)} distractors) -> "
          f"{os.path.relpath(BROAD)}")
    print(f"allow:   {len(HYDROLASES)} ids -> {os.path.relpath(ALLOW)}")
