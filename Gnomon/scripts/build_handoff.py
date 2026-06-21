"""One command that regenerates every number in the retool hand-off.

Runs the two seeded pipelines and refreshes the machine-readable artifacts the
hand-off docs cite:

  * ``gnomon.reproduce.run``  -> results/reproduction.json   (CP3 keep-set + verdict)
  * ``gnomon.reframe.run_reframe`` -> handoff/conditional_coverage.json  (CP4.1 reframe)

The Markdown hand-off (``RETOOL_HANDOFF.md``, ``handoff/CLAIMS_LEDGER.md``,
``docs/METHODS.md``, ``VERDICT.md``) is static prose whose numbers all trace to these
two JSON files. Re-run this, then the docs are current.

Usage:
    KMP_DUPLICATE_LIB_OK=TRUE PYTHONPATH=. python scripts/build_handoff.py            # full
    KMP_DUPLICATE_LIB_OK=TRUE PYTHONPATH=. python scripts/build_handoff.py --fast     # skip flow+real
"""
from __future__ import annotations

import os
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import argparse

from gnomon import reproduce, reframe


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true",
                    help="skip MAF training + OSIPI fetch (keep-set T4/T1 omitted)")
    ap.add_argument("--n-noise", type=int, default=None,
                    help="override noise realizations (default: manifest N_NOISE)")
    args = ap.parse_args()

    kw = {} if args.n_noise is None else {"n_noise": args.n_noise}

    print("== [1/2] reproduce: keep-set + verdict ==")
    rep = reproduce.run(run_flow=not args.fast, run_real=not args.fast, **kw)

    print("\n== [2/2] reframe: conditional-coverage table (both SD conventions) ==")
    rf = reframe.run_reframe(**kw)

    print("\n== hand-off regenerated ==")
    print(f"  verdict: {rep['verdict']}  (divergence: {rep.get('divergence')})")
    print(f"  reframe reconstruction_ok: {rf['reconstruction_ok']}")
    print("  artifacts: results/reproduction.json, handoff/conditional_coverage.json")


if __name__ == "__main__":
    main()
