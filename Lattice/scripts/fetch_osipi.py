"""CLI: fetch the external OSIPI TF2.4 reference DRO on demand + write provenance.

The OSIPI DRO is CC-BY-4.0 data hosted on Zenodo; it is downloaded on demand and
NEVER committed to this repository (Lattice ships synthetic data only). This
script wraps lattice.osipi.fetch and records a provenance manifest.

Usage:
    python scripts/fetch_osipi.py --dest data/osipi --fetched 2026-06-19T00:00:00Z
"""

import argparse
import json
from pathlib import Path

from lattice import osipi


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dest", default="data/osipi", help="download destination dir")
    ap.add_argument(
        "--fetched",
        default=None,
        help="ISO timestamp recorded in the provenance manifest (caller-supplied; "
        "this tool never reads the wall clock)",
    )
    args = ap.parse_args()

    print("OSIPI source:")
    print(json.dumps(osipi.OSIPI_SOURCE, indent=2))
    print(f"\nDestination: {Path(args.dest).resolve()}  (git-ignored; not redistributed)")

    try:
        path = osipi.fetch(args.dest, fetched=args.fetched)
        print(f"Fetched: {path}")
    except NotImplementedError as exc:
        print(f"\n[not wired] {exc}")
        print("Provenance record template for the manifest:")
        print(json.dumps(osipi.provenance_record(args.dest, fetched=args.fetched), indent=2))


if __name__ == "__main__":
    main()
