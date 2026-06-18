"""Gauge reproduction -- one command, fixed seed, numpy only.

Reproduces, on synthetic IVIM data, the headline phenomenon of the Gauge
manuscript (in review at MRM, 2026; pre-publication -- no publication DOI):

    marginal conformal restores pooled D* coverage, but the high-D* tercile
    stays under-covered (the identifiability wall), and group-conditional
    (Mondrian) correction buys per-tercile coverage back only by inflating width.

This is a qualitative reproduction using only caliper.conformal + caliper.metrics
-- NOT a published or independently validated result. See caliper.publication for
the paper's true (pre-publication) status.

Run:  python examples/gauge_repro.py        (no torch required)
"""
from __future__ import annotations

from caliper.repro_gauge import reproduce

_BANNER = (
    "NOTE: Gauge is pre-publication (in review at Magnetic Resonance in Medicine,\n"
    "      no publication DOI yet). The numbers below are a SYNTHETIC, qualitative\n"
    "      reproduction of the paper's phenomenon -- not a published result.\n"
)


def main() -> None:
    print(_BANNER)
    result = reproduce()
    print(result.format())


if __name__ == "__main__":
    main()
