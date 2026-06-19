"""PROVISIONAL stamping for ruler-dependent reference numbers.

The cardinal rule of Datum: **never present a ruler-dependent reference number as
final.** Any value produced by scoring through Fashion's (in-review) ruler is
wrapped/stamped PROVISIONAL so it cannot be silently mistaken for a locked number.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from datum.manifest import RULER, is_provisional

PROVISIONAL_BANNER = (
    "+-----------------------------------------------------------------------+\n"
    "|  PROVISIONAL -- numbers below are scored on Fashion's calibration      |\n"
    "|  ruler, which is IN REVIEW. They are NOT final reference values.       |\n"
    "|  Re-run `python revalidate.py` when the ruler locks (DOI assigned).    |\n"
    "+-----------------------------------------------------------------------+"
)


@dataclass
class Provisional:
    """A reference number that is provisional until the ruler locks."""
    value: Any
    metric: str
    reason: str = field(default="ruler-dependent: Fashion ruler in review")

    @property
    def provisional(self) -> bool:
        return is_provisional()

    def __repr__(self) -> str:  # always reads as provisional in logs/tables
        tag = "PROVISIONAL" if self.provisional else "final"
        return f"<{self.metric}={self.value!r} [{tag}]>"


def stamp(value: Any, metric: str, reason: str | None = None) -> Provisional:
    """Wrap a ruler-derived value as Provisional."""
    return Provisional(value=value, metric=metric,
                       reason=reason or f"scored on Fashion ruler {RULER['version']} "
                                        f"({RULER['manuscript_status']})")
