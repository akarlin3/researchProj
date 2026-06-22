"""Matrix — a synthetic-twin closed-loop harness for adaptive quantitative-MRI dosing.

Matrix is Keystone's *no-scanner* run-mode, built standalone and early to de-risk lab
access: it closes the loop ``scan -> posterior -> trust gate -> action gate -> dose
replan -> re-scan`` on a purely **synthetic digital twin**, with **no scanner** and **no
real patient data** anywhere.

The three components Matrix consumes are each stubbed behind a clean interface with a
clearly-labelled placeholder (NOT-Fashion / NOT-Minos / NOT-Forge); the real component
drops in without touching the loop. See ``ASSUMPTIONS.md``.

This is a *working closed-loop harness on synthetic data*, never a validated clinical
loop. Every result means "the loop closes + behaves sensibly on a synthetic twin."
"""
from __future__ import annotations

from .config import (MatrixConfig, DEFAULT, NORMAL, TUMOR, OAR,
                     SPARE, TREAT, ESCALATE, ACTION_NAMES)
from .twin import Twin
from .state import LoopState
from .loop import Interfaces, run_loop, run_iteration

__all__ = ["MatrixConfig", "DEFAULT", "NORMAL", "TUMOR", "OAR",
           "SPARE", "TREAT", "ESCALATE", "ACTION_NAMES",
           "Twin", "LoopState", "Interfaces", "run_loop", "run_iteration"]
