"""Ferry — a real-data substrate adapter for Matrix's synthetic-twin closed loop.

Ferry grounds Matrix on **real anatomy + real dose geometry** from a public RT dataset
(TCIA Pancreatic-CT-CBCT-SEG, CC BY 4.0), to pre-empt the "shown only on a pure synthetic
twin" objection. It is an **interface-swap only**: a :class:`GroundedTwin` drops into the
existing ``run_iteration`` engine, so ``loop.py`` is byte-unchanged (see ``verify_ferry_cp1.py``).

**Honest ceiling (read this).** Real anatomy + real dose geometry; **perfusion/IVIM stays
synthetic** — there is no scanner, hence no real diffusion data. A grounded result means
"the loop closes on *real geometry*", never a real-IVIM / clinical claim. See ``FERRY.md``.
"""
from __future__ import annotations

from .substrate import FerrySubstrate, GroundedTwin, rescale_dose_to_band
from .dataset import load_substrate, FerryDataUnavailable, COLLECTION, DOI, LICENSE
from .loop_grounded import run_grounded_loop

__all__ = ["FerrySubstrate", "GroundedTwin", "rescale_dose_to_band",
           "load_substrate", "FerryDataUnavailable", "run_grounded_loop",
           "COLLECTION", "DOI", "LICENSE"]
