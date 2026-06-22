"""CP3 gate — the Forge-shaped dose-replan stage consumes the gated decision.

Checks (all green => CP3 passes):

  (1) CONSUMES+REPLANS  the dose stage consumes the gated per-voxel decision and returns
                        a replan on the twin: TREAT voxels boost, SPARE de-escalate,
                        ESCALATE hold; the prescription stays within [dose_min, dose_max].
  (2) FORGE-DROP-IN     the interface is a clean drop-in: an *alternative* engine
                        implementing the same ``replan(current_dose, action, state, cfg)``
                        signature swaps in and the loop runs unchanged — proving loop.py
                        never needs editing when Forge's real engine lands.
  (3) ISOLATED+LABELLED the placeholder is analytic (no geometry/transport/ERE), clearly
                        labelled NOT-Forge, and isolated to the one object.

Run: <proteus python> Matrix/verify_cp3.py
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from matrix import MatrixConfig, Twin, Interfaces, run_loop, TREAT, SPARE, ESCALATE
from matrix.interfaces.dose import PlaceholderDoseEngine


def _hr(t): print("\n" + "=" * 74 + f"\n{t}\n" + "=" * 74)


class _ConstantDoseEngine:
    """An alternative engine with the SAME interface — used only to prove drop-in."""
    name = "constant-dose"
    label = "NOT-Forge (drop-in proof engine)"
    provisional = True

    def replan(self, current_dose, action, state, cfg):
        dose = np.asarray(current_dose, float).copy()
        return dict(dose=dose, delta=np.zeros_like(dose))


def main() -> int:
    print("CP3 verification — Forge-shaped dose-replan stage")
    cfg = MatrixConfig()

    _hr("CP3 check 1/3 — stage consumes the gated decision and returns a replan")
    twin, states = run_loop(cfg, Interfaces.placeholders(), n_iter=1)
    s = states[0]
    eng = PlaceholderDoseEngine()
    print(f"  dose engine = {eng.label}")
    d = s.delta_dose
    print(f"  TREAT voxels: n={int(np.sum(s.action==TREAT))}  all boosted>0: "
          f"{bool(np.all(d[s.action==TREAT] > 0))}")
    print(f"  SPARE voxels: n={int(np.sum(s.action==SPARE))}  all de-escalated<=0: "
          f"{bool(np.all(d[s.action==SPARE] <= 0))}")
    print(f"  ESCALATE voxels: n={int(np.sum(s.action==ESCALATE))}  all held==0: "
          f"{bool(np.all(d[s.action==ESCALATE] == 0))}")
    print(f"  replan within [dose_min,dose_max]: "
          f"{bool((s.replan.min()>=cfg.dose_min) and (s.replan.max()<=cfg.dose_max))}")
    assert np.all(d[s.action == TREAT] > 0)
    assert np.all(d[s.action == SPARE] <= 0)
    assert np.all(d[s.action == ESCALATE] == 0)
    assert s.replan.min() >= cfg.dose_min and s.replan.max() <= cfg.dose_max
    print("  CONSUMES+REPLANS: PASS")

    _hr("CP3 check 2/3 — Forge drop-in: swap the engine, loop.py untouched")
    swapped = Interfaces.placeholders()
    swapped.dose_engine = _ConstantDoseEngine()        # <-- the ONLY change
    twin2, states2 = run_loop(cfg, swapped, n_iter=3)
    assert all(st.is_complete() for st in states2), "loop broke after engine swap"
    assert np.all(states2[-1].delta_dose == 0), "swapped engine not in effect"
    print(f"  swapped in {swapped.dose_engine.label!r} by setting ifaces.dose_engine only;")
    print(f"  loop ran {len(states2)} iterations, every state complete, no loop.py edit.")
    print("  FORGE-DROP-IN: PASS")

    _hr("CP3 check 3/3 — placeholder is isolated + labelled (no Forge physics)")
    import ast
    from matrix.interfaces import dose as dose_mod
    tree = ast.parse(open(dose_mod.__file__).read())
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            imported.add((node.module or "").split(".")[0])
    # the dose interface depends only on numpy + the package itself (no Forge, no scipy,
    # no transport/geometry libs) — the physics is exactly what Forge will later provide.
    allowed = {"numpy", "", "matrix", "config", "__future__"}  # config = own package module
    leaks = imported - allowed
    print(f"  dose module imports: {sorted(i for i in imported if i)}")
    assert not leaks, f"dose interface leaks a heavy dependency: {leaks}"
    assert eng.provisional and "NOT-Forge" in eng.label
    print(f"  purely analytic, zero physics deps; label={eng.label!r}")
    print(f"  Forge's real Monte-Carlo dose engine is deferred to 2027 (see ASSUMPTIONS.md).")
    print("  ISOLATED+LABELLED: PASS")

    _hr("CP3 GATE: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
