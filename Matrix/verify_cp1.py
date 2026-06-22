"""CP1 gate — synthetic twin is reproducible AND the harness runs end-to-end (stubs).

Two falsifiable checks, both green => CP1 passes:

  (1) TWIN REPRODUCIBLE  two twins built from the same seed are bit-identical in every
                         ground-truth field; a different seed differs.
  (2) HARNESS END-TO-END  the full four-stage loop runs with pass-through stubs for all
                          n_iter iterations and every LoopState is complete (every stage
                          wrote its output, nothing raised).

Run:  <proteus python> Matrix/verify_cp1.py        # exit 0 = green
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from matrix import MatrixConfig, Twin, Interfaces, run_loop


def _hr(t): print("\n" + "=" * 74 + f"\n{t}\n" + "=" * 74)


def check_twin_reproducible() -> None:
    _hr("CP1 check 1/2 — synthetic twin reproducible from seed")
    cfg = MatrixConfig()
    a = Twin.build(cfg)
    b = Twin.build(cfg)
    for fld in ("D", "Dstar", "f", "labels", "highdstar", "lowsnr", "snr_map", "dose"):
        assert np.array_equal(getattr(a, fld), getattr(b, fld)), f"{fld} not reproducible"
    c = Twin.build(cfg.with_seed(cfg.seed + 1))
    assert not np.array_equal(a.f, c.f), "different seed gave identical twin"
    print(f"  V={a.n_voxels} voxels | tumour={int((a.labels==1).sum())} "
          f"OAR={int((a.labels==2).sum())} normal={int((a.labels==0).sum())} "
          f"| low-SNR(untrust) zone={int(a.lowsnr.sum())} | tumour∩low-SNR="
          f"{int(((a.labels==1)&a.lowsnr).sum())}")
    print(f"  same seed -> bit-identical (8/8 fields); seed+1 -> differs.  REPRODUCIBLE: PASS")


def check_harness_end_to_end() -> None:
    _hr("CP1 check 2/2 — harness runs end-to-end with pass-through stubs")
    cfg = MatrixConfig()
    twin, states = run_loop(cfg, Interfaces.passthrough())
    assert len(states) == cfg.n_iter, "wrong number of iterations"
    for s in states:
        assert s.is_complete(), f"iteration {s.iteration}: a stage left state incomplete"
    prov = states[0].components
    print(f"  ran {len(states)} iterations; every LoopState complete (all 4 stages wrote).")
    print(f"  stub provenance: ruler={prov['ruler']!r}")
    print(f"                   trust={prov['trust_gate']!r}")
    print(f"                   action={prov['action_gate']!r}")
    print(f"                   dose={prov['dose_engine']!r}")
    print("  HARNESS END-TO-END: PASS")


def main() -> int:
    print("CP1 verification — Matrix synthetic twin + harness skeleton")
    check_twin_reproducible()
    check_harness_end_to_end()
    _hr("CP1 GATE: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
