"""CP0 gate -- the identifiability object is built AND the kill test has a definite verdict.

Five independent checks, all green => CP0 gate passes:

  (1) WIRING        levy-core imports via _paths (read-only); Ouroboros reuse path resolves;
                    confirm Ouroboros has NO Fisher/CRLB symbol (Levy's layer is net-new).
  (2) UNIT TESTS    the levy-core pytest suite passes (forward Jacobian = finite-diff, Rician
                    info factor -> Gaussian limit, CRLB monotone, MLE recovers truth, etc.).
  (3) CRLB LIMIT    Rician FIM -> Gaussian FIM as SNR -> inf (the bound is correctly built).
  (4) WALL OBJECT   the recovery-collapse wall exists and is characterised at the headline
                    cell: a finite wall SNR* with a bootstrap CI, alpha-D degeneracy reported,
                    and it lands INSIDE the realistic SNR band.
  (5) VERDICT       the pre-registered REFUTE has a definite, honest answer (wall_stands xor
                    refuted), and the scope boundary (research_recovers) is reported.

Run:  <proteus python> verify_cp0.py            # FAST (smoke; default)
      <proteus python> verify_cp0.py --full     # full-N bootstrap
Exit 0 = CP0 gate green; nonzero = a check failed (the script names which).
"""
from __future__ import annotations

import os
import subprocess
import sys
import time

import numpy as np

import _paths  # noqa: E402  (same directory)

FULL = "--full" in sys.argv


def _hr(title: str) -> None:
    print("\n" + "=" * 78 + f"\n{title}\n" + "=" * 78)


def check_wiring() -> None:
    _hr("CP0 check 1/5 -- dependency wiring (read-only) + net-new confirmation")
    paths = _paths.add_all()
    for k, v in paths.items():
        print(f"  {k:<12} -> {v}")
    import levy  # noqa: F401
    from levy import fisher, forward, glreuse, identifiability, noise, wall  # noqa: F401
    print(f"  import levy        OK  ({levy.__file__})")
    # Reused GL layer (clean copy of Ouroboros operators) works without Ouroboros' heavy deps.
    assert abs(glreuse.gl_weights(0.5, 3)[0] - 1.0) < 1e-12
    print("  reused GL layer (levy.glreuse) OK -- clean copy, no pysindy dependency.")
    # Net-new confirmation: grep the Ouroboros SOURCE (do not import it -- it needs pysindy)
    # for any Fisher/CRLB/Cramer/FIM symbol. Absence => Levy's layer is net-new (GATE A).
    import re
    pat = re.compile(r"fisher|crlb|cram[ée]r|\bfim\b|information.matrix|profile.likelihood",
                     re.IGNORECASE)
    hits = []
    for f in sorted(_paths.OUROBOROS.glob("*.py")):
        for i, line in enumerate(f.read_text(errors="ignore").splitlines(), 1):
            if pat.search(line):
                hits.append(f"{f.name}:{i}: {line.strip()[:70]}")
    assert not hits, "unexpected Fisher/CRLB code in Ouroboros source:\n  " + "\n  ".join(hits)
    print(f"  grep of Ouroboros source ({len(list(_paths.OUROBOROS.glob('*.py')))} .py files): "
          "NO Fisher/CRLB/Cramer/FIM symbol -> Levy's CRLB layer is net-new.")
    print("  WIRING: PASS")


def check_unit_tests() -> None:
    _hr("CP0 check 2/5 -- levy-core unit test suite")
    core = _paths.LEVY_CORE
    t0 = time.time()
    r = subprocess.run([sys.executable, "-m", "pytest", "-q"], cwd=str(core),
                       capture_output=True, text=True)
    out = r.stdout + r.stderr
    passed = r.returncode == 0
    print(f"  pytest -> {'PASS' if passed else 'FAIL'} (rc={r.returncode}, {time.time()-t0:.1f}s)")
    if not passed:
        print("\n".join(out.splitlines()[-20:]))
        raise SystemExit("UNIT TESTS FAILED")
    print("  " + out.strip().splitlines()[-1])
    print("  UNIT TESTS: PASS")


def check_crlb_limit() -> None:
    _hr("CP0 check 3/5 -- Rician CRLB -> Gaussian CRLB at high SNR (bound correctly built)")
    from levy import fisher, forward, wall
    b = wall.default_b_design()
    theta = forward.StretchedExp().theta
    rels = {}
    for snr in (20, 500):
        fg = fisher.fisher_matrix(b, theta, snr, "gaussian")
        fr = fisher.fisher_matrix(b, theta, snr, "rician")
        rels[snr] = np.linalg.norm(fg - fr) / np.linalg.norm(fg)
        print(f"  SNR={snr:4d}: ||FIM_G - FIM_R||/||FIM_G|| = {rels[snr]:.4f}")
    assert rels[500] < 1e-2, "Rician FIM does not approach Gaussian at high SNR"
    assert rels[20] > rels[500], "Rician/Gaussian gap should shrink with SNR"
    print("  CRLB LIMIT: PASS")


def check_wall_object(rep) -> None:
    _hr("CP0 check 4/5 -- recovery-collapse wall exists and is characterised")
    from levy import wall
    h = rep.headline
    lo, hi = wall.REALISTIC_SNR_BAND
    print(f"  analytic CRLB wall SNR*  = {h.wall_snr:.1f}")
    print(f"  empirical wall SNR*      = {h.wall_snr_emp:.1f}  95% CI [{h.wall_ci[0]:.1f}, {h.wall_ci[1]:.1f}]")
    assert np.isfinite(h.wall_snr), "no analytic wall located"
    assert h.wall_exists, "wall_exists is False at the headline cell"
    assert lo <= h.wall_snr <= hi, f"analytic wall {h.wall_snr:.1f} not in realistic band [{lo},{hi}]"
    assert np.isfinite(h.wall_ci[0]) and np.isfinite(h.wall_ci[1]), "no bootstrap CI on the wall"
    assert h.wall_ci[0] < h.wall_ci[1], "degenerate wall CI"
    # alpha-D degeneracy must be reported (a number, finite)
    from levy import fisher
    rho = fisher.crlb(wall.default_b_design(b_max=3000.0, n_b=4),
                      np.array([1.0, 1.5e-3, 0.85]), 30, "rician").rho_alpha_D
    print(f"  alpha-D degeneracy (n_b=4,b_max=3000,SNR=30): rho_aD = {rho:+.3f}")
    assert np.isfinite(rho)
    print("  WALL OBJECT: PASS")


def check_verdict(rep) -> None:
    _hr("CP0 check 5/5 -- pre-registered REFUTE has a definite, honest verdict")
    print(f"  wall_stands={rep.wall_stands}  refuted={rep.refuted}  "
          f"research_recovers={rep.research_recovers}")
    # exactly one of wall_stands / refuted is true (definite verdict)
    assert rep.wall_stands != rep.refuted, "verdict is not definite (wall_stands == refuted)"
    if rep.refuted:
        print("  VERDICT: REFUTE TRIGGERED -- wedge dead. (Reported honestly; CP0 gate still PASSES")
        print("           because the kill test ran and produced a definite answer.)")
    else:
        print("  VERDICT: WALL STANDS (scoped to realistic clinical few-b acquisition).")
        print(f"           Scope boundary reported: research-dense acquisition recovers alpha "
              f"(research_recovers={rep.research_recovers}).")
    print("  VERDICT: PASS")


def main() -> int:
    print("CP0 verification -- fractional-order identifiability object + kill test")
    check_wiring()
    check_unit_tests()
    check_crlb_limit()

    from levy import seeding, wall
    n_boot = 300 if FULL else 80
    t0 = time.time()
    rep = wall.cp0_verdict(rng=seeding.make_rng(), n_boot=n_boot, do_ci=True)
    print(f"\n  (cp0_verdict computed in {time.time()-t0:.1f}s, n_boot={n_boot})")

    check_wall_object(rep)
    check_verdict(rep)

    _hr("CP0 GATE: PASS")
    print("  Fisher/CRLB identifiability object built (net-new); kill test ran with a definite,")
    print("  scoped verdict; load-bearing wall carries a bootstrap CI. See results/RESULTS_CP0.md.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
