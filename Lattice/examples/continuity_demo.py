"""Print the continuity-gate table: every family reduces to bi-exponential.

Run: python examples/continuity_demo.py
"""

from lattice.cohort import continuity_residual

EXACT = {"biexp", "dispersion_lognormal", "stretched", "triexp"}


def main() -> None:
    print("Continuity residuals  max|family@limit - biexp|  (n=2000, identical ground truth)")
    print("-" * 72)
    for fam in ["biexp", "dispersion_lognormal", "stretched", "triexp", "dispersion_gamma"]:
        r = continuity_residual(fam, n=2000)
        kind = "exact" if fam in EXACT else "asymptotic (k=1e8)"
        flag = "OK" if (r == 0.0 if fam in EXACT else r < 1e-6) else "FAIL"
        print(f"  {fam:22s} {r:.3e}   {kind:20s} [{flag}]")


if __name__ == "__main__":
    main()
