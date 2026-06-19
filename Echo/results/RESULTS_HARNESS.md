# CP1 method self-test (SOLID -- no upstream dependency)

Seed 0, n=20000, level=0.1 (90% intervals).

| check | value | target | pass |
|---|---|---|---|
| scaled interval recovers analytic repeat-coverage | 0.7585 | 0.7552 | True |
| **bias invariance** (precision-not-accuracy) | 0.7585 | 0.7552 | True |
| scale=0.5 coverage tracks analytic | 0.4355 | 0.4391 | True |
| scale=0.75 coverage tracks analytic | 0.6180 | 0.6170 | True |
| scale=1.0 coverage tracks analytic | 0.7627 | 0.7552 | True |
| scale=1.5 coverage tracks analytic | 0.9186 | 0.9190 | True |
| scale=2.0 coverage tracks analytic | 0.9795 | 0.9800 | True |
| distinct from Gauge: Spearman before/after rescale | 0.3728 / 0.3728 | equal | True |
| distinct from Gauge: coverage before/after rescale | 0.7591 / 0.4435 | differ | True |

**ALL_PASS: True**

Reading: a perfectly measurement-scaled 90% interval is *expected* to show ~76% test-retest coverage, not 90% -- this is the derivable gap between accuracy-coverage and repeat-coverage, and the reason Echo's signal cannot be read as accuracy. Bias invariance is the precision-not-accuracy guarantee. The Spearman/coverage split is the distinctness-from-Gauge proof.