/**
 * CP1 validation gate.
 *
 * Reproduce a known basin-map point from the existing characterization
 * (docs/CHIMERA_CHARACTERIZATION.md §3, the CP2 N-sweep at A=0.5, β=0.05):
 *
 *     N=8  → basin ≈ 63%      N=32 → basin ≈ 92%      (24 seeds, seed0=5000)
 *
 * This uses the OLD collapse definition (the probe's: a seed is "in basin" when
 * fracLive ≥ 0.8, where alive ⇔ max(R₁,R₂) > 0.9 and min(R₁,R₂) < 0.85), since
 * that is what produced those numbers. It re-implements the probe's exact loop
 * — transient 12 s discarded, then alive evaluated EVERY step over a 28 s window
 * — on top of the headless integrator's primitives, so this is an
 * apples-to-apples reproduction, not a re-derivation.
 *
 * Run:  node tools/chimera-campaign/validate.mjs
 */
import {
  mulberry32,
  seedChimera,
  orderParam,
  makeScratch,
  rk4StepInPlace,
  SYNC_HI,
  INCOH_LO,
} from './integrator.mjs';

const DT = 0.05;
const TRANSIENT = 12; // s, discarded
const ANALYZE = 28; // s, alive evaluated every step
const IN_BASIN = 0.8; // fracLive threshold for "in basin"
const SEED0 = 5000;
const SEEDS = 24;

/** fracLive of one seed: fraction of the analyze window that is a live chimera. */
function fracLive(Np, A, beta, seed) {
  const mu = (1 + A) / 2;
  const nu = (1 - A) / 2;
  const alpha = Math.PI / 2 - beta;
  const phases = seedChimera(Np, mulberry32(seed));
  const scratch = makeScratch(2 * Np);
  const nTrans = Math.round(TRANSIENT / DT);
  const nAna = Math.round(ANALYZE / DT);
  for (let s = 0; s < nTrans; s++)
    rk4StepInPlace(phases, Np, mu, nu, alpha, DT, scratch);
  let live = 0;
  for (let s = 0; s < nAna; s++) {
    rk4StepInPlace(phases, Np, mu, nu, alpha, DT, scratch);
    const hi = Math.max(
      orderParam(phases, 0, Np).R,
      orderParam(phases, Np, Np).R,
    );
    const lo = Math.min(
      orderParam(phases, 0, Np).R,
      orderParam(phases, Np, Np).R,
    );
    if (hi > SYNC_HI && lo < INCOH_LO) live++;
  }
  return live / nAna;
}

function basinFrac(Np, A, beta) {
  let inBasin = 0;
  for (let s = 0; s < SEEDS; s++) {
    if (fracLive(Np, A, beta, SEED0 + s) >= IN_BASIN) inBasin++;
  }
  return { inBasin, frac: inBasin / SEEDS };
}

// Binomial 95% tolerance band (normal approx) for n=24 around an expected p.
function band(p, n = SEEDS) {
  const sd = Math.sqrt((p * (1 - p)) / n);
  return [Math.max(0, p - 1.96 * sd), Math.min(1, p + 1.96 * sd)];
}

const A = 0.5;
const beta = 0.05;
const targets = [
  { Np: 8, expected: 0.63 },
  { Np: 32, expected: 0.92 },
];

console.log(
  'CP1 VALIDATION GATE — basin reproduction (old collapse definition)',
);
console.log(
  `A=${A}, β=${beta}, ${SEEDS} seeds (seed0=${SEED0}), transient ${TRANSIENT}s, analyze ${ANALYZE}s, dt=${DT}`,
);
console.log('');
console.log('  N    measured   expected   95% binom band     verdict');
let allPass = true;
for (const { Np, expected } of targets) {
  const { inBasin, frac } = basinFrac(Np, A, beta);
  const [lo, hi] = band(expected);
  const pass = frac >= lo && frac <= hi;
  if (!pass) allPass = false;
  console.log(
    `  ${String(Np).padStart(2)}   ${(frac * 100).toFixed(1).padStart(6)}%   ` +
      `${(expected * 100).toFixed(0).padStart(6)}%   ` +
      `[${(lo * 100).toFixed(0)}%, ${(hi * 100).toFixed(0)}%]`.padStart(14) +
      `   ${pass ? 'PASS' : 'FAIL'}  (${inBasin}/${SEEDS})`,
  );
}
console.log('');
console.log(allPass ? 'GATE: PASS ✅' : 'GATE: FAIL ❌ — diagnose before CP2');
process.exit(allPass ? 0 : 1);
