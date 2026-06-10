/**
 * Mechanism probe 2b — off-manifold Watanabe–Strogatz constants as the carrier
 * of the N-independent 3.2x prolongation (decomposition experiment).
 *
 * For identical oscillators the finite-N two-population dynamics IS the WS
 * dynamics with the seed's constants; the three-variable reduced flow is exact
 * on the uniform-constants (Poisson) submanifold. So compare, per campaign seed:
 *
 *   ACTUAL  — the seed's exact phases (seedChimera + mulberry32, verbatim).
 *   PARTNER — a "constants-uniformized" IC: per population, a Möbius push of
 *             the exact splay grid psi_j = 2*pi*j/Np (uniform constants — on the
 *             Poisson submanifold to the finite-N sampling floor) matched to the
 *             ACTUAL IC's complex Z_1. Construction: push along the fixed
 *             direction a = r*e^{i*pi/Np} (never antipode-aligned with a grid
 *             point for even Np, so |Z1| -> 1 as r -> 1), bisect r so
 *             |Z1(pushed grid)| = |Z1_target|, then apply the global rotation
 *             delta = arg(Z1_target) - arg(Z1(pushed grid)) to all phases of
 *             the block (a rotated splay is equally uniform-constants, and the
 *             closure defect |Z_m - Z_1^m| is rotation-invariant). This matches
 *             Z_1 exactly up to bisection resolution (gate: |dZ1| <= 1e-10).
 *             Same collective state, manifold-projected constants — the
 *             paired-IC pattern of cp4_intervention.mjs with the harmonic
 *             injection replaced by a manifold projection. (A plain 2-D Newton
 *             solve on `a` over the FIXED canonical grid fails for ~16% of
 *             sync populations: for even Np a grid point can sit at the exact
 *             antipode of `a`, capping the reachable |Z1| along that direction
 *             below the ~0.999 target — the rotation construction removes the
 *             problem exactly.)
 *
 * Both are integrated with the absorption campaign's machinery: identical RK4
 * core, identical streaming two-timescale Labeler (theta=0.85, W=5, T_v=120,
 * recThresh=0.8, recWin=5; dt=0.05, sampleStride=0.1, t_max=2000), so the
 * ACTUAL runs must reproduce t_abs_meas (reduced_runs.jsonl) bit-for-bit —
 * that is gate B. Gate A: the ACTUAL IC's (Rsync0, Rincoh0, dphi0) must
 * reproduce the recorded row values. Gate C: partner |dZ1| <= 1e-10 and
 * partner D0 at/below the sampling floor (poissonDistance, M=4).
 *
 * Per run also: T_b + capture phase in the measured-data convention
 * (breath.mjs detectPeaks on the smoothed pre-absorption window; Tb = median
 * peak spacing; phase = 2*pi*frac((t_abs - last peak)/Tb)) and
 * n_cycles = t_abs / T_b (the k_cyc variable).
 *
 * Output: manifold_results/ws_decomposition_N{N}.jsonl (one row per run, kinds
 * 'actual' and 'partner').
 *
 * Run:  node tools/manifold-probe/ws_decomposition.mjs --N 8 [--nseeds 200]
 */
import { readFileSync, writeFileSync, mkdirSync } from 'node:fs';
import { resolve, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import {
  mulberry32,
  seedChimera,
  orderParam,
  makeScratch,
  rk4StepInPlace,
  TAU,
} from '../chimera-campaign/integrator.mjs';
import { makeLabeler } from '../absorption-recampaign/labeling.mjs';
import { detectPeaks } from '../absorption-recampaign/breath.mjs';
import { circularMoments, poissonDistance, mobiusBlock, DEFAULT_M } from './moments.mjs';

const __dirname = dirname(fileURLToPath(import.meta.url));
const ROOT = resolve(__dirname, '../..');

// campaign / absorption settings (absorption.config.json, verbatim)
const A = 0.5;
const BETA = 0.05;
const DT = 0.05;
const SAMPLE_STRIDE = 0.1;
const T_MAX = 2000;
const M = DEFAULT_M; // 4
const MIN_CYCLES = 2;

const args = process.argv.slice(2);
const getArg = (k, dflt) => {
  const i = args.indexOf(k);
  return i >= 0 ? args[i + 1] : dflt;
};
const N = parseInt(getArg('--N', '8'), 10);
const NSEEDS = parseInt(getArg('--nseeds', '200'), 10);

const wrapPi = (d) => {
  let x = ((d % TAU) + TAU) % TAU;
  if (x > Math.PI) x -= TAU;
  return x;
};

// --------------------------------------------------------------------------
// Integrate from a supplied phase vector with the absorption labeler online.
// Verbatim adaptation of absorption-recampaign/tracer.mjs traceRun: the ONLY
// change is phases0 is supplied by the caller instead of seedChimera (gate B
// verifies the actual-IC path reproduces the recorded t_abs exactly).
// --------------------------------------------------------------------------
function traceFromPhases(phases0, Np) {
  const mu = (1 + A) / 2;
  const nu = (1 - A) / 2;
  const alpha = Math.PI / 2 - BETA;
  const n = 2 * Np;
  const phases = Float64Array.from(phases0);
  const scratch = makeScratch(n);

  const nSteps = Math.round(T_MAX / DT);
  const sampleEvery = Math.max(1, Math.round(SAMPLE_STRIDE / DT));
  const sampleDt = sampleEvery * DT;

  const labeler = makeLabeler(sampleDt); // labeling.mjs DEFAULTS (campaign)
  const Rincoh = [];

  const sampleAt = () => {
    const r1 = orderParam(phases, 0, Np).R;
    const r2 = orderParam(phases, Np, Np).R;
    const lo = r1 < r2 ? r1 : r2;
    Rincoh.push(lo);
    labeler.push(lo);
  };

  sampleAt();
  for (let step = 1; step <= nSteps; step++) {
    rk4StepInPlace(phases, Np, mu, nu, alpha, DT, scratch);
    if (step % sampleEvery === 0) {
      sampleAt();
      if (labeler.absorptionConfirmed) break;
    }
  }
  const lab = labeler.result(T_MAX);

  // measured-data convention: peaks on smoothed pre-absorption window
  const absIdx = lab.abs_censored
    ? Rincoh.length
    : Math.round(lab.t_abs / sampleDt);
  const pre = Rincoh.slice(0, Math.min(absIdx, Rincoh.length));
  let Tb = null;
  let phase = null;
  let nPeaks = 0;
  if (pre.length >= Math.round(6.0 / sampleDt)) {
    const { peaks } = detectPeaks(pre, sampleDt);
    nPeaks = peaks.length;
    if (peaks.length - 1 >= MIN_CYCLES) {
      const ptimes = peaks.map((p) => p * sampleDt);
      const iv = [];
      for (let i = 1; i < ptimes.length; i++) iv.push(ptimes[i] - ptimes[i - 1]);
      iv.sort((a, b) => a - b);
      const m = iv.length;
      Tb = m % 2 ? iv[(m - 1) / 2] : (iv[m / 2 - 1] + iv[m / 2]) / 2;
      if (!lab.abs_censored) {
        const frac = (((lab.t_abs - ptimes[ptimes.length - 1]) / Tb) % 1 + 1) % 1;
        phase = TAU * frac;
      }
    }
  }
  return {
    t_graze: lab.t_graze,
    graze_censored: lab.graze_censored,
    t_abs: lab.t_abs,
    abs_censored: lab.abs_censored,
    n_grazes_before_abs: lab.n_grazes_before_abs,
    T_b: Tb,
    n_breath_peaks: nPeaks,
    n_cycles: Tb ? lab.t_abs / Tb : null,
    capture_phase: phase,
  };
}

// --------------------------------------------------------------------------
// Möbius match: build a splay-constants phase block whose Z1 equals the target.
// Push direction fixed at phi0 = pi/Np (no grid point at its antipode for even
// Np, so |Z1| rises to 1 as r -> 1); bisect r for |Z1| = |target|; rotate the
// whole block to fix arg(Z1) exactly.
// --------------------------------------------------------------------------
function pushedGrid(r, Np) {
  const phi0 = Math.PI / Np;
  const out = new Float64Array(Np);
  mobiusBlock(out, 0, Np, { re: r * Math.cos(phi0), im: r * Math.sin(phi0) });
  return out;
}

function matchSplayBlock(z1t, Np) {
  const Rt = Math.hypot(z1t.re, z1t.im);
  let lo = 0;
  let hi = 1 - 1e-16;
  for (let k = 0; k < 200; k++) {
    const mid = (lo + hi) / 2;
    const z = circularMoments(pushedGrid(mid, Np), 0, Np, 1)[0];
    if (Math.hypot(z.re, z.im) < Rt) lo = mid;
    else hi = mid;
  }
  const r = (lo + hi) / 2;
  const block = pushedGrid(r, Np);
  const z = circularMoments(block, 0, Np, 1)[0];
  const delta = Math.atan2(z1t.im, z1t.re) - Math.atan2(z.im, z.re);
  for (let j = 0; j < Np; j++) block[j] = ((block[j] + delta) % TAU + TAU) % TAU;
  const zf = circularMoments(block, 0, Np, 1)[0];
  return {
    block,
    r,
    delta,
    err: Math.hypot(zf.re - z1t.re, zf.im - z1t.im),
  };
}

// --------------------------------------------------------------------------
// main
// --------------------------------------------------------------------------
const rrPath = resolve(ROOT, 'reduced_results/reduced_runs.jsonl');
const rows = readFileSync(rrPath, 'utf8')
  .split('\n')
  .filter((l) => l.trim())
  .map((l) => JSON.parse(l))
  .filter((r) => r.N === N && r.A === A)
  .sort((a, b) => a.seed - b.seed)
  .slice(0, NSEEDS);

console.log(`ws_decomposition: N=${N}, ${rows.length} seeds`);
const t0 = Date.now();
const out = [];

for (const row of rows) {
  const seed = row.seed;
  const phases = seedChimera(N, mulberry32(seed));

  // gate A: collective IC must reproduce the recorded row
  const op1 = orderParam(phases, 0, N);
  const op2 = orderParam(phases, N, N);
  const Rsync0 = Math.max(op1.R, op2.R);
  const Rincoh0 = Math.min(op1.R, op2.R);
  const dphi0 = wrapPi(op1.Phi - op2.Phi);
  const icDev = Math.max(
    Math.abs(Rsync0 - row.Rsync0),
    Math.abs(Rincoh0 - row.Rincoh0),
    Math.abs(dphi0 - row.dphi0),
  );

  // actual per-population Z1 and D0
  const z1p1 = circularMoments(phases, 0, N, 1)[0];
  const z1p2 = circularMoments(phases, N, N, 1)[0];
  const dAct1 = poissonDistance(phases, 0, N, M).D;
  const dAct2 = poissonDistance(phases, N, N, M).D;

  // partner: Möbius push of the splay grid matched to each Z1
  const s1 = matchSplayBlock(z1p1, N);
  const s2 = matchSplayBlock(z1p2, N);
  const partner = new Float64Array(2 * N);
  partner.set(s1.block, 0);
  partner.set(s2.block, N);
  const dPar1 = poissonDistance(partner, 0, N, M).D;
  const dPar2 = poissonDistance(partner, N, N, M).D;

  const ra = traceFromPhases(phases, N);
  const rp = traceFromPhases(partner, N);

  out.push(
    JSON.stringify({
      N,
      seed,
      kind: 'actual',
      ...ra,
      d0_pop1: dAct1,
      d0_pop2: dAct2,
      ic_dev: icDev,
      t_abs_meas: row.t_abs_meas,
      t_abs_match: ra.t_abs === row.t_abs_meas,
      t_capture_reduced: row.t_capture,
    }),
  );
  out.push(
    JSON.stringify({
      N,
      seed,
      kind: 'partner',
      ...rp,
      d0_pop1: dPar1,
      d0_pop2: dPar2,
      z1_err_pop1: s1.err,
      z1_err_pop2: s2.err,
      t_capture_reduced: row.t_capture,
    }),
  );
}

const outDir = resolve(ROOT, 'manifold_results');
mkdirSync(outDir, { recursive: true });
const outPath = resolve(outDir, `ws_decomposition_N${N}.jsonl`);
writeFileSync(outPath, out.join('\n') + '\n');
console.log(
  `Wrote ${outPath} (${out.length} rows, ${((Date.now() - t0) / 1000).toFixed(1)}s)`,
);
