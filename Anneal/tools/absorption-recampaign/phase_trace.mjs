/**
 * CP3 phase-clustering subset — re-traces the 100 lowest-id seeds per point
 * (N ∈ {8,16,32,64} × A ∈ {0.5,0.2}, matching PR #41) recording R_incoh(t)
 * through and past absorption, so the Python analysis can locate the breath peak
 * preceding each TRUE absorption and run the Rayleigh test (PR #41 machinery).
 *
 * Storage is bounded: a run that absorbs early-exits T_v after its absorption
 * crossing (R_incoh spans [0, t_abs+T_v]); only absorbed runs carry the R_incoh
 * array (censored runs have no absorption to phase-locate, so they store labels
 * only). The JS per-run T_b is included so the Python estimator can be
 * cross-checked against it.
 *
 * Output: absorption_results/phase_traces.jsonl
 * Run: node tools/absorption-recampaign/phase_trace.mjs
 */
import { readFileSync, writeFileSync, mkdirSync } from 'node:fs';
import { resolve, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { traceRun } from './tracer.mjs';

const __dirname = dirname(fileURLToPath(import.meta.url));
const cfg = JSON.parse(
  readFileSync(resolve(__dirname, 'absorption.config.json'), 'utf8'),
);
const ROOT = process.cwd();

const { beta, dt, sampleStride, t_max } = cfg.model;
const label = {
  theta: cfg.graze_criterion.theta,
  W: cfg.graze_criterion.W,
  T_v: cfg.absorption_criterion.T_v,
  recThresh: cfg.absorption_criterion.recoveryThreshold,
  recWin: cfg.absorption_criterion.recoveryWindowSec,
};
const breath = cfg.breath;
const points = [
  { Np: 8, A: 0.5 },
  { Np: 16, A: 0.5 },
  { Np: 32, A: 0.5 },
  { Np: 64, A: 0.5 },
  { Np: 8, A: 0.2 },
  { Np: 16, A: 0.2 },
  { Np: 32, A: 0.2 },
  { Np: 64, A: 0.2 },
];
const nSeeds = cfg.phase_subset.nSeedsPerPoint;
const seed0 = { 0.5: 100000, 0.2: 200000 };

const outDir = resolve(ROOT, cfg.output_dir);
mkdirSync(outDir, { recursive: true });
const outPath = resolve(outDir, 'phase_traces.jsonl');

const round4 = (x) => Number(x.toFixed(4));
const lines = [];
let totalAbsorbed = 0;
let totalCensored = 0;

for (const { Np, A } of points) {
  let absorbed = 0;
  let censored = 0;
  for (let s = 0; s < nSeeds; s++) {
    const seed = seed0[A] + s;
    const r = traceRun({
      Np,
      A,
      beta,
      seed,
      t_max,
      dt,
      sampleStride,
      label,
      breath,
      earlyExit: true,
      keepTrace: true,
      keepRsync: false,
    });
    const row = {
      N: Np,
      A,
      seed,
      t_graze: r.t_graze,
      graze_censored: r.graze_censored,
      t_abs: r.t_abs,
      abs_censored: r.abs_censored,
      n_grazes_before_abs: r.n_grazes_before_abs,
      T_b: r.T_b,
      n_breath_peaks: r.n_breath_peaks,
      breath_cycles: r.breath_cycles,
      sampleDt: round4(r.sampleDt),
      absIndex: r.absIndex,
      grazeIndex: r.grazeIndex,
    };
    if (!r.abs_censored) {
      row.n = r.R_incoh.length;
      row.R_incoh = Array.from(r.R_incoh, round4);
      absorbed++;
    } else {
      censored++;
    }
    lines.push(JSON.stringify(row));
  }
  totalAbsorbed += absorbed;
  totalCensored += censored;
  console.log(
    `phase-trace: N=${Np} A=${A} — absorbed ${absorbed}, censored ${censored} (traces stored for absorbed)`,
  );
}

writeFileSync(outPath, lines.join('\n') + '\n');
console.log(
  `\nWrote ${outPath} — ${lines.length} rows (${totalAbsorbed} absorbed w/ traces, ${totalCensored} censored labels-only).`,
);
