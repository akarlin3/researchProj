/**
 * CP2 pilot + CP1 sensitivity.
 *
 * Integrates the pilot subset (config.pilot.points × config.pilot.seeds) to the
 * FULL t_max with no early exit, so the whole R_incoh series is available to
 * re-label under every (T_v, recThresh) sensitivity setting from one integration.
 *
 * Reports:
 *   - t_abs censoring rate per point at the baseline t_max (the CP2 decision: is
 *     it ≲ censorTargetFrac? where is raising t_max impractical?).
 *   - how τ̂_abs (censored-exponential MLE) moves across T_v ∈ {60,120,240} and
 *     recThresh ∈ {0.75,0.80} — the CP1 robustness statement.
 *
 * Outputs (absorption_results/): pilot_summary.json, pilot_sensitivity.csv.
 * Run: node tools/absorption-recampaign/pilot.mjs [--seeds N] [--tmax T]
 */
import { readFileSync, writeFileSync, mkdirSync } from 'node:fs';
import { resolve, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { traceRun } from './tracer.mjs';
import { labelSeries } from './labeling.mjs';

const __dirname = dirname(fileURLToPath(import.meta.url));
const cfg = JSON.parse(
  readFileSync(resolve(__dirname, 'absorption.config.json'), 'utf8'),
);

function parseArgs(argv) {
  const a = {};
  for (let i = 2; i < argv.length; i++) {
    if (argv[i] === '--seeds') a.seeds = parseInt(argv[++i], 10);
    else if (argv[i] === '--tmax') a.tmax = parseFloat(argv[++i]);
  }
  return a;
}
const args = parseArgs(process.argv);

const { beta, dt, sampleStride } = cfg.model;
const t_max = args.tmax ?? cfg.model.t_max;
const theta = cfg.graze_criterion.theta;
const W = cfg.graze_criterion.W;
const recWin = cfg.absorption_criterion.recoveryWindowSec;
const seedsPerPoint = args.seeds ?? cfg.pilot.seeds;
const points = cfg.pilot.points;
const seed0 = { 0.5: 100000, 0.2: 200000 };

const T_v_set = cfg.sensitivity.T_v_set;
const rec_set = cfg.sensitivity.recoveryThreshold_set;

/** Censored-exponential MLE: τ̂ = Σ observed_time / #events; events=non-censored. */
function expMLE(pairs) {
  let total = 0;
  let d = 0;
  for (const [t, ev] of pairs) {
    total += t;
    if (ev) d++;
  }
  return { tau: d ? total / d : NaN, events: d, n: pairs.length, total };
}

const outDir = resolve(process.cwd(), cfg.output_dir);
mkdirSync(outDir, { recursive: true });

console.log(
  `Pilot: ${points.length} points × ${seedsPerPoint} seeds, t_max=${t_max}, no early-exit.`,
);
console.log(
  `Sensitivity grid: T_v∈${JSON.stringify(T_v_set)} × recThresh∈${JSON.stringify(rec_set)}\n`,
);

const perPoint = {};
const t0 = Date.now();

for (const { Np, A } of points) {
  const traces = [];
  for (let s = 0; s < seedsPerPoint; s++) {
    const seed = seed0[A] + s;
    const r = traceRun({
      Np,
      A,
      beta,
      seed,
      t_max,
      dt,
      sampleStride,
      label: {
        theta,
        W,
        T_v: cfg.absorption_criterion.T_v,
        recThresh: cfg.absorption_criterion.recoveryThreshold,
        recWin,
      },
      earlyExit: false,
      keepTrace: true,
    });
    traces.push({
      seed,
      R: r.R_incoh,
      sampleDt: r.sampleDt,
      t_graze: r.t_graze,
    });
  }
  // Re-label under every (T_v, recThresh).
  const grid = {};
  for (const T_v of T_v_set) {
    for (const recThresh of rec_set) {
      const pairs = traces.map((tr) => {
        const lab = labelSeries(
          tr.R,
          tr.sampleDt,
          { theta, W, T_v, recThresh, recWin },
          t_max,
        );
        return [lab.t_abs, lab.abs_censored ? 0 : 1, lab.n_grazes_before_abs];
      });
      const mle = expMLE(pairs.map((p) => [p[0], p[1]]));
      const censored = pairs.filter((p) => p[1] === 0).length;
      grid[`Tv${T_v}_rec${recThresh}`] = {
        T_v,
        recThresh,
        tau_abs: mle.tau,
        events: mle.events,
        censored,
        censored_frac: censored / pairs.length,
        median_grazes: median(pairs.map((p) => p[2])),
      };
    }
  }
  const key = `N${Np}_A${A}`;
  perPoint[key] = { Np, A, n: seedsPerPoint, grid };
  const base =
    grid[
      `Tv${cfg.absorption_criterion.T_v}_rec${cfg.absorption_criterion.recoveryThreshold}`
    ];
  console.log(
    `  N=${Np} A=${A}: baseline τ̂_abs=${fmt(base.tau_abs)} censored=${(base.censored_frac * 100).toFixed(0)}% ` +
      `(${base.censored}/${seedsPerPoint})  [${((Date.now() - t0) / 1000).toFixed(0)}s elapsed]`,
  );
}

function median(arr) {
  if (!arr.length) return NaN;
  const s = [...arr].sort((a, b) => a - b);
  const m = s.length;
  return m % 2 ? s[(m - 1) / 2] : (s[m / 2 - 1] + s[m / 2]) / 2;
}
function fmt(x) {
  return Number.isFinite(x) ? x.toFixed(1) : 'NA';
}

// --- CP1 sensitivity statement: τ̂_abs spread across the grid, per point ----- //
const sensitivityRows = [];
for (const key of Object.keys(perPoint)) {
  const { Np, A, grid } = perPoint[key];
  for (const gk of Object.keys(grid)) {
    const g = grid[gk];
    sensitivityRows.push([
      Np,
      A,
      g.T_v,
      g.recThresh,
      fmt(g.tau_abs),
      g.events,
      g.censored,
      (g.censored_frac * 100).toFixed(0),
      g.median_grazes,
    ]);
  }
}
const csvHeader =
  'N,A,T_v,recThresh,tau_abs,events,censored,censored_pct,median_grazes';
writeFileSync(
  resolve(outDir, 'pilot_sensitivity.csv'),
  csvHeader + '\n' + sensitivityRows.map((r) => r.join(',')).join('\n') + '\n',
);

// Per-point sensitivity: how much does τ̂_abs move across the grid (uncensored
// points only — where τ̂ is finite and meaningful)?
const robustness = {};
for (const key of Object.keys(perPoint)) {
  const { Np, A, grid } = perPoint[key];
  const taus = Object.values(grid)
    .map((g) => g.tau_abs)
    .filter((x) => Number.isFinite(x));
  const base =
    grid[
      `Tv${cfg.absorption_criterion.T_v}_rec${cfg.absorption_criterion.recoveryThreshold}`
    ];
  if (taus.length) {
    const lo = Math.min(...taus);
    const hi = Math.max(...taus);
    robustness[key] = {
      Np,
      A,
      baseline_tau_abs: base.tau_abs,
      baseline_censored_frac: base.censored_frac,
      tau_abs_min: lo,
      tau_abs_max: hi,
      spread_frac_of_baseline: Number.isFinite(base.tau_abs)
        ? (hi - lo) / base.tau_abs
        : null,
    };
  }
}

const summary = {
  config: {
    t_max,
    baseline_T_v: cfg.absorption_criterion.T_v,
    baseline_recThresh: cfg.absorption_criterion.recoveryThreshold,
    censorTargetFrac: cfg.pilot.censorTargetFrac,
    seedsPerPoint,
  },
  perPoint,
  robustness,
  t_max_decision: decideTmax(perPoint, cfg),
  wall_s: (Date.now() - t0) / 1000,
};
writeFileSync(
  resolve(outDir, 'pilot_summary.json'),
  JSON.stringify(summary, null, 2),
);

function decideTmax(perPoint, cfg) {
  const target = cfg.pilot.censorTargetFrac;
  const overs = [];
  for (const key of Object.keys(perPoint)) {
    const g =
      perPoint[key].grid[
        `Tv${cfg.absorption_criterion.T_v}_rec${cfg.absorption_criterion.recoveryThreshold}`
      ];
    if (g.censored_frac > target)
      overs.push({ point: key, censored_frac: g.censored_frac });
  }
  return {
    target,
    over_target_points: overs,
    note:
      overs.length === 0
        ? `All pilot points ≤ ${(target * 100).toFixed(0)}% t_abs censoring at t_max=${cfg.model.t_max}; keep t_max.`
        : `Points exceeding ${(target * 100).toFixed(0)}% censoring: ${overs
            .map((o) => `${o.point}=${(o.censored_frac * 100).toFixed(0)}%`)
            .join(
              ', ',
            )}. If these are stable-after-graze (churn that never absorbs), raising t_max will NOT reduce censoring — report as a finding.`,
  };
}

console.log(`\nPilot complete in ${summary.wall_s.toFixed(0)}s.`);
console.log('t_max decision:', summary.t_max_decision.note);
console.log(
  'Wrote absorption_results/pilot_summary.json, pilot_sensitivity.csv',
);
