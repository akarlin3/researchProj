/**
 * CP4 — supervisor over-trigger measurement (NO behavior change).
 *
 * Replays the SHIPPED supervisor's collapse detector — alive ⇔ max(R₁,R₂) > 0.9
 * ∧ min(R₁,R₂) < 0.85 (chimera.ts isChimeraAlive / SYNC_HI / INCOH_LO); a
 * "firing" = not-alive sustained ≥ COLLAPSE_HOLD_S = 2.0 s (chimeraVoice.ts) —
 * over the NATURAL (non-re-perturbed) traced trajectories, and classifies every
 * firing by whether the event SELF-RECOVERS under the absorption criterion
 * (recovery = R_incoh drops below recThresh sustained ≥ recWin AFTER the
 * episode). A firing on a self-recovering event is an over-trigger: the shipped
 * supervisor would have re-perturbed a graze that was going to reform on its own.
 *
 * The shipped supervisor and its constants are imported but UNCHANGED; we only
 * read its detector thresholds. We do not re-seed (that would change the
 * trajectory) — this is a pure measurement over the same dynamics the campaign
 * integrates.
 *
 * Output: absorption_results/supervisor_overtrigger.{json,csv}
 * Run: node tools/absorption-recampaign/supervisor_replay.mjs [--seeds N]
 */
import { readFileSync, writeFileSync, mkdirSync } from 'node:fs';
import { resolve, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { traceRun } from './tracer.mjs';
import { SYNC_HI, INCOH_LO } from '../chimera-campaign/integrator.mjs';

// Shipped supervisor hold (chimeraVoice.ts COLLAPSE_HOLD_S); imported as a
// constant value to keep this replay faithful without touching the voice.
const COLLAPSE_HOLD_S = 2.0;

const __dirname = dirname(fileURLToPath(import.meta.url));
const cfg = JSON.parse(
  readFileSync(resolve(__dirname, 'absorption.config.json'), 'utf8'),
);
const ROOT = process.cwd();

function parseArgs(argv) {
  const a = {};
  for (let i = 2; i < argv.length; i++)
    if (argv[i] === '--seeds') a.seeds = parseInt(argv[++i], 10);
  return a;
}
const args = parseArgs(process.argv);

const { beta, dt, sampleStride, t_max } = cfg.model;
const label = {
  theta: cfg.graze_criterion.theta,
  W: cfg.graze_criterion.W,
  T_v: cfg.absorption_criterion.T_v,
  recThresh: cfg.absorption_criterion.recoveryThreshold,
  recWin: cfg.absorption_criterion.recoveryWindowSec,
};
const recThresh = label.recThresh;
const recWin = label.recWin;
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
const nSeeds = args.seeds ?? cfg.phase_subset.nSeedsPerPoint;
const seed0 = { 0.5: 100000, 0.2: 200000 };

/** True if R_incoh has a run of ≥ need consecutive samples < thresh from `start`. */
function sustainedBelow(R, thresh, need, start) {
  let run = 0;
  for (let i = start; i < R.length; i++) {
    if (R[i] < thresh) {
      if (++run >= need) return true;
    } else run = 0;
  }
  return false;
}

/**
 * Replay the supervisor detector over one (R_incoh, R_sync) trace. Returns the
 * firing episodes (not-alive sustained ≥ hold), each classified self-recovering
 * (an over-trigger) or terminal.
 */
function replay(Rincoh, Rsync, sampleDt) {
  const need = Math.max(1, Math.round(COLLAPSE_HOLD_S / sampleDt));
  const recNeed = Math.max(1, Math.round(recWin / sampleDt));
  const n = Rincoh.length;
  // Maximal not-alive episodes.
  const episodes = [];
  let notAlive = 0;
  let start = -1;
  for (let i = 0; i < n; i++) {
    const alive = Rsync[i] > SYNC_HI && Rincoh[i] < INCOH_LO;
    if (!alive) {
      if (notAlive === 0) start = i;
      notAlive++;
    } else {
      if (notAlive > 0) episodes.push({ start, end: i - 1, dur: notAlive });
      notAlive = 0;
      start = -1;
    }
  }
  if (notAlive > 0)
    episodes.push({ start, end: n - 1, dur: notAlive, openEnd: true });

  let firings = 0;
  let overTriggers = 0;
  let terminal = 0;
  const wasted = [];
  for (const e of episodes) {
    if (e.dur < need) continue; // supervisor would not have fired
    firings++;
    // Self-recover under the absorption criterion: recovery AFTER this episode.
    const recovers = sustainedBelow(Rincoh, recThresh, recNeed, e.end + 1);
    if (recovers) {
      overTriggers++;
      // Time the chimera stayed not-alive past the 2 s hold before reforming —
      // the window the re-perturbation pre-empted unnecessarily.
      wasted.push((e.dur - need) * sampleDt);
    } else {
      terminal++;
    }
  }
  return { firings, overTriggers, terminal, wasted };
}

function median(a) {
  if (!a.length) return NaN;
  const s = [...a].sort((x, y) => x - y);
  const m = s.length;
  return m % 2 ? s[(m - 1) / 2] : (s[m / 2 - 1] + s[m / 2]) / 2;
}

const outDir = resolve(ROOT, cfg.output_dir);
mkdirSync(outDir, { recursive: true });

console.log(
  `CP4 supervisor replay: detector alive⇔maxR>${SYNC_HI}∧minR<${INCOH_LO}, fire=not-alive≥${COLLAPSE_HOLD_S}s; ` +
    `${points.length} points × ${nSeeds} seeds, full t_max=${t_max}.\n`,
);

const rows = [];
const t0 = Date.now();
for (const { Np, A } of points) {
  let firings = 0;
  let over = 0;
  let terminal = 0;
  let runHours = 0;
  const wastedAll = [];
  let runsWithFire = 0;
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
      earlyExit: false, // measure the full natural trajectory (incl. A=0.2 churn)
      keepTrace: true,
      keepRsync: true,
    });
    const rep = replay(r.R_incoh, r.R_sync, r.sampleDt);
    firings += rep.firings;
    over += rep.overTriggers;
    terminal += rep.terminal;
    wastedAll.push(...rep.wasted);
    runHours += (r.nSamples * r.sampleDt) / 3600;
    if (rep.firings > 0) runsWithFire++;
  }
  const row = {
    N: Np,
    A,
    n_runs: nSeeds,
    firings,
    over_triggers: over,
    terminal_fires: terminal,
    over_trigger_rate: firings ? over / firings : NaN,
    median_wasted_s: median(wastedAll),
    firings_per_run: firings / nSeeds,
    firings_per_runhour: runHours > 0 ? firings / runHours : NaN,
    runs_with_any_fire_frac: runsWithFire / nSeeds,
  };
  rows.push(row);
  console.log(
    `  N=${Np} A=${A}: firings=${firings} over-trigger=${(row.over_trigger_rate * 100).toFixed(0)}% ` +
      `(${over}/${firings}) median_wasted=${Number.isFinite(row.median_wasted_s) ? row.median_wasted_s.toFixed(1) : 'NA'}s ` +
      `firings/run=${row.firings_per_run.toFixed(1)}  [${((Date.now() - t0) / 1000).toFixed(0)}s]`,
  );
}

writeFileSync(
  resolve(outDir, 'supervisor_overtrigger.json'),
  JSON.stringify(
    {
      detector: {
        syncHi: SYNC_HI,
        incohLo: INCOH_LO,
        hold_s: COLLAPSE_HOLD_S,
        recThresh,
        recWin,
      },
      t_max,
      nSeeds,
      rows,
      wall_s: (Date.now() - t0) / 1000,
    },
    null,
    2,
  ),
);
const csvHeader =
  'N,A,n_runs,firings,over_triggers,terminal_fires,over_trigger_rate,median_wasted_s,firings_per_run,firings_per_runhour,runs_with_any_fire_frac';
writeFileSync(
  resolve(outDir, 'supervisor_overtrigger.csv'),
  csvHeader +
    '\n' +
    rows
      .map((r) =>
        [
          r.N,
          r.A,
          r.n_runs,
          r.firings,
          r.over_triggers,
          r.terminal_fires,
          r.over_trigger_rate.toFixed(4),
          Number.isFinite(r.median_wasted_s)
            ? r.median_wasted_s.toFixed(2)
            : 'NA',
          r.firings_per_run.toFixed(3),
          Number.isFinite(r.firings_per_runhour)
            ? r.firings_per_runhour.toFixed(2)
            : 'NA',
          r.runs_with_any_fire_frac.toFixed(3),
        ].join(','),
      )
      .join('\n') +
    '\n',
);
console.log(
  `\nWrote absorption_results/supervisor_overtrigger.{json,csv} in ${((Date.now() - t0) / 1000).toFixed(0)}s.`,
);
