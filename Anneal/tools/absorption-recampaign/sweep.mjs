/**
 * CP2 — absorption-grade re-campaign driver.
 *
 * Re-runs the full collapse-time campaign (same logged seeds as the published
 * campaign) with the absorption-grade tracer, producing ONE JSONL row per run
 * carrying BOTH labels:
 *
 *   { N, A, beta, dt, seed, t_graze, graze_censored, t_abs, abs_censored,
 *     n_grazes_before_abs, T_b, n_breath_peaks, breath_cycles, theta, W, T_v,
 *     recThresh, recWin, t_max, git_hash, runner_version, wall_ms }
 *
 * Worker pool (one per core), resumable (skips (params,seed) already present),
 * append-only. Early-exits absorbed runs at t_abs + T_v; censored runs (no
 * absorption by t_max) integrate the full horizon — that is intentional and is
 * what makes the t_abs censoring rate an honest measurement.
 *
 * Usage:
 *   node tools/absorption-recampaign/sweep.mjs                 # full campaign
 *   node tools/absorption-recampaign/sweep.mjs --seeds 20      # cap seeds/point
 *   node tools/absorption-recampaign/sweep.mjs --workers 4
 *   node tools/absorption-recampaign/sweep.mjs --out path.jsonl
 *   node tools/absorption-recampaign/sweep.mjs --tmax 4000     # override t_max
 */
import {
  Worker,
  isMainThread,
  parentPort,
  workerData,
} from 'node:worker_threads';
import { fileURLToPath } from 'node:url';
import { execSync } from 'node:child_process';
import {
  appendFileSync,
  readFileSync,
  existsSync,
  mkdirSync,
  openSync,
  closeSync,
} from 'node:fs';
import { dirname, resolve } from 'node:path';
import { availableParallelism } from 'node:os';
import { traceRun } from './tracer.mjs';
import { RUNNER_VERSION } from '../chimera-campaign/integrator.mjs';

const __dirname = dirname(fileURLToPath(import.meta.url));

function labelOf(cfg) {
  return {
    theta: cfg.graze_criterion.theta,
    W: cfg.graze_criterion.W,
    T_v: cfg.absorption_criterion.T_v,
    recThresh: cfg.absorption_criterion.recoveryThreshold,
    recWin: cfg.absorption_criterion.recoveryWindowSec,
  };
}

// --------------------------------------------------------------------------- //
// Worker
// --------------------------------------------------------------------------- //
if (!isMainThread) {
  const { label, breath } = workerData;
  parentPort.on('message', (job) => {
    if (job === null) process.exit(0);
    const t0 = Date.now();
    const r = traceRun({
      Np: job.Np,
      A: job.A,
      beta: job.beta,
      seed: job.seed,
      t_max: job.t_max,
      dt: job.dt,
      sampleStride: job.sampleStride,
      label,
      breath,
      earlyExit: true,
    });
    parentPort.postMessage({
      row: {
        N: job.Np,
        A: job.A,
        beta: job.beta,
        dt: job.dt,
        seed: job.seed,
        t_graze: r.t_graze,
        graze_censored: r.graze_censored,
        t_abs: r.t_abs,
        abs_censored: r.abs_censored,
        n_grazes_before_abs: r.n_grazes_before_abs,
        T_b: r.T_b,
        n_breath_peaks: r.n_breath_peaks,
        breath_cycles: r.breath_cycles,
        theta: label.theta,
        W: label.W,
        T_v: label.T_v,
        recThresh: label.recThresh,
        recWin: label.recWin,
        t_max: job.t_max,
        git_hash: workerData.gitHash,
        runner_version: RUNNER_VERSION,
        wall_ms: Date.now() - t0,
      },
      steps: r.stepsIntegrated,
    });
  });
}

// --------------------------------------------------------------------------- //
// Main
// --------------------------------------------------------------------------- //
if (isMainThread) await main();

function parseArgs(argv) {
  const a = { config: resolve(__dirname, 'absorption.config.json') };
  for (let i = 2; i < argv.length; i++) {
    const k = argv[i];
    if (k === '--seeds') a.seedsCap = parseInt(argv[++i], 10);
    else if (k === '--workers') a.workers = parseInt(argv[++i], 10);
    else if (k === '--config') a.config = resolve(argv[++i]);
    else if (k === '--out') a.out = resolve(argv[++i]);
    else if (k === '--tmax') a.tmax = parseFloat(argv[++i]);
  }
  return a;
}

function gitHash() {
  try {
    return execSync('git rev-parse --short HEAD', { encoding: 'utf8' }).trim();
  } catch {
    return 'unknown';
  }
}

function jobKey(j) {
  return `${j.A}|${j.beta}|${j.Np}|${j.seed}|${j.t_max}`;
}

function loadDoneKeys(outPath) {
  const done = new Set();
  if (!existsSync(outPath)) return done;
  for (const line of readFileSync(outPath, 'utf8').split('\n')) {
    if (!line.trim()) continue;
    try {
      const r = JSON.parse(line);
      done.add(`${r.A}|${r.beta}|${r.N}|${r.seed}|${r.t_max}`);
    } catch {
      /* tolerate a partial trailing line */
    }
  }
  return done;
}

function buildJobs(cfg, seedsCap, tmaxOverride) {
  const { beta, dt, sampleStride } = cfg.model;
  const t_max = tmaxOverride ?? cfg.model.t_max;
  const jobs = [];
  for (const sweep of cfg.sweeps) {
    const seeds = seedsCap ? Math.min(seedsCap, sweep.seeds) : sweep.seeds;
    for (const Np of sweep.Ns) {
      for (let s = 0; s < seeds; s++) {
        jobs.push({
          Np,
          A: sweep.A,
          beta,
          dt,
          sampleStride,
          t_max,
          seed: sweep.seed0 + s,
        });
      }
    }
  }
  return jobs;
}

function fmtDur(ms) {
  const s = ms / 1000;
  if (s < 90) return `${s.toFixed(0)}s`;
  const m = s / 60;
  if (m < 90) return `${m.toFixed(1)}min`;
  return `${(m / 60).toFixed(1)}h`;
}

async function main() {
  const args = parseArgs(process.argv);
  const cfg = JSON.parse(readFileSync(args.config, 'utf8'));
  const outPath = args.out ? args.out : resolve(process.cwd(), cfg.output);
  mkdirSync(dirname(outPath), { recursive: true });

  const nWorkers = args.workers ?? Math.max(1, availableParallelism());
  const gh = gitHash();
  const label = labelOf(cfg);
  const breath = cfg.breath;

  const allJobs = buildJobs(cfg, args.seedsCap, args.tmax);
  const done = loadDoneKeys(outPath);
  const jobs = allJobs.filter((j) => !done.has(jobKey(j)));

  console.log(`Absorption re-campaign: ${cfg.name} v${cfg.version}`);
  console.log(`Output: ${outPath} (git ${gh}, ${RUNNER_VERSION})`);
  console.log(
    `Workers: ${nWorkers} | total: ${allJobs.length} | done: ${allJobs.length - jobs.length} | to run: ${jobs.length}`,
  );
  console.log(
    `Label: θ=${label.theta} W=${label.W} | absorption: T_v=${label.T_v} recThresh=${label.recThresh} recWin=${label.recWin} | t_max=${args.tmax ?? cfg.model.t_max}`,
  );
  if (jobs.length === 0) {
    console.log('Nothing to do — campaign already complete for this config.');
    return;
  }

  const fd = openSync(outPath, 'a');
  const t0 = Date.now();
  let nextJob = 0;
  let completed = 0;
  let stepsTotal = 0;
  const total = jobs.length;
  const logEvery = Math.max(1, Math.floor(total / 50));

  await new Promise((resolveAll) => {
    let alive = Math.min(nWorkers, total);
    const assign = (w) => {
      if (nextJob >= total) {
        w.postMessage(null);
        return;
      }
      w.postMessage(jobs[nextJob++]);
    };
    for (let i = 0; i < Math.min(nWorkers, total); i++) {
      const w = new Worker(fileURLToPath(import.meta.url), {
        workerData: { gitHash: gh, label, breath },
      });
      w.on('message', (msg) => {
        appendFileSync(fd, JSON.stringify(msg.row) + '\n');
        completed++;
        stepsTotal += msg.steps;
        if (completed % logEvery === 0 || completed === total) {
          const elapsed = Date.now() - t0;
          const rate = completed / elapsed;
          const eta = (total - completed) / rate;
          process.stdout.write(
            `\r  ${completed}/${total} (${((completed / total) * 100).toFixed(0)}%)  ` +
              `elapsed ${fmtDur(elapsed)}  ETA ${fmtDur(eta)}  ${(rate * 1000).toFixed(1)} runs/s    `,
          );
        }
        assign(w);
      });
      w.on('exit', () => {
        if (--alive === 0) resolveAll();
      });
      w.on('error', (e) => console.error('\nWorker error:', e));
      assign(w);
    }
  });

  closeSync(fd);
  const elapsed = Date.now() - t0;
  process.stdout.write('\n');
  console.log(
    `Done: ${completed} runs in ${fmtDur(elapsed)} on ${nWorkers} workers ` +
      `(${((completed / elapsed) * 1000).toFixed(1)} runs/s; ${(stepsTotal / 1e6).toFixed(1)}M RK4 steps).`,
  );
}
