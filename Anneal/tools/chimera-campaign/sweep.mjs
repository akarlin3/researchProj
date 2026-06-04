/**
 * CP3 — campaign sweep driver.
 *
 * Runs the full collapse-time campaign defined by campaign.config.json across a
 * worker pool (one worker per core), appends one JSONL row per run, and is
 * RESUMABLE — on restart it skips (params, seed) combinations already present in
 * the output file. Prints progress with an ETA and, at the end, the estimated
 * total compute time for the full campaign on this machine.
 *
 * Each output row:
 *   { N, A, beta, dt, seed, lifetime, censored, theta, W, t_max,
 *     git_hash, runner_version, wall_ms }
 * (lifetime carries t_max when censored; censored is the boolean flag.)
 *
 * Usage:
 *   node tools/chimera-campaign/sweep.mjs                 # full campaign
 *   node tools/chimera-campaign/sweep.mjs --seeds 20      # cap seeds/point (smoke)
 *   node tools/chimera-campaign/sweep.mjs --workers 8
 *   node tools/chimera-campaign/sweep.mjs --config path/to/config.json
 *
 * Deterministic: every row is reproducible from (seed, params, dt, theta, W).
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
import { runUntilCollapse } from './runner.mjs';
import { RUNNER_VERSION } from './integrator.mjs';

const __dirname = dirname(fileURLToPath(import.meta.url));

// ---------------------------------------------------------------------------
// Worker: run one job and post the resulting row back.
// ---------------------------------------------------------------------------
if (!isMainThread) {
  parentPort.on('message', (job) => {
    if (job === null) {
      process.exit(0);
    }
    const t0 = Date.now();
    const r = runUntilCollapse({
      Np: job.Np,
      A: job.A,
      beta: job.beta,
      seed: job.seed,
      t_max: job.t_max,
      dt: job.dt,
      theta: job.theta,
      W: job.W,
      sampleStride: job.sampleStride,
    });
    parentPort.postMessage({
      row: {
        N: job.Np,
        A: job.A,
        beta: job.beta,
        dt: job.dt,
        seed: job.seed,
        lifetime: r.lifetime,
        censored: r.censored,
        theta: job.theta,
        W: job.W,
        t_max: job.t_max,
        git_hash: workerData.gitHash,
        runner_version: RUNNER_VERSION,
        wall_ms: Date.now() - t0,
      },
      stepsIntegrated: r.stepsIntegrated,
    });
  });
}

// ---------------------------------------------------------------------------
// Main driver.
// ---------------------------------------------------------------------------
if (isMainThread) {
  await main();
}

function parseArgs(argv) {
  const a = { config: resolve(__dirname, 'campaign.config.json') };
  for (let i = 2; i < argv.length; i++) {
    const k = argv[i];
    if (k === '--seeds') a.seedsCap = parseInt(argv[++i], 10);
    else if (k === '--workers') a.workers = parseInt(argv[++i], 10);
    else if (k === '--config') a.config = resolve(argv[++i]);
    else if (k === '--out') a.out = resolve(argv[++i]);
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
  // Identity of a run for resumability (params + seed + criterion + dt).
  return `${j.A}|${j.beta}|${j.Np}|${j.seed}|${j.theta}|${j.W}|${j.dt}|${j.t_max}`;
}

function loadDoneKeys(outPath) {
  const done = new Set();
  if (!existsSync(outPath)) return done;
  const text = readFileSync(outPath, 'utf8');
  for (const line of text.split('\n')) {
    if (!line.trim()) continue;
    try {
      const r = JSON.parse(line);
      done.add(
        `${r.A}|${r.beta}|${r.N}|${r.seed}|${r.theta}|${r.W}|${r.dt}|${r.t_max}`,
      );
    } catch {
      /* tolerate a partially-written trailing line */
    }
  }
  return done;
}

function buildJobs(cfg, seedsCap) {
  const { beta, dt, sampleStride, t_max } = cfg.model;
  const { theta, W } = cfg.criterion;
  const jobs = [];
  for (const sweep of cfg.sweeps) {
    const seeds = seedsCap ? Math.min(seedsCap, sweep.seeds) : sweep.seeds;
    for (const Np of sweep.Ns) {
      for (let s = 0; s < seeds; s++) {
        jobs.push({
          label: sweep.label,
          Np,
          A: sweep.A,
          beta,
          dt,
          sampleStride,
          t_max,
          theta,
          W,
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

  const allJobs = buildJobs(cfg, args.seedsCap);
  const done = loadDoneKeys(outPath);
  const jobs = allJobs.filter((j) => !done.has(jobKey(j)));

  console.log(`Campaign: ${cfg.name} v${cfg.version}`);
  console.log(`Config: ${args.config}`);
  console.log(`Output: ${outPath}  (git ${gh}, ${RUNNER_VERSION})`);
  console.log(
    `Workers: ${nWorkers} | total runs: ${allJobs.length} | already done: ${done.size ? allJobs.length - jobs.length : 0} | to run: ${jobs.length}`,
  );
  if (args.seedsCap) console.log(`(seeds capped at ${args.seedsCap}/point)`);
  if (jobs.length === 0) {
    console.log('Nothing to do — campaign already complete for this config.');
    return;
  }

  // Open the output once for append; flush each row immediately (resumable even
  // if killed mid-run).
  const fd = openSync(outPath, 'a');

  const t0 = Date.now();
  let nextJob = 0;
  let completed = 0;
  let stepsTotal = 0;
  const total = jobs.length;
  const logEvery = Math.max(1, Math.floor(total / 50));

  const workers = [];
  await new Promise((resolveAll) => {
    let alive = nWorkers;

    const assign = (w) => {
      if (nextJob >= total) {
        w.postMessage(null); // tell worker to exit
        return;
      }
      w.postMessage(jobs[nextJob++]);
    };

    for (let i = 0; i < Math.min(nWorkers, total); i++) {
      const w = new Worker(fileURLToPath(import.meta.url), {
        workerData: { gitHash: gh },
      });
      workers.push(w);
      w.on('message', (msg) => {
        appendFileSync(fd, JSON.stringify(msg.row) + '\n');
        completed++;
        stepsTotal += msg.stepsIntegrated;
        if (completed % logEvery === 0 || completed === total) {
          const elapsed = Date.now() - t0;
          const rate = completed / elapsed; // jobs/ms
          const eta = (total - completed) / rate;
          process.stdout.write(
            `\r  ${completed}/${total} (${((completed / total) * 100).toFixed(0)}%)  ` +
              `elapsed ${fmtDur(elapsed)}  ETA ${fmtDur(eta)}  ` +
              `${(rate * 1000).toFixed(1)} runs/s    `,
          );
        }
        assign(w);
      });
      w.on('exit', () => {
        alive--;
        if (alive === 0) resolveAll();
      });
      w.on('error', (e) => {
        console.error('\nWorker error:', e);
      });
      assign(w); // kick off the worker's first job
    }
  });

  closeSync(fd);
  const elapsed = Date.now() - t0;
  process.stdout.write('\n');
  console.log('');
  console.log(
    `Done: ${completed} runs in ${fmtDur(elapsed)} on ${nWorkers} workers.`,
  );
  console.log(
    `Throughput: ${((completed / elapsed) * 1000).toFixed(1)} runs/s; ${(stepsTotal / 1e6).toFixed(1)}M RK4 steps.`,
  );

  // Full-campaign compute estimate (extrapolate from this run's measured rate to
  // the FULL job set, in case --seeds capped it).
  const fullRuns = buildJobs(cfg, undefined).length;
  const perRun = elapsed / completed;
  const fullEst = perRun * fullRuns;
  console.log('');
  console.log(`Full-campaign estimate on this machine (${nWorkers} cores):`);
  console.log(`  full job count: ${fullRuns} runs`);
  console.log(
    `  measured per-run wall: ${perRun.toFixed(1)} ms (amortized across ${nWorkers} workers)`,
  );
  console.log(`  ⇒ estimated total compute: ${fmtDur(fullEst)}`);
}
