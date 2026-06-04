/**
 * CP4 — timestep sanity check.
 *
 * At N=8 and N=32 (primary A=0.5, β=0.05), rerun 40 seeds at dt/4 and compare
 * the lifetime distributions against the shipped dt=0.05. Emits one JSONL row
 * per run (tagged with dt) to campaign_results/cp4_dt.jsonl; the rigorous
 * comparison (KM overlay + two-sample test + τ CIs + PASS/FAIL verdict) is done
 * by the Python analysis (analysis.py --cp4), which reads this file.
 *
 * The shipped dt is part of what the paper characterizes, so this checks that
 * the measured collapse statistics are not an artifact of the integration step.
 *
 * Run:  node tools/chimera-campaign/cp4_dt.mjs
 */
import { appendFileSync, openSync, closeSync, mkdirSync } from 'node:fs';
import { execSync } from 'node:child_process';
import { resolve } from 'node:path';
import { runUntilCollapse } from './runner.mjs';
import { RUNNER_VERSION } from './integrator.mjs';

const A = 0.5;
const BETA = 0.05;
const THETA = 0.85;
const W = 5.0;
const T_MAX = 2000;
const SEED0 = 300000;
const SEEDS = 40;
const NS = [8, 32];
const SHIPPED_DT = 0.05;
const DTS = [SHIPPED_DT, SHIPPED_DT / 4]; // 0.05 and 0.0125

const gitHash = (() => {
  try {
    return execSync('git rev-parse --short HEAD', { encoding: 'utf8' }).trim();
  } catch {
    return 'unknown';
  }
})();

mkdirSync('campaign_results', { recursive: true });
const outPath = resolve(process.cwd(), 'campaign_results/cp4_dt.jsonl');
const fd = openSync(outPath, 'w'); // fresh each run

console.log('CP4 — timestep sanity check (dt vs dt/4)');
console.log(
  `A=${A}, β=${BETA}, θ=${THETA}, W=${W}, t_max=${T_MAX}, ${SEEDS} seeds (seed0=${SEED0})`,
);
console.log('');
console.log('  N    dt        collapsed   τ̂(MLE,censored)   mean(uncens)');

for (const Np of NS) {
  for (const dt of DTS) {
    const lifes = [];
    let censored = 0;
    let totalObs = 0;
    let observed = 0;
    for (let s = 0; s < SEEDS; s++) {
      const seed = SEED0 + s;
      const r = runUntilCollapse({
        Np,
        A,
        beta: BETA,
        seed,
        t_max: T_MAX,
        dt,
        theta: THETA,
        W,
      });
      appendFileSync(
        fd,
        JSON.stringify({
          N: Np,
          A,
          beta: BETA,
          dt,
          seed,
          lifetime: r.lifetime,
          censored: r.censored,
          theta: THETA,
          W,
          t_max: T_MAX,
          git_hash: gitHash,
          runner_version: RUNNER_VERSION,
        }) + '\n',
      );
      totalObs += r.lifetime; // lifetime carries t_max when censored
      if (r.censored) censored++;
      else {
        lifes.push(r.lifetime);
        observed++;
      }
    }
    // Exponential MLE with right-censoring: τ̂ = total observed time / # collapses.
    const tauMLE = observed > 0 ? totalObs / observed : NaN;
    const meanUncens = lifes.length
      ? lifes.reduce((a, b) => a + b, 0) / lifes.length
      : NaN;
    console.log(
      `  ${String(Np).padStart(2)}   ${dt.toFixed(4)}    ${String(SEEDS - censored).padStart(2)}/${SEEDS}        ` +
        `${tauMLE.toFixed(1).padStart(7)}s          ${meanUncens.toFixed(1)}s`,
    );
  }
}
closeSync(fd);
console.log('');
console.log(
  `Wrote ${outPath}. Run the Python analysis for the KM overlay + two-sample test + verdict.`,
);
