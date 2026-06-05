/**
 * CP1 — traced re-run sample for the breath-phase / absorption study.
 *
 * For each configured (N, A) point we take the 100 lowest-id campaign seeds present
 * for that point (= the "first 100 logged seeds"; seed0+0..99), skip any
 * right-censored run, and RE-INTEGRATE each surviving seed with the manifold tracer
 * at the campaign stride, recording R_incoh(t) = min(R1,R2) until the collapse
 * criterion fires plus a 60 s tail (for the CP3 absorption check).
 *
 * Determinism gate: every re-run reuses the campaign seed and the shipped-identical
 * RK4 + the same min(R1,R2)>theta sustained-for-W criterion at sampleStride=0.1, so
 * the traced lifetime must equal the LOGGED lifetime bit-for-bit. We assert this and
 * report the max absolute deviation; a single mismatch fails the gate.
 *
 * Output (new files only; campaign JSONL is read-only input):
 *   phase_results/cp1_traces.jsonl    — one compact row per traced run:
 *       { N, A, seed, lifetime, loggedLifetime, censored, collapseIndex,
 *         sampleDt, n, R_incoh:[...] }   (t[i] = i*sampleDt, omitted to save space)
 *   phase_results/cp1_determinism.json — gate verdict + per-point match counts.
 *
 * Run:  node tools/breath-phase/cp1_trace.mjs
 */
import { readFileSync, writeFileSync, mkdirSync } from 'node:fs';
import { resolve, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { traceFromSeed } from '../manifold-probe/trace.mjs';
import { DEFAULT_M } from '../manifold-probe/moments.mjs';

const __dirname = dirname(fileURLToPath(import.meta.url));
const cfg = JSON.parse(
  readFileSync(resolve(__dirname, 'phase.config.json'), 'utf8'),
);
const M = cfg.cp1.M ?? DEFAULT_M;
const { beta, dt, sampleStride, t_max, theta, W } = cfg.model;
const { points, nSeedsPerPoint, tailAfterCollapse } = cfg.cp1;

const campaignPath = resolve(process.cwd(), cfg.inputs.campaign_jsonl);
const rows = readFileSync(campaignPath, 'utf8')
  .split('\n')
  .filter((l) => l.trim())
  .map((l) => JSON.parse(l));

const outDir = resolve(process.cwd(), cfg.output_dir);
mkdirSync(outDir, { recursive: true });
const outPath = resolve(outDir, 'cp1_traces.jsonl');
const gatePath = resolve(outDir, 'cp1_determinism.json');

/** Lowest-`k` campaign seeds for (N,A), with their logged lifetime/censored. */
function pickLogged(N, A, k) {
  return rows
    .filter((r) => r.N === N && r.A === A)
    .sort((a, b) => a.seed - b.seed)
    .slice(0, k);
}

const round4 = (x) => Number(x.toFixed(4));
const lines = [];
const perPoint = [];
let worstAbsDev = 0;
let worstCase = null;
let totalTraced = 0;
let totalSkippedCensored = 0;

for (const { Np, A } of points) {
  const logged = pickLogged(Np, A, nSeedsPerPoint);
  let traced = 0;
  let skipped = 0;
  let mismatches = 0;
  let pointWorst = 0;
  for (const lr of logged) {
    if (lr.censored) {
      skipped++;
      continue;
    }
    const tr = traceFromSeed({
      Np,
      A,
      beta,
      seed: lr.seed,
      t_max,
      dt,
      theta,
      W,
      M,
      sampleStride,
      tailAfterCollapse,
    });
    // Determinism gate: traced lifetime vs logged lifetime.
    const dev = Math.abs(tr.lifetime - lr.lifetime);
    if (dev > pointWorst) pointWorst = dev;
    if (dev > worstAbsDev) {
      worstAbsDev = dev;
      worstCase = {
        N: Np,
        A,
        seed: lr.seed,
        traced: tr.lifetime,
        logged: lr.lifetime,
      };
    }
    // exact (same RK4, same criterion, same stride) — tolerate only fp noise.
    if (dev > 1e-9) mismatches++;
    lines.push(
      JSON.stringify({
        N: Np,
        A,
        seed: lr.seed,
        lifetime: tr.lifetime,
        loggedLifetime: lr.lifetime,
        censored: tr.censored,
        collapseIndex: tr.collapseIndex,
        sampleDt: round4(tr.sampleEvery * tr.dt),
        n: tr.R_incoh.length,
        R_incoh: Array.from(tr.R_incoh, round4),
      }),
    );
    traced++;
    totalTraced++;
  }
  totalSkippedCensored += skipped;
  perPoint.push({
    N: Np,
    A,
    traced,
    skippedCensored: skipped,
    mismatches,
    worstAbsDev: pointWorst,
  });
  console.log(
    `CP1: N=${Np} A=${A} — traced ${traced}, skipped(censored) ${skipped}, mismatches ${mismatches}, worstDev ${pointWorst.toExponential(2)}`,
  );
}

writeFileSync(outPath, lines.join('\n') + '\n');

const gate = {
  passed: worstAbsDev <= 1e-9,
  worstAbsDev,
  worstCase,
  totalTraced,
  totalSkippedCensored,
  sampleStride,
  theta,
  W,
  perPoint,
  _note:
    'Gate passes iff every traced re-run reproduced its logged campaign lifetime to within 1e-9 s. Same seed, same shipped RK4, same min(R1,R2)>theta sustained-for-W criterion at sampleStride=0.1.',
};
writeFileSync(gatePath, JSON.stringify(gate, null, 2));
console.log(
  `\nWrote ${outPath} (${lines.length} traces)\nWrote ${gatePath} — gate ${gate.passed ? 'PASSED' : 'FAILED'} (worst |dev| ${worstAbsDev.toExponential(2)} s)`,
);
