/**
 * Determinism gate — the absorption re-campaign's t_graze must reproduce the
 * published campaign's lifetime bit-for-bit (same seeds, same shipped RK4, same
 * min(R₁,R₂)>θ sustained-for-W criterion at sampleStride=0.1). A single mismatch
 * beyond fp noise fails the gate.
 *
 * The published campaign was logged at t_max=2000; we compare every (A,N,seed)
 * present in BOTH files. Runs that were graze-censored in the original at
 * t_max=2000 are bit-identical only if the re-run used the same t_max — they are
 * compared on (lifetime, censored) jointly. (Runs re-run at a larger t_max would
 * legitimately resolve a few additional grazes; that is reported, not failed.)
 *
 * Output: absorption_results/determinism_gate.json
 * Run: node tools/absorption-recampaign/determinism.mjs
 */
import { readFileSync, writeFileSync, mkdirSync } from 'node:fs';
import { resolve, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const __dirname = dirname(fileURLToPath(import.meta.url));
const cfg = JSON.parse(
  readFileSync(resolve(__dirname, 'absorption.config.json'), 'utf8'),
);
const ROOT = process.cwd();

const orig = readJsonl(resolve(ROOT, cfg.inputs.campaign_jsonl));
const rerun = readJsonl(resolve(ROOT, cfg.output));

function readJsonl(p) {
  return readFileSync(p, 'utf8')
    .split('\n')
    .filter((l) => l.trim())
    .map((l) => JSON.parse(l));
}
function key(r) {
  return `${r.A}|${r.beta}|${r.N}|${r.seed}`;
}

const origByKey = new Map(orig.map((r) => [key(r), r]));

let compared = 0;
let mismatches = 0;
let worstDev = 0;
let worstCase = null;
let tmaxMismatchCensoring = 0; // runs where t_max differs and censoring legitimately changed
const perPoint = new Map();

for (const r of rerun) {
  const o = origByKey.get(key(r));
  if (!o) continue;
  const pk = `N${r.N}_A${r.A}`;
  if (!perPoint.has(pk))
    perPoint.set(pk, { N: r.N, A: r.A, compared: 0, mismatches: 0, worst: 0 });
  const pp = perPoint.get(pk);

  // Compare t_graze (re-run) to lifetime (original).
  const sameTmax = r.t_max === o.t_max;
  const bothUncensored = !r.graze_censored && !o.censored;

  if (bothUncensored || sameTmax) {
    const dev = Math.abs(r.t_graze - o.lifetime);
    const censoredMatch = !!r.graze_censored === !!o.censored;
    compared++;
    pp.compared++;
    if (dev > worstDev) {
      worstDev = dev;
      worstCase = {
        key: key(r),
        t_graze: r.t_graze,
        lifetime: o.lifetime,
        re_cens: r.graze_censored,
        orig_cens: o.censored,
      };
    }
    if (dev > pp.worst) pp.worst = dev;
    // Fail if uncensored-vs-uncensored deviate, or censoring flags disagree at
    // the same t_max.
    if ((bothUncensored && dev > 1e-9) || (sameTmax && !censoredMatch)) {
      mismatches++;
      pp.mismatches++;
    }
  } else {
    // Different t_max AND at least one censored: censoring may legitimately
    // change. Only count if the original was uncensored (then t_graze must still
    // match regardless of t_max).
    if (!o.censored) {
      const dev = Math.abs(r.t_graze - o.lifetime);
      compared++;
      pp.compared++;
      if (dev > worstDev) {
        worstDev = dev;
        worstCase = { key: key(r), t_graze: r.t_graze, lifetime: o.lifetime };
      }
      if (dev > 1e-9) {
        mismatches++;
        pp.mismatches++;
      }
    } else {
      tmaxMismatchCensoring++;
    }
  }
}

const gate = {
  passed: mismatches === 0,
  compared,
  mismatches,
  worstAbsDev: worstDev,
  worstCase,
  tmaxMismatchCensoring,
  rerun_t_max: rerun.length ? rerun[0].t_max : null,
  orig_t_max: orig.length ? orig[0].t_max : null,
  perPoint: [...perPoint.values()].sort((a, b) => a.A - b.A || a.N - b.N),
  _note:
    'Gate passes iff every (A,N,seed) shared with the published campaign reproduced its logged lifetime as t_graze to within 1e-9 s (and censoring flags agree at equal t_max). Same shipped RK4, same min(R1,R2)>θ sustained-for-W criterion, same stride.',
};

mkdirSync(resolve(ROOT, cfg.output_dir), { recursive: true });
writeFileSync(
  resolve(ROOT, cfg.output_dir, 'determinism_gate.json'),
  JSON.stringify(gate, null, 2),
);

console.log(
  `Determinism gate: ${gate.passed ? 'PASSED ✅' : 'FAILED ❌'} — compared ${compared}, ` +
    `mismatches ${mismatches}, worst |Δ| ${worstDev.toExponential(2)} s`,
);
if (worstCase) console.log('worst case:', JSON.stringify(worstCase));
console.log('Wrote absorption_results/determinism_gate.json');
