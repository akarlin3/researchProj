/**
 * CP0/CP1/CP2 data extraction for the deterministic-transient tests.
 *
 * Single source of truth for peak detection: PR #41's canonical breath detector
 * (../absorption-recampaign/breath.mjs detectPeaks), run on the smoothed
 * pre-absorption min(R1,R2). Initial conditions are recomputed from the logged
 * seed via seedChimera (../chimera-campaign/integrator.mjs) with ZERO time
 * integration — the IC is a pure function of the seed. The A=0.2 persistent-run
 * contrast requires re-tracing (those runs were never trace-dumped); that is the
 * only integration this script performs, done with the shipped-identical tracer
 * and gated for bit-for-bit reproduction of the logged campaign labels.
 *
 * Emits (into transient_results/, all deterministic):
 *   cp1_mk_a05.jsonl    per-run M_k sequence for traced A=0.5 runs (from disk)
 *   cp1_mk_a02.jsonl    per-run M_k sequence for re-traced A=0.2 persistent runs
 *   cp2_features.jsonl  per-run t=0 collective IC + first-cycle features + t_abs
 *   determinism_gate.json   re-trace vs logged-label check
 *   coverage.json       CP0 light-audit: trace coverage per (N,A)
 *
 * Usage: node tools/transient-tests/extract.mjs
 */
import {
  readFileSync,
  writeFileSync,
  mkdirSync,
  createWriteStream,
} from 'node:fs';
import { resolve, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import {
  mulberry32,
  seedChimera,
  orderParam,
} from '../chimera-campaign/integrator.mjs';
import {
  detectPeaks,
  movingAverage,
} from '../absorption-recampaign/breath.mjs';
import { traceRun } from '../absorption-recampaign/tracer.mjs';

const __dirname = dirname(fileURLToPath(import.meta.url));
const ROOT = process.cwd();
const cfg = JSON.parse(
  readFileSync(resolve(__dirname, 'transient.config.json'), 'utf8'),
);
const OUT = resolve(ROOT, cfg.output_dir);
mkdirSync(OUT, { recursive: true });

const TAU = 2 * Math.PI;
const wrapPi = (d) => {
  // wrap to (-pi, pi]
  let x = ((d % TAU) + TAU) % TAU;
  if (x > Math.PI) x -= TAU;
  return x;
};

/** Read a JSONL file into an array of parsed objects. */
function readJsonl(path) {
  const txt = readFileSync(resolve(ROOT, path), 'utf8');
  const rows = [];
  for (const line of txt.split('\n')) {
    const s = line.trim();
    if (s) rows.push(JSON.parse(s));
  }
  return rows;
}

/**
 * Collective IC features for a logged seed — pure function of the seed, NO
 * integration. pop1 is the seeded synchronized cluster, pop2 incoherent.
 */
function icFeatures(Np, seed) {
  const ph = seedChimera(Np, mulberry32(seed));
  const op1 = orderParam(ph, 0, Np);
  const op2 = orderParam(ph, Np, Np);
  const Rsync0 = Math.max(op1.R, op2.R);
  const Rincoh0 = Math.min(op1.R, op2.R);
  const dphi0 = wrapPi(op1.Phi - op2.Phi);
  return { Rsync0, Rincoh0, dphi0, absdphi0: Math.abs(dphi0) };
}

/**
 * M_k sequence from a pre-absorption min(R1,R2) prefix using the canonical
 * detector. Returns smoothed peak heights M_k, raw heights, peak times, and the
 * first-cycle features (M1, Tb1). cyclesCompleted = nPeaks - 1.
 */
function mkSequence(rPre, sampleDt) {
  const { peaks } = detectPeaks(rPre, sampleDt, cfg.breath);
  const w = Math.max(1, Math.round(cfg.breath.smoothWindowSec / sampleDt));
  const sm = movingAverage(rPre, w);
  const Mk = peaks.map((p) => sm[p]);
  const MkRaw = peaks.map((p) => rPre[p]);
  const ptimes = peaks.map((p) => p * sampleDt);
  const M1 = Mk.length >= 1 ? Mk[0] : null;
  const Tb1 = ptimes.length >= 2 ? ptimes[1] - ptimes[0] : null;
  return {
    nPeaks: peaks.length,
    cycles: Math.max(0, peaks.length - 1),
    peakTimes: ptimes,
    Mk,
    MkRaw,
    M1,
    Tb1,
  };
}

// ---------------------------------------------------------------------------
// CP0 — light audit: trace coverage per (N, A).
// ---------------------------------------------------------------------------
const traces = readJsonl(cfg.inputs.phase_traces_jsonl);
const campaign = readJsonl(cfg.inputs.campaign_jsonl);

const coverage = {};
for (const t of traces) {
  const k = `N${t.N}_A${t.A}`;
  coverage[k] ??= {
    N: t.N,
    A: t.A,
    records: 0,
    with_full_trace: 0,
    absorbed_among_traced: 0,
  };
  coverage[k].records++;
  if (Array.isArray(t.R_incoh) && t.R_incoh.length) {
    coverage[k].with_full_trace++;
    if (!t.abs_censored) coverage[k].absorbed_among_traced++;
  }
}
const campaignCounts = {};
for (const r of campaign) {
  const k = `N${r.N}_A${r.A}`;
  campaignCounts[k] ??= { N: r.N, A: r.A, total: 0, persistent: 0 };
  campaignCounts[k].total++;
  if (r.abs_censored) campaignCounts[k].persistent++;
}
writeFileSync(
  resolve(OUT, 'coverage.json'),
  JSON.stringify(
    {
      generated: new Date().toISOString(),
      note: 'CP0 light audit. with_full_trace = trace records carrying the R_incoh array. A=0.5 traces cover N in {8,16,32,64} (the phase subset, 100 lowest-id seeds each); A=0.2 only the absorbed members of that subset were trace-dumped, so persistent A=0.2 runs are re-traced for CP1.',
      trace_coverage: coverage,
      campaign_counts: campaignCounts,
      ic_extraction:
        'IC for any logged seed = seedChimera(Np, mulberry32(seed)); R_incoh0=min(R1,R2), R_sync0=max(R1,R2), dphi0=Phi1-Phi2 at t=0. No integration required.',
      peak_detector:
        'tools/absorption-recampaign/breath.mjs detectPeaks (PR #41), on movingAverage(min(R1,R2), smoothWindowSec) over the pre-absorption prefix R_incoh[0:absIndex].',
    },
    null,
    2,
  ),
);
console.log('CP0 coverage written.');

// ---------------------------------------------------------------------------
// CP1 (A=0.5) + CP2 features — from the existing on-disk traces.
// ---------------------------------------------------------------------------
// Index campaign rows by (N,A,seed) for t_abs lookup.
const campByKey = new Map();
for (const r of campaign) campByKey.set(`${r.N}_${r.A}_${r.seed}`, r);

const theta = cfg.boundary.theta;
const mkA05 = createWriteStream(resolve(OUT, 'cp1_mk_a05.jsonl'));
const featRows = []; // accumulate CP2 features (written after IC pass)
const firstCycleBySeed = new Map(); // (N_seed) -> {M1,Tb1} from traces

for (const t of traces) {
  if (t.A !== 0.5) continue;
  if (!(Array.isArray(t.R_incoh) && t.R_incoh.length)) continue;
  if (t.abs_censored) continue; // all A=0.5 absorbed, but guard anyway
  const absIdx = t.absIndex > 0 ? t.absIndex : t.R_incoh.length;
  const pre = t.R_incoh.slice(0, Math.min(absIdx, t.R_incoh.length));
  const seq = mkSequence(pre, t.sampleDt);
  mkA05.write(
    JSON.stringify({
      N: t.N,
      A: t.A,
      seed: t.seed,
      t_abs: t.t_abs,
      abs_censored: t.abs_censored,
      sampleDt: t.sampleDt,
      theta,
      n_peaks: seq.nPeaks,
      cycles: seq.cycles,
      peak_times: seq.peakTimes,
      Mk: seq.Mk,
      Mk_raw: seq.MkRaw,
    }) + '\n',
  );
  firstCycleBySeed.set(`${t.N}_${t.seed}`, { M1: seq.M1, Tb1: seq.Tb1 });
}
await new Promise((r) => mkA05.end(r));
console.log('CP1 A=0.5 M_k sequences written.');

// CP2 features for EVERY A=0.5 campaign run (IC from seed; first-cycle where traced).
for (const r of campaign) {
  if (r.A !== 0.5) continue;
  const ic = icFeatures(r.N, r.seed);
  const fc = firstCycleBySeed.get(`${r.N}_${r.seed}`) || {
    M1: null,
    Tb1: null,
  };
  featRows.push({
    N: r.N,
    A: r.A,
    seed: r.seed,
    t_abs: r.t_abs,
    abs_censored: r.abs_censored,
    Rincoh0: ic.Rincoh0,
    Rsync0: ic.Rsync0,
    dphi0: ic.dphi0,
    absdphi0: ic.absdphi0,
    M1: fc.M1,
    Tb1: fc.Tb1,
    traced: fc.M1 !== null || fc.Tb1 !== null,
  });
}
writeFileSync(
  resolve(OUT, 'cp2_features.jsonl'),
  featRows.map((x) => JSON.stringify(x)).join('\n') + '\n',
);
console.log(`CP2 features written (${featRows.length} A=0.5 runs).`);

// ---------------------------------------------------------------------------
// CP1 (A=0.2 contrast) — re-trace persistent (abs_censored) runs over t_max.
// ---------------------------------------------------------------------------
const label = {
  theta: cfg.graze_criterion.theta,
  W: cfg.graze_criterion.W,
  T_v: cfg.absorption_criterion.T_v,
  recThresh: cfg.absorption_criterion.recoveryThreshold,
  recWin: cfg.absorption_criterion.recoveryWindowSec,
};
const { beta, dt, sampleStride, t_max } = cfg.model;

const a02 = createWriteStream(resolve(OUT, 'cp1_mk_a02.jsonl'));
const detGate = []; // determinism: re-trace vs logged label
let a02Count = 0;
for (const Np of cfg.cp1.a02_Ns) {
  // persistent seeds for this N, lowest-id first
  const persistent = campaign
    .filter((r) => r.A === 0.2 && r.N === Np && r.abs_censored)
    .sort((a, b) => a.seed - b.seed)
    .slice(0, cfg.cp1.a02_persistent_per_point);
  for (const r of persistent) {
    const tr = traceRun({
      Np,
      A: 0.2,
      beta,
      seed: r.seed,
      t_max,
      dt,
      label,
      breath: cfg.breath,
      sampleStride,
      earlyExit: false,
      keepTrace: true,
    });
    // Determinism gate: re-traced labels must match the logged campaign row.
    detGate.push({
      N: Np,
      seed: r.seed,
      logged_t_graze: r.t_graze,
      retrace_t_graze: tr.t_graze,
      logged_graze_censored: r.graze_censored,
      retrace_graze_censored: tr.graze_censored,
      logged_abs_censored: r.abs_censored,
      retrace_abs_censored: tr.abs_censored,
      graze_match:
        r.graze_censored === tr.graze_censored &&
        (r.graze_censored || Math.abs(r.t_graze - tr.t_graze) < 1e-6),
      abs_match: r.abs_censored === tr.abs_censored,
    });
    // Persistent ⇒ pre window is the whole trace (no absorbing cycle to exclude).
    const seq = mkSequence(tr.R_incoh, tr.sampleDt);
    a02.write(
      JSON.stringify({
        N: Np,
        A: 0.2,
        seed: r.seed,
        abs_censored: tr.abs_censored,
        sampleDt: tr.sampleDt,
        theta,
        window_t: t_max,
        n_peaks: seq.nPeaks,
        cycles: seq.cycles,
        peak_times: seq.peakTimes,
        Mk: seq.Mk,
        Mk_raw: seq.MkRaw,
      }) + '\n',
    );
    a02Count++;
  }
}
await new Promise((r) => a02.end(r));
console.log(`CP1 A=0.2 re-traced ${a02Count} persistent runs.`);

const grazeOk = detGate.filter((d) => d.graze_match).length;
const absOk = detGate.filter((d) => d.abs_match).length;
writeFileSync(
  resolve(OUT, 'determinism_gate.json'),
  JSON.stringify(
    {
      note: 'Re-traced A=0.2 persistent runs vs logged absorption_campaign rows. graze_match requires identical censor flag and (if uncensored) t_graze within 1e-6.',
      n: detGate.length,
      graze_match: `${grazeOk}/${detGate.length}`,
      abs_censor_match: `${absOk}/${detGate.length}`,
      pass: grazeOk === detGate.length && absOk === detGate.length,
      sample: detGate.slice(0, 8),
    },
    null,
    2,
  ),
);
console.log(
  `Determinism gate: graze ${grazeOk}/${detGate.length}, abs-censor ${absOk}/${detGate.length}.`,
);
console.log('extract.mjs done.');
