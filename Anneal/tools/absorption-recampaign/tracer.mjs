/**
 * Absorption-grade tracer — integrates one trajectory with the SHIPPED-identical
 * RK4 core (chimera-campaign/integrator.mjs primitives, verbatim) and drives the
 * single-source-of-truth two-timescale `Labeler` online, so one integration
 * yields BOTH the legacy graze label (t_graze, == published campaign) and the new
 * absorption label (t_abs). Also computes the per-run breath period T_b from the
 * pre-absorption R_incoh.
 *
 * Adds nothing to the dynamics: same mulberry32 seed, same seedChimera, same
 * rk4StepInPlace, same min(R₁,R₂) sampling cadence as the campaign — so t_graze
 * reproduces the logged campaign lifetime bit-for-bit (the determinism gate).
 *
 * Two modes:
 *   traceRun({ earlyExit:true })  — campaign mode: early-exit the moment the
 *       labeler confirms absorption (≈ t_abs + T_v of compute), keep only the
 *       prefix needed for T_b. Used by the full sweep.
 *   traceRun({ earlyExit:false }) — pilot mode: integrate to t_max regardless so
 *       the whole R_incoh series is available to re-label under every (T_v,
 *       recThresh) sensitivity setting.
 *   traceRun({ keepRsync:true })  — also record R_sync = max(R₁,R₂) (for the CP3
 *       phase subset + CP4 supervisor replay).
 */
import {
  mulberry32,
  seedChimera,
  orderParam,
  makeScratch,
  rk4StepInPlace,
  DEFAULT_DT,
} from '../chimera-campaign/integrator.mjs';
import { makeLabeler, DEFAULTS as LABEL_DEFAULTS } from './labeling.mjs';
import { breathPeriod, BREATH_DEFAULTS } from './breath.mjs';

/**
 * @param {object} o
 * @param {number} o.Np
 * @param {number} o.A
 * @param {number} o.beta
 * @param {number} o.seed
 * @param {number} o.t_max
 * @param {number} [o.dt]
 * @param {object} [o.label]   overrides for the labeler {theta,W,T_v,recThresh,recWin}
 * @param {object} [o.breath]  overrides for the breath estimator
 * @param {number} [o.sampleStride]
 * @param {boolean} [o.earlyExit]  default true (campaign); false integrates full t_max
 * @param {boolean} [o.keepRsync]  default false; record R_sync=max(R1,R2) too
 * @param {boolean} [o.keepTrace]  default false; return the full R_incoh (+R_sync) arrays
 */
export function traceRun({
  Np,
  A,
  beta,
  seed,
  t_max,
  dt = DEFAULT_DT,
  label = {},
  breath = {},
  sampleStride = 0.1,
  earlyExit = true,
  keepRsync = false,
  keepTrace = false,
}) {
  const mu = (1 + A) / 2;
  const nu = (1 - A) / 2;
  const alpha = Math.PI / 2 - beta;
  const n = 2 * Np;

  const phases = seedChimera(Np, mulberry32(seed));
  const scratch = makeScratch(n);

  const nSteps = Math.round(t_max / dt);
  const sampleEvery = Math.max(1, Math.round(sampleStride / dt));
  const sampleDt = sampleEvery * dt;

  const labeler = makeLabeler(sampleDt, label);
  const Rincoh = [];
  const Rsync = keepRsync ? [] : null;
  let stepsIntegrated = 0;

  const sampleAt = () => {
    const r1 = orderParam(phases, 0, Np).R;
    const r2 = orderParam(phases, Np, Np).R;
    const lo = r1 < r2 ? r1 : r2;
    Rincoh.push(lo);
    if (Rsync) Rsync.push(r1 < r2 ? r2 : r1);
    labeler.push(lo);
  };

  sampleAt(); // initial condition (step 0)
  for (let step = 1; step <= nSteps; step++) {
    rk4StepInPlace(phases, Np, mu, nu, alpha, dt, scratch);
    stepsIntegrated = step;
    if (step % sampleEvery === 0) {
      sampleAt();
      if (earlyExit && labeler.absorptionConfirmed) break;
    }
  }

  const lab = labeler.result(t_max);

  // Breath period over the pre-absorption window (excludes the absorbing cycle).
  const absIdx = lab.abs_censored
    ? Rincoh.length
    : Math.round(lab.t_abs / sampleDt);
  const pre = Rincoh.slice(0, Math.min(absIdx, Rincoh.length));
  const bp = breathPeriod(pre, sampleDt, breath);

  const out = {
    N: Np,
    A,
    seed,
    t_graze: lab.t_graze,
    graze_censored: lab.graze_censored,
    t_abs: lab.t_abs,
    abs_censored: lab.abs_censored,
    n_grazes_before_abs: lab.n_grazes_before_abs,
    T_b: bp.Tb,
    n_breath_peaks: bp.nPeaks,
    breath_cycles: bp.cyclesCompleted,
    sampleDt,
    dt,
    stepsIntegrated,
    nSamples: Rincoh.length,
  };
  if (keepTrace) {
    out.R_incoh = Rincoh;
    if (Rsync) out.R_sync = Rsync;
    out.absIndex = lab.abs_censored ? -1 : absIdx;
    out.grazeIndex = lab.graze_censored
      ? -1
      : Math.round(lab.t_graze / sampleDt);
  }
  return out;
}

export { LABEL_DEFAULTS, BREATH_DEFAULTS };
