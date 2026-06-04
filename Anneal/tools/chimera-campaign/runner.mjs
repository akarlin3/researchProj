/**
 * Online collapse-time runner — integrates a trajectory step-by-step, evaluates
 * the CP2 collapse criterion as it goes, and EARLY-EXITS the moment a sustained
 * collapse is confirmed. Collapsing runs therefore cost ~lifetime of compute,
 * not the full t_max, which is what makes the 200-seed campaign tractable.
 *
 * Numerically identical to `runTrajectory` + `collapse.lifetimeOf` (same RK4,
 * same min(R₁,R₂) > θ sustained-for-W rule) — it just doesn't store the whole
 * series. Exactly reproducible from (Np, A, beta, seed, dt, theta, W, t_max).
 */
import {
  mulberry32,
  seedChimera,
  orderParam,
  makeScratch,
  rk4StepInPlace,
  DEFAULT_DT,
} from './integrator.mjs';
import { DEFAULT_THETA, DEFAULT_W } from './collapse.mjs';

/**
 * @param {object} o
 * @param {number} o.Np
 * @param {number} o.A
 * @param {number} o.beta
 * @param {number} o.seed
 * @param {number} o.t_max                 horizon / censoring time (model s)
 * @param {number} [o.dt]                  default shipped 0.05
 * @param {number} [o.theta]               default 0.85
 * @param {number} [o.W]                   default 5.0 s
 * @param {number} [o.sampleStride]        model-time s between criterion samples (default 0.1)
 * @returns {{ lifetime:number, censored:boolean, dt:number, theta:number,
 *             W:number, stepsIntegrated:number }}
 */
export function runUntilCollapse({
  Np,
  A,
  beta,
  seed,
  t_max,
  dt = DEFAULT_DT,
  theta = DEFAULT_THETA,
  W = DEFAULT_W,
  sampleStride = 0.1,
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
  const need = Math.max(1, Math.round(W / sampleDt)); // samples to sustain above θ

  let runStartStep = -1; // integration step where the current above-θ run began
  let consecutive = 0;

  // Evaluate the criterion on a given sample (step index `step`).
  const evalSample = (step) => {
    const r1 = orderParam(phases, 0, Np).R;
    const r2 = orderParam(phases, Np, Np).R;
    const rIncoh = r1 < r2 ? r1 : r2;
    if (rIncoh > theta) {
      if (consecutive === 0) runStartStep = step;
      consecutive++;
      if (consecutive >= need) {
        return runStartStep * dt; // collapse, dated to the run's start
      }
    } else {
      consecutive = 0;
      runStartStep = -1;
    }
    return null;
  };

  // Sample the initial condition.
  let life = evalSample(0);
  if (life !== null) {
    return {
      lifetime: life,
      censored: false,
      dt,
      theta,
      W,
      stepsIntegrated: 0,
    };
  }

  for (let step = 1; step <= nSteps; step++) {
    rk4StepInPlace(phases, Np, mu, nu, alpha, dt, scratch);
    if (step % sampleEvery === 0) {
      life = evalSample(step);
      if (life !== null) {
        return {
          lifetime: life,
          censored: false,
          dt,
          theta,
          W,
          stepsIntegrated: step,
        };
      }
    }
  }
  return {
    lifetime: t_max,
    censored: true,
    dt,
    theta,
    W,
    stepsIntegrated: nSteps,
  };
}
