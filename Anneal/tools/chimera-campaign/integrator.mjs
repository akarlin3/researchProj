/**
 * Headless deterministic two-population Sakaguchi–Kuramoto integrator.
 *
 * Self-contained ESM (no node_modules, no audio graph) so the collapse-time
 * campaign runs at maximum speed in plain `node` — on this machine or on the
 * Origin desktop over SSH. This is a *faithful, verbatim* port of the shipped
 * dynamics core in `src/audio/chimera.ts` (and its offline reference
 * `examples/probes/chimera_probe.mjs`): same mean-field RK4 update, same
 * mulberry32 PRNG, same canonical chimera seed, same order parameter.
 *
 * A vitest cross-check (`tools/chimera-campaign/integrator.crosscheck.test.ts`)
 * asserts this module reproduces `src/audio/chimera.ts` bit-for-bit from the
 * same seed, so there is zero drift from the shipped voice's math.
 *
 * The campaign measures DYNAMICS ONLY — per-population order parameters R₁, R₂.
 * No fusion law, no centroid, no audio rendering anywhere in this path.
 *
 * Model (Abrams, Mirollo, Strogatz & Wiley 2008). Two equal populations
 * σ ∈ {1,2} of `Np` oscillators, identical natural frequency ω = 0 (rotating
 * frame — the identical-ω requirement, δω = 0 exactly):
 *
 *   dθᵢ^σ/dt = μ·R_σ·sin(Φ_σ − θᵢ^σ − α) + ν·R_σ'·sin(Φ_σ' − θᵢ^σ − α)
 *
 * with μ=(1+A)/2, ν=(1−A)/2, α=π/2−β. Integrated with RK4 at dt (default 0.05,
 * the shipped control-rate step).
 */

export const TAU = 2 * Math.PI;

/** Shipped default phase lag β (α = π/2 − β); see chimera.ts DEFAULT_BETA. */
export const DEFAULT_BETA = 0.02;
/** Shipped control-rate timestep (orchestrator DRIFT_DT). */
export const DEFAULT_DT = 0.05;
/** Half-width (rad) of the seeded synchronized cluster at t=0 (chimera.ts). */
export const SEED_JITTER = 0.25;
/** A locked population has order parameter above this (chimera.ts SYNC_HI). */
export const SYNC_HI = 0.9;
/** An incoherent population has order parameter below this (chimera.ts INCOH_LO). */
export const INCOH_LO = 0.85;

/**
 * mulberry32 — verbatim from src/audio/analysis/__tests__/redistribution.test.ts
 * and the probe. Integer seed → deterministic stream in [0,1).
 */
export function mulberry32(seed) {
  let a = seed >>> 0;
  return () => {
    a |= 0;
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/** Wrap a phase to [0, 2π). Matches chimera.ts `wrap`. */
function wrap(v) {
  return ((v % TAU) + TAU) % TAU;
}

/**
 * Complex order parameter R·e^{iΦ} of a contiguous block of `count` phases
 * starting at `start`. Returns { R, Phi }. Verbatim from chimera.ts.
 */
export function orderParam(phases, start, count) {
  let c = 0;
  let s = 0;
  for (let k = 0; k < count; k++) {
    const th = phases[start + k];
    c += Math.cos(th);
    s += Math.sin(th);
  }
  c /= count;
  s /= count;
  return { R: Math.sqrt(c * c + s * s), Phi: Math.atan2(s, c) };
}

/**
 * Canonical chimera seed (chimera.ts `seedChimera`): pop 1 a tight synchronized
 * cluster about a random anchor (R₁≈1, ±jitter rad), pop 2 incoherent (uniform).
 * RNG consumed in the exact order anchor → Np jitter → Np incoherent, so the
 * runner and the shipped voice reproduce identical trajectories from one seed.
 * Returns a Float64Array of length 2·Np.
 */
export function seedChimera(Np, rng, jitter = SEED_JITTER) {
  const n = 2 * Np;
  const phases = new Float64Array(n);
  const anchor = rng() * TAU;
  for (let i = 0; i < Np; i++) {
    phases[i] = wrap(anchor + jitter * (rng() - 0.5));
  }
  for (let i = Np; i < n; i++) {
    phases[i] = rng() * TAU;
  }
  return phases;
}

/** Mean-field derivative for both populations (ω = 0) into `out`. chimera.ts. */
function deriv(phases, Np, mu, nu, alpha, out) {
  const op1 = orderParam(phases, 0, Np);
  const op2 = orderParam(phases, Np, Np);
  for (let i = 0; i < Np; i++) {
    const th = phases[i];
    out[i] =
      mu * op1.R * Math.sin(op1.Phi - th - alpha) +
      nu * op2.R * Math.sin(op2.Phi - th - alpha);
  }
  for (let i = Np; i < 2 * Np; i++) {
    const th = phases[i];
    out[i] =
      mu * op2.R * Math.sin(op2.Phi - th - alpha) +
      nu * op1.R * Math.sin(op1.Phi - th - alpha);
  }
}

/**
 * Advance the two-population system one RK4 step IN PLACE. Mutates `phases`
 * (the runner owns a single buffer per trajectory for speed). Numerically
 * identical to chimera.ts `chimeraStep`, which returns a fresh array instead.
 */
export function rk4StepInPlace(phases, Np, mu, nu, alpha, dt, scratch) {
  const n = phases.length;
  const { k1, k2, k3, k4, tmp } = scratch;
  deriv(phases, Np, mu, nu, alpha, k1);
  for (let i = 0; i < n; i++) tmp[i] = phases[i] + 0.5 * dt * k1[i];
  deriv(tmp, Np, mu, nu, alpha, k2);
  for (let i = 0; i < n; i++) tmp[i] = phases[i] + 0.5 * dt * k2[i];
  deriv(tmp, Np, mu, nu, alpha, k3);
  for (let i = 0; i < n; i++) tmp[i] = phases[i] + dt * k3[i];
  deriv(tmp, Np, mu, nu, alpha, k4);
  for (let i = 0; i < n; i++) {
    const v = phases[i] + (dt / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    phases[i] = wrap(v);
  }
}

export function makeScratch(n) {
  return {
    k1: new Float64Array(n),
    k2: new Float64Array(n),
    k3: new Float64Array(n),
    k4: new Float64Array(n),
    tmp: new Float64Array(n),
  };
}

/**
 * Integrate one trajectory and return the per-population order-parameter time
 * series. Exactly reproducible from its arguments.
 *
 * @param {object} o
 * @param {number} o.Np         oscillators per population
 * @param {number} o.A          coupling disparity
 * @param {number} o.beta       phase lag β (α = π/2 − β)
 * @param {number} o.seed       integer seed for mulberry32
 * @param {number} o.t_max      model-time horizon (units; 1 unit ≔ 1 s)
 * @param {number} [o.dt]       integration step (default shipped 0.05)
 * @param {number} [o.sampleStride] sample R₁,R₂ every this many MODEL TIME units
 *                                   (default 0.1; not every dt step)
 * @returns {{ t:Float64Array, R1:Float64Array, R2:Float64Array, dt:number,
 *             nSteps:number, sampleEvery:number }}
 */
export function runTrajectory({
  Np,
  A,
  beta,
  seed,
  t_max,
  dt = DEFAULT_DT,
  sampleStride = 0.1,
}) {
  const mu = (1 + A) / 2;
  const nu = (1 - A) / 2;
  const alpha = Math.PI / 2 - beta;
  const n = 2 * Np;

  const rng = mulberry32(seed);
  const phases = seedChimera(Np, rng);
  const scratch = makeScratch(n);

  const nSteps = Math.round(t_max / dt);
  // Sample every `sampleEvery` integration steps (≥1) so the cadence is the
  // requested model-time stride regardless of dt.
  const sampleEvery = Math.max(1, Math.round(sampleStride / dt));
  const nSamp = Math.floor(nSteps / sampleEvery) + 1;

  const t = new Float64Array(nSamp);
  const R1 = new Float64Array(nSamp);
  const R2 = new Float64Array(nSamp);

  let si = 0;
  // Sample the initial condition (step 0) too.
  const record = (step) => {
    t[si] = step * dt;
    R1[si] = orderParam(phases, 0, Np).R;
    R2[si] = orderParam(phases, Np, Np).R;
    si++;
  };
  record(0);
  for (let step = 1; step <= nSteps; step++) {
    rk4StepInPlace(phases, Np, mu, nu, alpha, dt, scratch);
    if (step % sampleEvery === 0) record(step);
  }

  return {
    t: t.subarray(0, si),
    R1: R1.subarray(0, si),
    R2: R2.subarray(0, si),
    dt,
    nSteps,
    sampleEvery,
  };
}

/** Runner version — bump when the integrator or sampling contract changes. */
export const RUNNER_VERSION = 'chimera-campaign/1.0.0';
