/**
 * CP2 — paper-grade, dynamics-level collapse criterion (decoupled from the
 * shipped supervisor).
 *
 * Collapse of a seeded chimera is the loss of the incoherent population: the two
 * populations merge to global sync, so the *incoherent* population's order
 * parameter R rises. We track the incoherent population as the WEAKER of the two
 * at each instant, R_incoh(t) = min(R₁, R₂) — robust to role swaps. Collapse is
 * declared at the first time R_incoh sustainedly exceeds a threshold θ for a
 * window W:
 *
 *   collapse_time = first t₀ such that R_incoh(t) > θ for all t ∈ [t₀, t₀+W].
 *
 * Lifetime = that first-crossing time t₀. If no such sustained excursion occurs
 * by t_max, the lifetime is RIGHT-CENSORED at t_max.
 *
 * Defaults (motivated by CP0):
 *   θ = 0.85  — the shipped supervisor's INCOH_LO: "the incoherent population is
 *               no longer incoherent." (SYNC_HI=0.9, INCOH_LO=0.85.)
 *   W = 5.0 s — ≫ the breath's fast phase (breathing period ~20–70 s, the fast
 *               relaxation return is a few s), so momentary breathing dips that
 *               push R_incoh up do NOT false-trigger. (The supervisor itself
 *               uses a shorter 2 s hold because it can afford to re-perturb;
 *               the paper criterion is deliberately more conservative.)
 */

export const DEFAULT_THETA = 0.85;
export const DEFAULT_W = 5.0;

/**
 * Compute the collapse lifetime from a sampled incoherent-R series.
 *
 * @param {Float64Array|number[]} rIncoh  R_incoh sampled on a uniform grid
 * @param {number} sampleDt               model-time seconds between samples
 * @param {number} theta                  collapse threshold θ
 * @param {number} W                      sustain window (s)
 * @param {number} tMax                   horizon; censoring time
 * @returns {{ lifetime:number, censored:boolean }}
 */
export function collapseTime(rIncoh, sampleDt, theta, W, tMax) {
  const need = Math.max(1, Math.round(W / sampleDt)); // samples that must stay above θ
  let runStart = -1; // index where the current above-θ run began
  for (let i = 0; i < rIncoh.length; i++) {
    if (rIncoh[i] > theta) {
      if (runStart < 0) runStart = i;
      // Sustained for the full window ⇒ collapse, dated to the run's start.
      if (i - runStart + 1 >= need) {
        return { lifetime: runStart * sampleDt, censored: false };
      }
    } else {
      runStart = -1; // dropped back below θ; reset
    }
  }
  return { lifetime: tMax, censored: true };
}

/** Pointwise incoherent-population R: min(R₁, R₂). */
export function incoherentR(R1, R2) {
  const out = new Float64Array(R1.length);
  for (let i = 0; i < R1.length; i++) out[i] = Math.min(R1[i], R2[i]);
  return out;
}

/**
 * Lifetime of one trajectory result (from runTrajectory) under (θ, W).
 * @param {{R1:Float64Array,R2:Float64Array,t:Float64Array}} run
 */
export function lifetimeOf(run, theta, W, tMax) {
  const sampleDt = run.t.length > 1 ? run.t[1] - run.t[0] : 0;
  const rIncoh = incoherentR(run.R1, run.R2);
  return collapseTime(rIncoh, sampleDt, theta, W, tMax);
}
