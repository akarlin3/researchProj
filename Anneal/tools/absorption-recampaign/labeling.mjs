/**
 * CP1 — two-timescale labeling of an R_incoh(t) trace.
 *
 * PR #41 showed the published campaign's "collapse" criterion — the first
 * W-sustained θ-crossing of R_incoh = min(R₁,R₂) — is NOT an absorption time: up
 * to 98% of such crossings recover (the chimera reforms). This module derives
 * BOTH labels from one trace so they are directly comparable per run:
 *
 *   t_graze — the published criterion, UNCHANGED: the run-start time of the first
 *             maximal above-θ run that lasts ≥ W. (collapse.mjs verbatim.)
 *
 *   t_abs   — ABSORPTION-grade: the run-start time of the first W-sustained
 *             θ-crossing that is NOT followed by RECOVERY within a verification
 *             horizon T_v. Recovery = R_incoh drops back below `recThresh`
 *             sustained for ≥ `recWin` seconds after the crossing is confirmed
 *             (PR #41 CP3's recovery definition). A θ-crossing that DOES recover
 *             is a graze; we count it (n_grazes_before_abs) and keep searching.
 *             A crossing whose full T_v window is not yet available at t_max
 *             leaves the run RIGHT-CENSORED for t_abs.
 *
 * Dating convention is identical to the campaign runner (runner.mjs): a crossing
 * is dated to the integration step where its above-θ run *began*. Because samples
 * are uniform (index·sampleDt = step·dt), labelling the sampled series by
 * index·sampleDt reproduces the runner's runStartStep·dt bit-for-bit — this is
 * what makes the t_graze determinism gate pass.
 *
 * The `Labeler` is the single source of truth: the online campaign tracer feeds
 * it sample-by-sample and early-exits when `.absorptionConfirmed` is true, while
 * `labelSeries` runs it over a whole array for the unit tests. One code path.
 *
 * Pure dynamics post-processing; no audio, no integrator state.
 */

export const DEFAULTS = {
  theta: 0.85, // graze threshold θ (collapse.mjs DEFAULT_THETA)
  W: 5.0, // graze sustain window (collapse.mjs DEFAULT_W)
  T_v: 120.0, // absorption verification horizon (s)
  recThresh: 0.8, // recovery threshold on R_incoh
  recWin: 5.0, // recovery sustain window (s)
};

/**
 * Construct a streaming two-timescale labeler.
 *
 * @param {number} sampleDt  model-time seconds between pushed samples
 * @param {object} [opts]    overrides for {theta, W, T_v, recThresh, recWin}
 */
export function makeLabeler(sampleDt, opts = {}) {
  const { theta, W, T_v, recThresh, recWin } = { ...DEFAULTS, ...opts };
  const need = Math.max(1, Math.round(W / sampleDt)); // samples above θ ⇒ crossing
  const recNeed = Math.max(1, Math.round(recWin / sampleDt)); // samples below recThresh ⇒ recovery
  const vHorizon = Math.max(1, Math.round(T_v / sampleDt)); // samples post-confirmation to certify absorption

  let idx = -1;

  // SEARCHING-phase state: detect the next W-sustained above-θ run.
  let consec = 0;
  let runStart = -1;

  // First crossing ⇒ t_graze.
  let firstCrossStartIdx = -1;

  // VERIFYING-phase state for the current pending crossing.
  let phase = 'searching'; // 'searching' | 'verifying'
  let crossStartIdx = -1; // run-start index of the pending crossing
  let confIdx = -1; // confirmation index (run-start + need - 1)
  let belowConsec = 0; // consecutive samples below recThresh since confirmation
  let grazeCount = 0; // resolved grazes so far

  // Absorption result.
  let absStartIdx = -1;
  let absGrazes = -1;
  let done = false;

  function push(r) {
    idx++;
    if (done) return;

    if (phase === 'searching') {
      if (r > theta) {
        if (consec === 0) runStart = idx;
        consec++;
        if (consec >= need) {
          crossStartIdx = runStart;
          if (firstCrossStartIdx < 0) firstCrossStartIdx = runStart;
          confIdx = idx;
          belowConsec = 0;
          phase = 'verifying';
        }
      } else {
        consec = 0;
        runStart = -1;
      }
      return;
    }

    // phase === 'verifying': wait for recovery (graze) or T_v elapse (absorb).
    if (r < recThresh) {
      belowConsec++;
      if (belowConsec >= recNeed) {
        // Recovery → the pending crossing was a graze. Resume searching.
        grazeCount++;
        phase = 'searching';
        consec = 0;
        runStart = -1;
        return;
      }
    } else {
      belowConsec = 0;
    }
    // Recovery is checked first (any recovery inside the closed T_v window wins);
    // only if none has occurred and the full horizon has elapsed do we absorb.
    if (idx - confIdx >= vHorizon) {
      absStartIdx = crossStartIdx;
      absGrazes = grazeCount;
      done = true;
    }
  }

  /**
   * Finalize. `tMax` is the censoring time recorded when a label is censored
   * (carried as the lifetime, matching the campaign's convention).
   */
  function result(tMax) {
    const grazeCensored = firstCrossStartIdx < 0;
    const t_graze = grazeCensored ? tMax : firstCrossStartIdx * sampleDt;
    let t_abs, absCensored, nGraze;
    if (absStartIdx >= 0) {
      t_abs = absStartIdx * sampleDt;
      absCensored = false;
      nGraze = absGrazes;
    } else {
      t_abs = tMax;
      absCensored = true;
      nGraze = grazeCount; // grazes observed before censoring (lower bound)
    }
    return {
      t_graze,
      graze_censored: grazeCensored,
      t_abs,
      abs_censored: absCensored,
      n_grazes_before_abs: nGraze,
    };
  }

  return {
    push,
    result,
    get absorptionConfirmed() {
      return done;
    },
    get nSamples() {
      return idx + 1;
    },
  };
}

/**
 * Label a whole sampled R_incoh series. Convenience wrapper over `makeLabeler`
 * used by the unit tests; the campaign tracer drives `makeLabeler` online.
 *
 * @param {Float64Array|number[]} R  R_incoh sampled on a uniform grid
 * @param {number} sampleDt          seconds between samples
 * @param {object} [opts]            label overrides
 * @param {number} [tMax]            censoring time (default last sample time)
 */
export function labelSeries(R, sampleDt, opts = {}, tMax = null) {
  const L = makeLabeler(sampleDt, opts);
  for (let i = 0; i < R.length; i++) L.push(R[i]);
  const horizon = tMax == null ? (R.length - 1) * sampleDt : tMax;
  return L.result(horizon);
}
