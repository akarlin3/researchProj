/**
 * Per-run breath-period estimator — JS port of PR #41's Python estimator
 * (tools/breath-phase/analysis.py: moving_average / auto_period / detect_peaks),
 * so T_b can be written into every campaign row instead of requiring the full
 * R_incoh trace downstream. A numerical cross-check against the Python estimator
 * is reported by the pilot (breath_crosscheck).
 *
 * T_b = median peak-to-peak interval of R_incoh maxima over the pre-collapse
 * window, EXCLUDING the final (collapsing) cycle — the diff of detected peak
 * times naturally drops the trailing partial cycle. Peaks are detected on
 * R_incoh smoothed with a centered ~2 s moving average; a candidate survives if
 * its topographic prominence exceeds `promFrac` of the pre-collapse range and it
 * is at least `sepFrac · P_auto` from a taller neighbour, where P_auto is the
 * dominant autocorrelation period (self-tuning; no breath period hard-coded).
 *
 * find_peaks parity: SciPy applies the distance filter (keep the taller of any
 * two peaks closer than `distance`) BEFORE the prominence filter, and computes
 * prominence against the full signal. This port does the same.
 */

export const BREATH_DEFAULTS = {
  smoothWindowSec: 2.0,
  minProminenceFrac: 0.1,
  minPeakSepFrac: 0.5,
  minCyclesForPhase: 2,
  autoPeriodFloorSec: 8.0,
};

/** Centered moving average, window `w` samples (forced odd), edge-padded. */
export function movingAverage(x, w) {
  w = Math.max(1, w | 0);
  if (w % 2 === 0) w += 1; // force odd (matches numpy `int(w)|1`)
  if (w === 1) return Array.from(x);
  const pad = (w - 1) / 2;
  const n = x.length;
  const out = new Array(n);
  for (let i = 0; i < n; i++) {
    let s = 0;
    for (let k = -pad; k <= pad; k++) {
      let j = i + k;
      if (j < 0)
        j = 0; // edge padding
      else if (j >= n) j = n - 1;
      s += x[j];
    }
    out[i] = s / w;
  }
  return out;
}

/**
 * Dominant period (s) from the unbiased autocorrelation of mean-subtracted x —
 * the first local maximum with correlation > 0.1 in [2 s, 120 s]. Only used to
 * self-tune the peak min-separation, so a rough value is fine.
 */
export function autoPeriod(x, dt) {
  const n = x.length;
  if (n < Math.round(6.0 / dt)) return NaN;
  let mean = 0;
  for (let i = 0; i < n; i++) mean += x[i];
  mean /= n;
  const xc = new Array(n);
  for (let i = 0; i < n; i++) xc[i] = x[i] - mean;
  // Unbiased autocorrelation at lag k (k = 0..n-1).
  const ac0 = (() => {
    let s = 0;
    for (let i = 0; i < n; i++) s += xc[i] * xc[i];
    return s / n;
  })();
  if (ac0 <= 0) return NaN;
  const lagMin = Math.max(1, Math.round(2.0 / dt));
  const lagMax = Math.min(n - 2, Math.round(120.0 / dt));
  if (lagMax <= lagMin) return NaN;
  const ac = (k) => {
    let s = 0;
    for (let i = 0; i < n - k; i++) s += xc[i] * xc[i + k];
    return s / (n - k) / ac0; // unbiased, normalized by ac[0]
  };
  let prev = ac(lagMin - 1);
  let cur = ac(lagMin);
  for (let k = lagMin; k < lagMax; k++) {
    const next = ac(k + 1);
    if (cur > prev && cur >= next && cur > 0.1) return k * dt;
    prev = cur;
    cur = next;
  }
  return NaN;
}

/** All strict/plateau local maxima indices (plateau ⇒ left edge, like SciPy). */
function localMaxima(x) {
  const peaks = [];
  const n = x.length;
  let i = 1;
  while (i < n - 1) {
    if (x[i - 1] < x[i]) {
      let ahead = i + 1;
      while (ahead < n - 1 && x[ahead] === x[i]) ahead++; // walk plateau
      if (x[ahead] < x[i]) {
        peaks.push((i + ahead - 1) >> 1); // plateau midpoint
        i = ahead;
        continue;
      }
    }
    i++;
  }
  return peaks;
}

/** Topographic prominence of each peak against the full signal (SciPy wlen=None). */
function prominences(x, peaks) {
  const n = x.length;
  return peaks.map((p) => {
    const h = x[p];
    // left base: walk left to the nearest sample >= h (or the array start),
    // tracking the minimum in between.
    let leftMin = h;
    for (let j = p - 1; j >= 0; j--) {
      if (x[j] > h) break;
      if (x[j] < leftMin) leftMin = x[j];
    }
    let rightMin = h;
    for (let j = p + 1; j < n; j++) {
      if (x[j] > h) break;
      if (x[j] < rightMin) rightMin = x[j];
    }
    return h - Math.max(leftMin, rightMin);
  });
}

/** SciPy distance filter: keep taller peaks, drop any within `distance`. */
function selectByDistance(peaks, heights, distance) {
  const keep = new Array(peaks.length).fill(true);
  // priority = ascending height; iterate tallest first.
  const order = [...peaks.keys()].sort((a, b) => heights[a] - heights[b]);
  for (let oi = order.length - 1; oi >= 0; oi--) {
    const j = order[oi];
    if (!keep[j]) continue;
    // suppress neighbours within `distance` on both sides
    let k = j - 1;
    while (k >= 0 && peaks[j] - peaks[k] < distance) {
      keep[k] = false;
      k--;
    }
    k = j + 1;
    while (k < peaks.length && peaks[k] - peaks[j] < distance) {
      keep[k] = false;
      k++;
    }
  }
  return peaks.filter((_, i) => keep[i]);
}

/**
 * Detect breath maxima on the pre-collapse R_incoh. Returns { peaks, pAuto }.
 * Mirrors analysis.py detect_peaks.
 */
export function detectPeaks(rPre, dt, opts = {}) {
  const o = { ...BREATH_DEFAULTS, ...opts };
  if (rPre.length < Math.round(6.0 / dt)) return { peaks: [], pAuto: NaN };
  const sm = movingAverage(rPre, Math.round(o.smoothWindowSec / dt));
  let mn = Infinity,
    mx = -Infinity;
  for (const v of sm) {
    if (v < mn) mn = v;
    if (v > mx) mx = v;
  }
  const rng = mx - mn;
  if (rng <= 1e-6) return { peaks: [], pAuto: NaN };
  const pAuto = autoPeriod(sm, dt);
  const distance = Number.isFinite(pAuto)
    ? Math.max(1, Math.round((o.minPeakSepFrac * pAuto) / dt))
    : Math.max(1, Math.round(o.autoPeriodFloorSec / dt));

  let peaks = localMaxima(sm);
  if (peaks.length === 0) return { peaks: [], pAuto };
  // distance filter BEFORE prominence (SciPy order)
  const heights = peaks.map((p) => sm[p]);
  peaks = selectByDistance(peaks, heights, distance);
  // prominence filter
  const proms = prominences(sm, peaks);
  const promFloor = o.minProminenceFrac * rng;
  peaks = peaks.filter((_, i) => proms[i] >= promFloor);
  return { peaks, pAuto };
}

/**
 * Per-run breath period. `rPre` is R_incoh over the pre-collapse window
 * (excludes the absorbing/collapsing cycle by passing the pre-event prefix).
 *
 * @returns {{ Tb:number|null, nPeaks:number, cyclesCompleted:number,
 *             intervals:number[], pAuto:number|null }}
 */
export function breathPeriod(rPre, dt, opts = {}) {
  const o = { ...BREATH_DEFAULTS, ...opts };
  const { peaks, pAuto } = detectPeaks(rPre, dt, o);
  const nPeaks = peaks.length;
  const cyclesCompleted = Math.max(0, nPeaks - 1);
  if (cyclesCompleted < o.minCyclesForPhase) {
    return {
      Tb: null,
      nPeaks,
      cyclesCompleted,
      intervals: [],
      pAuto: Number.isFinite(pAuto) ? pAuto : null,
      peaks,
    };
  }
  const ptimes = peaks.map((p) => p * dt);
  const intervals = [];
  for (let i = 1; i < ptimes.length; i++)
    intervals.push(ptimes[i] - ptimes[i - 1]);
  const sorted = [...intervals].sort((a, b) => a - b);
  const m = sorted.length;
  const Tb =
    m % 2 ? sorted[(m - 1) / 2] : (sorted[m / 2 - 1] + sorted[m / 2]) / 2;
  return {
    Tb,
    nPeaks,
    cyclesCompleted,
    intervals,
    pAuto: Number.isFinite(pAuto) ? pAuto : null,
    peaks,
  };
}
