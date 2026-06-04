import { describe, expect, it } from 'vitest';
import { collapseTime, incoherentR } from './collapse.mjs';

/**
 * Collapse-criterion unit tests on synthetic R(t) traces — the criterion logic
 * (sustained-excursion detection, censoring, run reset) is verified in isolation
 * from the integrator.
 */
describe('collapseTime — sustained-excursion detection on synthetic R(t)', () => {
  const sampleDt = 0.1;
  const theta = 0.85;
  const W = 5.0;
  const tMax = 100;

  // Build a trace sampled at sampleDt that is `low` until t=onset, then `high`.
  function stepTrace(onset, low = 0.4, high = 0.97, total = tMax) {
    const n = Math.round(total / sampleDt) + 1;
    const r = new Float64Array(n);
    for (let i = 0; i < n; i++) r[i] = i * sampleDt >= onset ? high : low;
    return r;
  }

  it('censors a trace that never crosses θ', () => {
    const r = new Float64Array(1001).fill(0.4);
    const { lifetime, censored } = collapseTime(r, sampleDt, theta, W, tMax);
    expect(censored).toBe(true);
    expect(lifetime).toBe(tMax);
  });

  it('dates the lifetime to the START of the sustained excursion', () => {
    const r = stepTrace(30);
    const { lifetime, censored } = collapseTime(r, sampleDt, theta, W, tMax);
    expect(censored).toBe(false);
    expect(lifetime).toBeCloseTo(30, 6);
  });

  it('does NOT trigger on a brief spike shorter than W', () => {
    // spike above θ for 2 s (< W=5 s) then back down forever
    const n = Math.round(tMax / sampleDt) + 1;
    const r = new Float64Array(n).fill(0.4);
    for (let i = 0; i < n; i++) {
      const t = i * sampleDt;
      if (t >= 20 && t < 22) r[i] = 0.97;
    }
    const { censored } = collapseTime(r, sampleDt, theta, W, tMax);
    expect(censored).toBe(true);
  });

  it('a spike resets the run, so a later sustained excursion is timed correctly', () => {
    const n = Math.round(tMax / sampleDt) + 1;
    const r = new Float64Array(n).fill(0.4);
    for (let i = 0; i < n; i++) {
      const t = i * sampleDt;
      if (t >= 10 && t < 12) r[i] = 0.97; // brief spike, must be ignored
      if (t >= 50) r[i] = 0.97; // real collapse
    }
    const { lifetime, censored } = collapseTime(r, sampleDt, theta, W, tMax);
    expect(censored).toBe(false);
    expect(lifetime).toBeCloseTo(50, 6);
  });

  it('θ and W move the lifetime monotonically as expected', () => {
    // A ramp that slowly rises: higher θ ⇒ later crossing.
    const n = Math.round(tMax / sampleDt) + 1;
    const r = new Float64Array(n);
    for (let i = 0; i < n; i++) r[i] = Math.min(1, 0.3 + (i * sampleDt) / 100);
    const lo = collapseTime(r, sampleDt, 0.8, 1, tMax).lifetime;
    const hi = collapseTime(r, sampleDt, 0.9, 1, tMax).lifetime;
    expect(hi).toBeGreaterThan(lo);
  });

  it('incoherentR takes the pointwise minimum (role-swap robust)', () => {
    const R1 = new Float64Array([0.9, 0.2, 0.95]);
    const R2 = new Float64Array([0.3, 0.97, 0.1]);
    const inc = incoherentR(R1, R2);
    expect(Array.from(inc)).toEqual([0.3, 0.2, 0.1]);
  });
});
