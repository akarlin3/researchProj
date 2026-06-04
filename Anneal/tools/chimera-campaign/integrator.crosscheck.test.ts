import { describe, expect, it } from 'vitest';
import {
  chimeraStep,
  isChimeraAlive,
  orderParam as shippedOrderParam,
  seedChimera as shippedSeed,
} from '@/audio/chimera';
import {
  mulberry32,
  orderParam as campaignOrderParam,
  runTrajectory,
  seedChimera as campaignSeed,
  makeScratch,
  rk4StepInPlace,
} from './integrator.mjs';

/**
 * Zero-drift guarantee: the headless campaign integrator
 * (`tools/chimera-campaign/integrator.mjs`) must reproduce the SHIPPED dynamics
 * core (`src/audio/chimera.ts`) bit-for-bit from the same integer seed. If this
 * ever fails, the campaign is no longer measuring the shipped voice.
 *
 * Both consume mulberry32 in the same order, so a shared seed → identical seed
 * phases → identical RK4 trajectory → identical order parameters.
 */
describe('campaign integrator matches shipped src/audio/chimera.ts', () => {
  const Np = 16;
  const A = 0.5;
  const beta = 0.05;
  const dt = 0.05;

  it('seedChimera produces identical initial phases', () => {
    for (const seed of [1, 42, 5000, 9001]) {
      const a = shippedSeed(Np, mulberry32(seed));
      const b = campaignSeed(Np, mulberry32(seed));
      expect(b.length).toBe(a.length);
      for (let i = 0; i < a.length; i++) {
        expect(b[i]).toBe(a[i]);
      }
    }
  });

  it('RK4 trajectory matches the shipped chimeraStep exactly over 4000 steps', () => {
    const seed = 5000;
    // Shipped: pure, returns a fresh array each step.
    let shipped = shippedSeed(Np, mulberry32(seed));
    // Campaign: in-place on a single buffer.
    const campaign = campaignSeed(Np, mulberry32(seed));
    const scratch = makeScratch(2 * Np);
    const mu = (1 + A) / 2;
    const nu = (1 - A) / 2;
    const alpha = Math.PI / 2 - beta;

    for (let step = 0; step < 4000; step++) {
      const s = chimeraStep(shipped, { Np, A, beta }, dt);
      shipped = s.phases;
      rk4StepInPlace(campaign, Np, mu, nu, alpha, dt, scratch);
      if (step % 250 === 0 || step === 3999) {
        for (let i = 0; i < campaign.length; i++) {
          expect(campaign[i]).toBe(shipped[i]);
        }
      }
    }
  });

  it('order parameter agrees and runTrajectory tracks the shipped R₁,R₂', () => {
    const seed = 9001;
    const t_max = 60;
    const run = runTrajectory({
      Np,
      A,
      beta,
      seed,
      t_max,
      dt,
      sampleStride: 0.1,
    });

    // Re-derive the same series via the shipped step, sampling at the same cadence.
    let shipped = shippedSeed(Np, mulberry32(seed));
    const sampleEvery = run.sampleEvery;
    const nSteps = Math.round(t_max / dt);
    // index 0 = initial condition
    let si = 0;
    const cmp = (phases) => {
      const op1 = shippedOrderParam(phases, 0, Np);
      const op2 = shippedOrderParam(phases, Np, Np);
      // campaign order parameter is the same function — sanity check it too
      const cop1 = campaignOrderParam(phases, 0, Np);
      expect(cop1.R).toBe(op1.R);
      expect(run.R1[si]).toBeCloseTo(op1.R, 12);
      expect(run.R2[si]).toBeCloseTo(op2.R, 12);
      si++;
    };
    cmp(shipped);
    for (let step = 1; step <= nSteps; step++) {
      shipped = chimeraStep(shipped, { Np, A, beta }, dt).phases;
      if (step % sampleEvery === 0) cmp(shipped);
    }
    expect(si).toBe(run.R1.length);
  });

  it('isChimeraAlive thresholds match (SYNC_HI=0.9, INCOH_LO=0.85)', () => {
    expect(isChimeraAlive({ R: 0.95, Phi: 0 }, { R: 0.5, Phi: 1 })).toBe(true);
    expect(isChimeraAlive({ R: 0.95, Phi: 0 }, { R: 0.88, Phi: 1 })).toBe(
      false,
    );
    expect(isChimeraAlive({ R: 0.8, Phi: 0 }, { R: 0.5, Phi: 1 })).toBe(false);
  });
});
