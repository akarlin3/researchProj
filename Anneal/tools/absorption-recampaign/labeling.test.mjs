/**
 * CP1 unit tests — two-timescale labeling on synthetic R_incoh(t).
 *
 * Run: node --test tools/absorption-recampaign/labeling.test.mjs
 *
 * Covers the five required shapes plus a shallow-dip edge case (a sub-θ dip that
 * never reaches the recovery band must NOT count as recovery).
 */
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { labelSeries, DEFAULTS } from './labeling.mjs';

const DT = 0.1;
const HI = 0.95; // above θ=0.85
const SHALLOW = 0.82; // below θ but above recThresh=0.80 (not a recovery)
const LO = 0.3; // below recThresh ⇒ recovery when sustained

/** Concatenate constant-value segments specified in seconds. */
function build(segments) {
  const out = [];
  for (const [val, secs] of segments) {
    const n = Math.round(secs / DT);
    for (let i = 0; i < n; i++) out.push(val);
  }
  return out;
}

const near = (a, b, tol = DT + 1e-9) => Math.abs(a - b) <= tol;

test('never-collapse: both labels censored', () => {
  const R = build([[LO, 300]]);
  const r = labelSeries(R, DT, {}, 300);
  assert.equal(r.graze_censored, true);
  assert.equal(r.abs_censored, true);
  assert.equal(r.n_grazes_before_abs, 0);
});

test('immediate-absorb: t_graze≈0, t_abs≈0, no grazes', () => {
  const R = build([[HI, 300]]);
  const r = labelSeries(R, DT, {}, 300);
  assert.equal(r.graze_censored, false);
  assert.equal(r.abs_censored, false);
  assert.ok(near(r.t_graze, 0), `t_graze=${r.t_graze}`);
  assert.ok(near(r.t_abs, 0), `t_abs=${r.t_abs}`);
  assert.equal(r.n_grazes_before_abs, 0);
});

test('graze-then-recover (stays low): t_graze set, t_abs censored, 1 graze', () => {
  const R = build([
    [LO, 20],
    [HI, 10], // graze
    [LO, 300], // recovers and never absorbs
  ]);
  const r = labelSeries(R, DT, {}, 330);
  assert.equal(r.graze_censored, false);
  assert.ok(near(r.t_graze, 20), `t_graze=${r.t_graze}`);
  assert.equal(r.abs_censored, true);
  assert.equal(r.n_grazes_before_abs, 1);
});

test('graze-then-absorb: t_graze at graze, t_abs at second crossing, 1 graze', () => {
  const R = build([
    [LO, 20],
    [HI, 10], // graze 1
    [LO, 10], // recovery
    [HI, 200], // absorbs
  ]);
  const r = labelSeries(R, DT, {}, 240);
  assert.ok(near(r.t_graze, 20), `t_graze=${r.t_graze}`);
  assert.equal(r.abs_censored, false);
  assert.ok(near(r.t_abs, 40), `t_abs=${r.t_abs}`);
  assert.equal(r.n_grazes_before_abs, 1);
});

test('churn-then-absorb: 3 grazes then permanent merge', () => {
  const R = build([
    [LO, 10],
    [HI, 8], // cross 1  (t=10)
    [LO, 10], // recovery
    [HI, 8], // cross 2  (t=28)
    [LO, 10], // recovery
    [HI, 8], // cross 3  (t=46)
    [LO, 10], // recovery
    [HI, 200], // absorbs (t=64)
  ]);
  const r = labelSeries(R, DT, {}, 304);
  assert.ok(near(r.t_graze, 10), `t_graze=${r.t_graze}`);
  assert.ok(near(r.t_abs, 64), `t_abs=${r.t_abs}`);
  assert.equal(r.n_grazes_before_abs, 3);
});

test('shallow sub-θ dip is not recovery: absorbs at first crossing', () => {
  const R = build([
    [LO, 20],
    [HI, 30], // crossing
    [SHALLOW, 40], // dips below θ but never below recThresh ⇒ NOT recovery
    [HI, 120], // back up; whole post-crossing window has no recovery
  ]);
  const r = labelSeries(R, DT, {}, 210);
  assert.ok(near(r.t_graze, 20), `t_graze=${r.t_graze}`);
  assert.equal(r.abs_censored, false);
  assert.ok(near(r.t_abs, 20), `t_abs=${r.t_abs}`);
  assert.equal(r.n_grazes_before_abs, 0);
});

test('censoring: crossing too late for full T_v window ⇒ t_abs censored', () => {
  // Crossing confirmed at ~25s; only ~60s of no-recovery data before t_max=100
  // (< T_v=120) ⇒ cannot certify absorption.
  const R = build([
    [LO, 20],
    [HI, 80], // never recovers, but window < T_v
  ]);
  const r = labelSeries(R, DT, { T_v: 120 }, 100);
  assert.ok(near(r.t_graze, 20), `t_graze=${r.t_graze}`);
  assert.equal(r.abs_censored, true);
});

test('T_v sensitivity: shorter horizon certifies the same crossing as absorbed', () => {
  const R = build([
    [LO, 20],
    [HI, 80],
  ]);
  const short = labelSeries(R, DT, { T_v: 60 }, 100);
  assert.equal(short.abs_censored, false);
  assert.ok(near(short.t_abs, 20), `t_abs=${short.t_abs}`);
});

test('recThresh sensitivity: looser 0.75 ignores a shallow 0.78 dip', () => {
  const R = build([
    [LO, 20],
    [HI, 10], // crossing
    [0.78, 10], // dip: recovery under 0.80 but NOT under 0.75
    [HI, 200], // stays merged
  ]);
  const strict = labelSeries(R, DT, { recThresh: 0.8 }, 240); // 0.78<0.80 ⇒ graze
  const loose = labelSeries(R, DT, { recThresh: 0.75 }, 240); // 0.78>0.75 ⇒ no recovery
  assert.equal(strict.n_grazes_before_abs, 1);
  assert.ok(near(strict.t_abs, 40), `strict t_abs=${strict.t_abs}`);
  assert.equal(loose.n_grazes_before_abs, 0);
  assert.ok(near(loose.t_abs, 20), `loose t_abs=${loose.t_abs}`);
});

test('defaults are the documented absorption-grade constants', () => {
  assert.equal(DEFAULTS.theta, 0.85);
  assert.equal(DEFAULTS.W, 5.0);
  assert.equal(DEFAULTS.T_v, 120.0);
  assert.equal(DEFAULTS.recThresh, 0.8);
  assert.equal(DEFAULTS.recWin, 5.0);
});
