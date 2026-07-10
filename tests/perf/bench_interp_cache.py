#!/usr/bin/env python3
"""Benchmark one-time caching of linear-interpolant index/weight data.

Decision tool for the question: since the surrogate's native time /
shifted-mean-anomaly grids are fixed and only the interpolated VALUES change
per call, does precomputing the (index, weight) pair per output sample once
and reusing it beat the current per-sample index arithmetic
(`pos = (q0 + i*dq - g0) / dg; idx = int(pos); w = pos - idx`) inside the
fused kernels?

Scope notes (what a cache could even apply to):
* The t-grid interpolation queries are the uniform output grid; its idx/w
  depend on (M, t_start, delta_t, num_samples) -- reusable only across calls
  that share the full grid configuration.
* The l-grid interpolation queries are `l_s`, itself the output of the t-grid
  interpolation of `lt_relation`, which depends on (q, e_ref, l_ref) -- NOT
  cacheable across calls under any configuration.
So the cache can at most replace the t-grid index arithmetic. This script
measures (a) the isolated interp stage, (b) the full fused circular kernel,
and (c) the full fused eccentric kernel, each in "current" (index arithmetic)
vs "cached" (precomputed idx/w arrays) form, at the short/near-merger and
long-inspiral output sizes, plus the one-time cache build cost.
"""

from __future__ import annotations

import argparse
import time

import numpy as np
from numba import njit


# --- cache construction (the one-time cost) ---------------------------------


@njit(fastmath=True, cache=True)
def build_cache(g0, dg, ntab, q0, dq, n_out):
    """Precompute clamped (idx, w) per output sample for a uniform table grid.

    Encodes np.interp boundary clamping the same way the fused kernels do:
    out-of-range queries get w=0 at idx=0 or idx=last-1, w=1.0."""
    last = ntab - 1
    inv_dg = 1.0 / dg
    idx = np.empty(n_out, dtype=np.int64)
    w = np.empty(n_out, dtype=np.float64)
    for i in range(n_out):
        pos = (q0 + i * dq - g0) * inv_dg
        if pos <= 0.0:
            idx[i] = 0
            w[i] = 0.0
        elif pos >= last:
            idx[i] = last - 1
            w[i] = 1.0
        else:
            j = int(pos)
            idx[i] = j
            w[i] = pos - j
    return idx, w


# --- (a) isolated interp stage ----------------------------------------------


@njit(fastmath=True, cache=True)
def interp_current(g0, dg, table, q0, dq, n_out):
    last = table.shape[0] - 1
    inv_dg = 1.0 / dg
    out = np.empty(n_out, dtype=np.float64)
    for i in range(n_out):
        pos = (q0 + i * dq - g0) * inv_dg
        if pos <= 0.0:
            out[i] = table[0]
        elif pos >= last:
            out[i] = table[last]
        else:
            j = int(pos)
            w = pos - j
            out[i] = table[j] + w * (table[j + 1] - table[j])
    return out


@njit(fastmath=True, cache=True)
def interp_cached(idx, w, table, n_out):
    out = np.empty(n_out, dtype=np.float64)
    for i in range(n_out):
        j = idx[i]
        wi = w[i]
        out[i] = table[j] + wi * (table[j + 1] - table[j])
    return out


# Control variant: same on-the-fly index arithmetic as "current" but with the
# np.interp boundary clamping done branchlessly (min/max), no cached data.
# Distinguishes "the cache data helps" from "removing the branches helps".


@njit(fastmath=True, cache=True)
def interp_branchless(g0, dg, table, q0, dq, n_out):
    last = table.shape[0] - 1
    inv_dg = 1.0 / dg
    out = np.empty(n_out, dtype=np.float64)
    for i in range(n_out):
        pos = max((q0 + i * dq - g0) * inv_dg, 0.0)
        j = min(int(pos), last - 1)
        w = min(pos - j, 1.0)
        out[i] = table[j] + w * (table[j + 1] - table[j])
    return out


# --- (b) fused circular mode kernel (as in _fused_interp_mode_uniform) ------


@njit(fastmath=True, cache=True)
def circ_mode_current(g0, dg, amp_table, phase_table, q0, dq, n_out, amp_scale):
    last = amp_table.shape[0] - 1
    inv_dg = 1.0 / dg
    out = np.empty(n_out, dtype=np.complex128)
    for i in range(n_out):
        pos = (q0 + i * dq - g0) * inv_dg
        if pos <= 0.0:
            a = amp_table[0]
            ph = phase_table[0]
        elif pos >= last:
            a = amp_table[last]
            ph = phase_table[last]
        else:
            j = int(pos)
            w = pos - j
            a = amp_table[j] + w * (amp_table[j + 1] - amp_table[j])
            ph = phase_table[j] + w * (phase_table[j + 1] - phase_table[j])
        a *= amp_scale
        out[i] = a * (np.cos(ph) - 1j * np.sin(ph))
    return out


@njit(fastmath=True, cache=True)
def circ_mode_cached(idx, w, amp_table, phase_table, n_out, amp_scale):
    out = np.empty(n_out, dtype=np.complex128)
    for i in range(n_out):
        j = idx[i]
        wi = w[i]
        a = amp_table[j] + wi * (amp_table[j + 1] - amp_table[j])
        ph = phase_table[j] + wi * (phase_table[j + 1] - phase_table[j])
        a *= amp_scale
        out[i] = a * (np.cos(ph) - 1j * np.sin(ph))
    return out


@njit(fastmath=True, cache=True)
def circ_mode_branchless(g0, dg, amp_table, phase_table, q0, dq, n_out, amp_scale):
    last = amp_table.shape[0] - 1
    inv_dg = 1.0 / dg
    out = np.empty(n_out, dtype=np.complex128)
    for i in range(n_out):
        pos = max((q0 + i * dq - g0) * inv_dg, 0.0)
        j = min(int(pos), last - 1)
        w = min(pos - j, 1.0)
        a = amp_table[j] + w * (amp_table[j + 1] - amp_table[j])
        ph = phase_table[j] + w * (phase_table[j + 1] - phase_table[j])
        a *= amp_scale
        out[i] = a * (np.cos(ph) - 1j * np.sin(ph))
    return out


# --- (c) fused eccentric kernel (as in _fused_ecc_mode_uniform) --------------
# Only the t-grid (lt_relation) lookup can use the cache; the l-grid lookups
# depend on the per-call l_s values and stay index-arithmetic in both variants.


@njit(fastmath=True, cache=True)
def ecc_mode_current(
    tg0, tdg, lt_relation, q0, dq, n_out, lg0, ldg, damp, dphase, rcp,
    amp0, phase0, amp_scale,
):
    last_t = lt_relation.shape[0] - 1
    inv_tdg = 1.0 / tdg
    last_l = damp.shape[0] - 1
    inv_ldg = 1.0 / ldg
    out = np.empty(n_out, dtype=np.complex128)
    for i in range(n_out):
        pos = (q0 + i * dq - tg0) * inv_tdg
        if pos <= 0.0:
            ls = lt_relation[0]
        elif pos >= last_t:
            ls = lt_relation[last_t]
        else:
            j = int(pos)
            w = pos - j
            ls = lt_relation[j] + w * (lt_relation[j + 1] - lt_relation[j])
        lpos = (ls - lg0) * inv_ldg
        if lpos <= 0.0:
            dA = damp[0]
            dP = dphase[0]
            rcpi = rcp[0]
        elif lpos >= last_l:
            dA = damp[last_l]
            dP = dphase[last_l]
            rcpi = rcp[last_l]
        else:
            j = int(lpos)
            w = lpos - j
            dA = damp[j] + w * (damp[j + 1] - damp[j])
            dP = dphase[j] + w * (dphase[j + 1] - dphase[j])
            rcpi = rcp[j] + w * (rcp[j + 1] - rcp[j])
        amp = amp0[i] + amp_scale * dA
        ph = phase0[i] + rcpi + dP
        out[i] = amp * (np.cos(ph) - 1j * np.sin(ph))
    return out


@njit(fastmath=True, cache=True)
def ecc_mode_cached_t(
    t_idx, t_w, lt_relation, n_out, lg0, ldg, damp, dphase, rcp,
    amp0, phase0, amp_scale,
):
    last_l = damp.shape[0] - 1
    inv_ldg = 1.0 / ldg
    out = np.empty(n_out, dtype=np.complex128)
    for i in range(n_out):
        j = t_idx[i]
        wi = t_w[i]
        ls = lt_relation[j] + wi * (lt_relation[j + 1] - lt_relation[j])
        lpos = (ls - lg0) * inv_ldg
        if lpos <= 0.0:
            dA = damp[0]
            dP = dphase[0]
            rcpi = rcp[0]
        elif lpos >= last_l:
            dA = damp[last_l]
            dP = dphase[last_l]
            rcpi = rcp[last_l]
        else:
            j = int(lpos)
            w = lpos - j
            dA = damp[j] + w * (damp[j + 1] - damp[j])
            dP = dphase[j] + w * (dphase[j + 1] - dphase[j])
            rcpi = rcp[j] + w * (rcp[j + 1] - rcp[j])
        amp = amp0[i] + amp_scale * dA
        ph = phase0[i] + rcpi + dP
        out[i] = amp * (np.cos(ph) - 1j * np.sin(ph))
    return out


@njit(fastmath=True, cache=True)
def ecc_mode_branchless(
    tg0, tdg, lt_relation, q0, dq, n_out, lg0, ldg, damp, dphase, rcp,
    amp0, phase0, amp_scale,
):
    last_t = lt_relation.shape[0] - 1
    inv_tdg = 1.0 / tdg
    last_l = damp.shape[0] - 1
    inv_ldg = 1.0 / ldg
    out = np.empty(n_out, dtype=np.complex128)
    for i in range(n_out):
        pos = max((q0 + i * dq - tg0) * inv_tdg, 0.0)
        j = min(int(pos), last_t - 1)
        w = min(pos - j, 1.0)
        ls = lt_relation[j] + w * (lt_relation[j + 1] - lt_relation[j])
        lpos = max((ls - lg0) * inv_ldg, 0.0)
        j = min(int(lpos), last_l - 1)
        w = min(lpos - j, 1.0)
        dA = damp[j] + w * (damp[j + 1] - damp[j])
        dP = dphase[j] + w * (dphase[j + 1] - dphase[j])
        rcpi = rcp[j] + w * (rcp[j + 1] - rcp[j])
        amp = amp0[i] + amp_scale * dA
        ph = phase0[i] + rcpi + dP
        out[i] = amp * (np.cos(ph) - 1j * np.sin(ph))
    return out


def _time(fn, reps):
    fn()  # warm (JIT + caches)
    fn()
    t0 = time.perf_counter()
    for _ in range(reps):
        fn()
    return (time.perf_counter() - t0) / reps


def _max_abs(a, b):
    return float(np.max(np.abs(a - b)))


def run_case(name, n_out, n_tab_t, n_tab_l, reps):
    rng = np.random.default_rng(1234)
    # Uniform native t grid spanning the query range with margin (start buffer
    # of ~5 samples, as the truncation index leaves).
    tg0, t_end = -1.0, 0.0
    tdg = (t_end - tg0) / (n_tab_t - 1) * (1.0 + 12.0 / n_tab_t)
    q0 = tg0 + 6.0 * tdg
    dq = (t_end - q0) / n_out * 0.999
    amp_table = rng.random(n_tab_t) + 0.5
    phase_table = np.linspace(0.0, 6.0e4, n_tab_t) + rng.random(n_tab_t)
    # l-grid pieces: lt_relation maps t->l monotonically over the l table span
    lg0, ldg = 0.0, 1.0 / (n_tab_l - 1)
    lt_relation = np.linspace(lg0, lg0 + (n_tab_l - 1) * ldg, n_tab_t)
    damp = rng.random(n_tab_l) * 0.01
    dphase = rng.random(n_tab_l) * 0.1
    rcp = rng.random(n_tab_l) * 0.1
    amp0 = rng.random(n_out) + 0.5
    phase0 = np.linspace(0.0, 6.0e4, n_out)
    amp_scale = 1.7

    idx, w = build_cache(tg0, tdg, n_tab_t, q0, dq, n_out)
    t_build = _time(lambda: build_cache(tg0, tdg, n_tab_t, q0, dq, n_out), reps)

    rows = []

    a = interp_current(tg0, tdg, amp_table, q0, dq, n_out)
    b = interp_cached(idx, w, amp_table, n_out)
    c = interp_branchless(tg0, tdg, amp_table, q0, dq, n_out)
    rows.append(("interp-only current   ", _time(
        lambda: interp_current(tg0, tdg, amp_table, q0, dq, n_out), reps), None))
    rows.append(("interp-only cached    ", _time(
        lambda: interp_cached(idx, w, amp_table, n_out), reps), _max_abs(a, b)))
    rows.append(("interp-only branchless", _time(
        lambda: interp_branchless(tg0, tdg, amp_table, q0, dq, n_out), reps),
        _max_abs(a, c)))

    a = circ_mode_current(tg0, tdg, amp_table, phase_table, q0, dq, n_out, amp_scale)
    b = circ_mode_cached(idx, w, amp_table, phase_table, n_out, amp_scale)
    c = circ_mode_branchless(
        tg0, tdg, amp_table, phase_table, q0, dq, n_out, amp_scale)
    rows.append(("circ fused  current   ", _time(
        lambda: circ_mode_current(
            tg0, tdg, amp_table, phase_table, q0, dq, n_out, amp_scale), reps), None))
    rows.append(("circ fused  cached    ", _time(
        lambda: circ_mode_cached(
            idx, w, amp_table, phase_table, n_out, amp_scale), reps), _max_abs(a, b)))
    rows.append(("circ fused  branchless", _time(
        lambda: circ_mode_branchless(
            tg0, tdg, amp_table, phase_table, q0, dq, n_out, amp_scale), reps),
        _max_abs(a, c)))

    a = ecc_mode_current(
        tg0, tdg, lt_relation, q0, dq, n_out, lg0, ldg, damp, dphase, rcp,
        amp0, phase0, amp_scale)
    b = ecc_mode_cached_t(
        idx, w, lt_relation, n_out, lg0, ldg, damp, dphase, rcp,
        amp0, phase0, amp_scale)
    c = ecc_mode_branchless(
        tg0, tdg, lt_relation, q0, dq, n_out, lg0, ldg, damp, dphase, rcp,
        amp0, phase0, amp_scale)
    rows.append(("ecc fused   current   ", _time(
        lambda: ecc_mode_current(
            tg0, tdg, lt_relation, q0, dq, n_out, lg0, ldg, damp, dphase, rcp,
            amp0, phase0, amp_scale), reps), None))
    rows.append(("ecc fused   cached    ", _time(
        lambda: ecc_mode_cached_t(
            idx, w, lt_relation, n_out, lg0, ldg, damp, dphase, rcp,
            amp0, phase0, amp_scale), reps), _max_abs(a, b)))
    rows.append(("ecc fused   branchless", _time(
        lambda: ecc_mode_branchless(
            tg0, tdg, lt_relation, q0, dq, n_out, lg0, ldg, damp, dphase, rcp,
            amp0, phase0, amp_scale), reps), _max_abs(a, c)))

    print(f"\n{name}: n_out={n_out}, n_tab_t={n_tab_t}, n_tab_l={n_tab_l}, reps={reps}")
    print(f"  cache build (one-time): {t_build*1e3:9.4f} ms, "
          f"size {(idx.nbytes + w.nbytes)/1e6:.1f} MB")
    for label, sec, diff in rows:
        extra = "" if diff is None else f"  max_abs_diff={diff:.3e}"
        print(f"  {label}: {sec*1e3:9.4f} ms{extra}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--reps-small", type=int, default=200)
    p.add_argument("--reps-large", type=int, default=30)
    args = p.parse_args()
    # short/near-merger regime (t_start=-1.0 at M=10, 4096 Hz)
    run_case("short", n_out=4_097, n_tab_t=9_000, n_tab_l=1_000, reps=args.reps_small)
    # long-inspiral regime (t_start=-136.0 at M=10, 4096 Hz); table sizes match
    # the circular (1.21M) / eccentric (277k) native t grids and l grid (125k)
    run_case("long-circ", n_out=557_057, n_tab_t=1_214_525, n_tab_l=125_000,
             reps=args.reps_large)
    run_case("long-ecc", n_out=557_057, n_tab_t=277_000, n_tab_l=125_000,
             reps=args.reps_large)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
