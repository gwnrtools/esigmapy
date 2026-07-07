#!/usr/bin/env python3
"""Standalone benchmark for the fused final-stage kernel (Part A step A1).

Motivation (from A0 evidence): replacing ``np.interp`` with a standalone
index-arithmetic interpolation is NOT a win -- ``np.interp`` is already
well-optimized C, and manual uniform-grid variants are the same speed or
slower. The only lever left on the final uniform-time-grid stage is *fusion*:
doing the two interpolations, the amplitude scaling, the optional initial-phase
subtraction, and the complex assembly ``amp * exp(-1j*phase)`` in a single pass
over the output grid, so the large output array is written once and the native
tables are read once, instead of materializing several full-size temporaries.

This script measures the CURRENT pipeline

    amp   = np.interp(new_t_grid, t_grid_sur, amp_native)
    phase = np.interp(new_t_grid, t_grid_sur, phase_native)
    amp  *= amp_scale
    phase -= phase[0]            # optional
    mode  = amp * np.exp(-1j*phase)

against a single fused numba kernel, at the two surrogate-relevant regimes.
Both the native grid and the output grid are uniform (verified for the real
surrogate data), so the kernel indexes by arithmetic. It reproduces
``np.interp`` boundary clamping exactly.

Accuracy is checked against the current pipeline; the A1 gate is
max_abs_diff <= 1e-9 / max_rel_diff <= 1e-8 (fastmath carve-out, since the
kernel uses fastmath sincos) with the standard 1e-12/1e-10 reported too.
"""

from __future__ import annotations

import argparse
import time

import numpy as np
from numba import njit


# --- current pipeline -------------------------------------------------------


def _current_mode(new_t_grid, t_grid_sur, amp_native, phase_native, amp_scale,
                  remove_phase0):
    amp = np.interp(new_t_grid, t_grid_sur, amp_native)
    phase = np.interp(new_t_grid, t_grid_sur, phase_native)
    amp = amp * amp_scale
    if remove_phase0:
        phase = phase - phase[0]
    return amp * np.exp(-1j * phase)


# --- fused kernel -----------------------------------------------------------


@njit(fastmath=True, cache=True)
def _fused_interp_mode_uniform(g0, dg, amp_table, phase_table, q0, dq, n_out,
                               amp_scale, phase0):
    """One-pass uniform-grid interp + scale + complex assembly.

    g0, dg  : first value and spacing of the native (source) grid.
    amp_table, phase_table : native tables (same length, uniform grid).
    q0, dq  : first value and spacing of the uniform output/query grid.
    phase0  : value subtracted from every phase (0.0 to disable).

    Reproduces np.interp clamping: queries at/below g0 take table[0], queries
    at/above the last node take table[-1].
    """
    ntab = amp_table.shape[0]
    out = np.empty(n_out, dtype=np.complex128)
    inv_dg = 1.0 / dg
    last = ntab - 1
    for i in range(n_out):
        tq = q0 + i * dq
        pos = (tq - g0) * inv_dg
        if pos <= 0.0:
            a = amp_table[0]
            ph = phase_table[0]
        elif pos >= last:
            a = amp_table[last]
            ph = phase_table[last]
        else:
            idx = int(pos)
            w = pos - idx
            a = amp_table[idx] + w * (amp_table[idx + 1] - amp_table[idx])
            ph = phase_table[idx] + w * (phase_table[idx + 1] - phase_table[idx])
        a *= amp_scale
        ph -= phase0
        out[i] = a * (np.cos(ph) - 1j * np.sin(ph))
    return out


@njit(fastmath=True, cache=True)
def _fused_interp_amp_phase_uniform(g0, dg, amp_table, phase_table, q0, dq,
                                    n_out, amp_scale, phase0):
    ntab = amp_table.shape[0]
    amp = np.empty(n_out, dtype=np.float64)
    phase = np.empty(n_out, dtype=np.float64)
    inv_dg = 1.0 / dg
    last = ntab - 1
    for i in range(n_out):
        tq = q0 + i * dq
        pos = (tq - g0) * inv_dg
        if pos <= 0.0:
            a = amp_table[0]
            ph = phase_table[0]
        elif pos >= last:
            a = amp_table[last]
            ph = phase_table[last]
        else:
            idx = int(pos)
            w = pos - idx
            a = amp_table[idx] + w * (amp_table[idx + 1] - amp_table[idx])
            ph = phase_table[idx] + w * (phase_table[idx + 1] - phase_table[idx])
        amp[i] = a * amp_scale
        phase[i] = ph - phase0
    return amp, phase


def _bench_case(n_tab, n_out, seed, repetitions):
    rng = np.random.default_rng(seed)
    g0, dg = -136.0, 0.0004925508729
    t_grid_sur = g0 + dg * np.arange(n_tab)
    amp_native = np.abs(rng.standard_normal(n_tab)) + 0.1
    phase_native = np.cumsum(np.abs(rng.standard_normal(n_tab))) * 1e-3

    # Uniform output grid spanning slightly inside the native range.
    q0 = g0 + 3.3 * dg
    dq = (t_grid_sur[-1] - q0) / (n_out - 1)
    new_t_grid = q0 + dq * np.arange(n_out)
    amp_scale = 1.234
    remove_phase0 = True
    phase0 = phase_native[0] if False else 0.0  # phase0 for kernel set below

    ref = _current_mode(new_t_grid, t_grid_sur, amp_native, phase_native,
                        amp_scale, remove_phase0)

    # For the fused kernel, phase0 must equal the interpolated phase at sample 0.
    ph0 = np.interp(new_t_grid[0], t_grid_sur, phase_native) if remove_phase0 else 0.0

    # warmup / compile
    got = _fused_interp_mode_uniform(g0, dg, amp_native, phase_native, q0, dq,
                                    n_out, amp_scale, ph0)

    diff = np.abs(ref - got)
    denom = np.maximum(np.abs(ref), 1e-30)
    max_abs = float(diff.max())
    max_rel = float((diff / denom).max())

    def time_fn(fn, *a):
        t0 = time.perf_counter()
        for _ in range(repetitions):
            fn(*a)
        return (time.perf_counter() - t0) / repetitions

    t_cur = time_fn(_current_mode, new_t_grid, t_grid_sur, amp_native,
                    phase_native, amp_scale, remove_phase0)
    t_fused = time_fn(_fused_interp_mode_uniform, g0, dg, amp_native,
                     phase_native, q0, dq, n_out, amp_scale, ph0)

    print(f"n_tab={n_tab} n_out={n_out}")
    print(f"  current pipeline  avg_sec={t_cur:.6e}")
    print(f"  fused kernel      avg_sec={t_fused:.6e}  speedup={t_cur/t_fused:.2f}x")
    print(f"  max_abs_diff={max_abs:.3e}  max_rel_diff={max_rel:.3e}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--repetitions", type=int, default=50)
    p.add_argument("--seed", type=int, default=7)
    args = p.parse_args()
    # near-merger: small output; long inspiral: ~1.2M output
    _bench_case(n_tab=280_000, n_out=5_000, seed=args.seed, repetitions=args.repetitions)
    _bench_case(n_tab=1_214_525, n_out=1_100_000, seed=args.seed, repetitions=args.repetitions)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
