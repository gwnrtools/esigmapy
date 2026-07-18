#!/usr/bin/env python3
"""Assess float32 (JAX default precision) for the JAX surrogate backend.

``surrogate_jax`` runs in float64 because importing ``TPI_jax`` enables JAX's
x64 mode. This decision tool measures what NOT forcing float64 would cost in
accuracy and buy in speed: ``--mode run`` evaluates a set of waveforms (short
and long, several parameter draws) in one precision and saves them; run it
once per precision (the x64 flag is process-wide, so the two precisions need
two processes), then ``--mode compare`` reports per-case amp/phase/mode
differences and waveform mismatches, the metric that decides usability.

Usage:
    python bench_jax_float32.py --mode run --precision f64 --out /tmp/f64.npz
    python bench_jax_float32.py --mode run --precision f32 --out /tmp/f32.npz
    python bench_jax_float32.py --mode compare /tmp/f64.npz /tmp/f32.npz

``--time`` in run mode also reports jitted single-call and vmapped-batch
timings for that precision.
"""

from __future__ import annotations

import argparse
import os
import statistics
import time

import numpy as np

DELTA_T = 2.44140625e-04  # 1/4096 s
M_TOT = 10.0
CASES = {  # name -> t_start (s); M=10 so -136 s is the full surrogate length
    "short": -1.0,
    "long": -136.0,
}
PARAM_DRAWS = 8
BATCH = 25


def _draws(n, seed=37):
    rng = np.random.default_rng(seed)
    return (
        rng.uniform(1.0, 6.0, n),
        rng.uniform(0.01, 0.43, n),
        rng.uniform(0.0, 2 * np.pi, n),
    )


def run(precision, out_path, do_time):
    import esigmapy.surrogate.surrogate_jax as sj  # imports TPI_jax, enables x64
    import jax

    if precision == "f32":
        # x64 was enabled by the TPI_jax import above; turn it off again BEFORE
        # the surrogate is constructed, so every device array (data tables and
        # traced values alike) is float32/int32.
        jax.config.update("jax_enable_x64", False)

    import jax.numpy as jnp

    data_dir = os.environ["ESIGMASUR_DATA_PATH"]
    sur = sj.EccentricSurrogateJAX(
        ecc_data_dir=os.path.join(data_dir, "ecc_sur_data"),
        circ_data_dir=os.path.join(data_dir, "circ_sur_data"),
    )

    q_arr, e_arr, l_arr = _draws(PARAM_DRAWS)
    saved = {"q": q_arr, "e": e_arr, "l": l_arr}
    for case, t_start in CASES.items():
        t_grid, fn = sur.parameter_space_evaluator(
            M=M_TOT, delta_t=DELTA_T, t_start=t_start, output="amp_phase"
        )
        fn_j = jax.jit(fn)
        amps, phases = [], []
        for i in range(PARAM_DRAWS):
            amp, phase = fn_j(q_arr[i], e_arr[i], l_arr[i])
            amps.append(np.asarray(amp, dtype=np.float64))
            phases.append(np.asarray(phase, dtype=np.float64))
        saved[f"{case}_t"] = np.asarray(t_grid, dtype=np.float64)
        saved[f"{case}_amp"] = np.asarray(amps)
        saved[f"{case}_phase"] = np.asarray(phases)

        if do_time:
            _, fn_m = sur.parameter_space_evaluator(
                M=M_TOT, delta_t=DELTA_T, t_start=t_start
            )
            single = jax.jit(fn_m)
            batched = jax.jit(jax.vmap(fn_m))
            qb, eb, lb = (
                jnp.asarray(a[:1].repeat(BATCH)) for a in (q_arr, e_arr, l_arr)
            )
            jax.block_until_ready(single(q_arr[0], e_arr[0], l_arr[0]))
            jax.block_until_ready(batched(qb, eb, lb))
            reps = 50 if case == "short" else 10
            ts = []
            for _ in range(reps):
                t0 = time.perf_counter()
                jax.block_until_ready(single(q_arr[0], e_arr[0], l_arr[0]))
                ts.append(time.perf_counter() - t0)
            t_single = statistics.median(ts)
            ts = []
            for _ in range(max(reps // 5, 3)):
                t0 = time.perf_counter()
                jax.block_until_ready(batched(qb, eb, lb))
                ts.append(time.perf_counter() - t0)
            t_batch = statistics.median(ts) / BATCH
            print(
                f"[{precision}] {case}: single {t_single*1e3:8.3f} ms/call, "
                f"batched(x{BATCH}) {t_batch*1e3:8.3f} ms/waveform"
            )

    np.savez_compressed(out_path, precision=precision, **saved)
    print(f"[{precision}] saved {out_path}")


def _mismatch(h_a, h_b):
    overlap = np.abs(np.vdot(h_a, h_b)) / np.sqrt(
        np.vdot(h_a, h_a).real * np.vdot(h_b, h_b).real
    )
    return 1.0 - overlap


def compare(path_ref, path_head):
    ref = np.load(path_ref)
    head = np.load(path_head)
    print(f"reference: {ref['precision']}, head: {head['precision']}")
    for case in CASES:
        a_r, p_r = ref[f"{case}_amp"], ref[f"{case}_phase"]
        a_h, p_h = head[f"{case}_amp"], head[f"{case}_phase"]
        amp_rel = np.max(np.abs(a_h - a_r) / np.max(np.abs(a_r), axis=1, keepdims=True))
        dphase = np.abs(p_h - p_r)
        mms = []
        for i in range(a_r.shape[0]):
            h_r = a_r[i] * np.exp(-1j * p_r[i])
            h_h = a_h[i] * np.exp(-1j * p_h[i])
            mms.append(_mismatch(h_r, h_h))
        mms = np.asarray(mms)
        print(
            f"{case}: max rel amp diff {amp_rel:.2e}, "
            f"max |dphase| {np.max(dphase):.2e} rad, "
            f"mismatch median {np.median(mms):.2e} / max {np.max(mms):.2e}"
        )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=("run", "compare"), required=True)
    p.add_argument("--precision", choices=("f64", "f32"))
    p.add_argument("--out")
    p.add_argument("--time", action="store_true")
    p.add_argument("paths", nargs="*")
    args = p.parse_args()
    if args.mode == "run":
        if not (args.precision and args.out):
            p.error("--mode run needs --precision and --out")
        run(args.precision, args.out, args.time)
    else:
        if len(args.paths) != 2:
            p.error("--mode compare needs two npz paths (ref head)")
        compare(*args.paths)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
