#!/usr/bin/env python3
"""Benchmark the JAX IMR pipeline against the NumPy IMR backend.

Reports:

* per-call steady-state wall time of the NumPy ``get_imr_esigmasur_mode``, the
  JAX host path (``IMRESIGMASurJAX.__call__``, eager kernels) and the jitted
  traced evaluator (``imr_parameter_space_evaluator``), for a short and a
  full-length inspiral (JAX warmup/compile reported separately);
* vmap-batched evaluator throughput per waveform vs batch size (CPU, and GPU
  when available);
* the cost multiple of a full parameter-space gradient over a forward eval;
* the (2,2) waveform mismatch between the two backends (accuracy context).

Run in the esigmapy_jax environment with ESIGMASUR_DATA_PATH set; add
XLA_PYTHON_CLIENT_PREALLOCATE=false on small GPUs. JAX_PLATFORMS=cpu forces a
CPU-only run.
"""

from __future__ import annotations

import os
import statistics
import time

import numpy as np

from esigmapy.surrogate.generator_jax import IMRESIGMASurJAX

import jax
import jax.numpy as jnp

from esigmapy.surrogate.generator import get_imr_esigmasur_mode

_DT = 0.000244140625
_PARAMS = [
    dict(m1=41.4, m2=18.6, e_ref=0.25, l_ref=1.3, t_start=-40.0, tag="short (40 s)"),
    dict(m1=41.4, m2=18.6, e_ref=0.25, l_ref=1.3, t_start=None, tag="full length"),
]


def _median_ms(fn, block, n):
    fn()  # warmup / compile
    samples = []
    for _ in range(n):
        t0 = time.perf_counter()
        out = fn()
        block(out)
        samples.append((time.perf_counter() - t0) * 1e3)
    return statistics.median(samples)


def _mismatch(h_np, h_j):
    try:
        from pycbc.filter import match
        from pycbc.types import TimeSeries
    except Exception:
        return float("nan")
    ln = 2 ** int(np.ceil(np.log2(max(len(h_np), len(h_j)))))
    tsn = TimeSeries(np.real(np.pad(h_np, (0, ln - len(h_np)))), delta_t=_DT)
    tsj = TimeSeries(np.real(np.pad(h_j, (0, ln - len(h_j)))), delta_t=_DT)
    mm, _ = match(tsn, tsj)
    return 1.0 - mm


def main():
    data = os.environ["ESIGMASUR_DATA_PATH"]
    imr = IMRESIGMASurJAX(
        os.path.join(data, "ecc_sur_data"), os.path.join(data, "circ_sur_data")
    )

    print("=== single evaluation: NumPy driver vs JAX host path vs jitted fn ===")
    for p in _PARAMS:
        n = 10 if p["t_start"] == -40.0 else 3
        M, q = p["m1"] + p["m2"], p["m1"] / p["m2"]
        np_kw = dict(
            mass1=p["m1"],
            mass2=p["m2"],
            delta_t=_DT,
            reference_eccentricity=p["e_ref"],
            reference_mean_anomaly=p["l_ref"],
            t_start=p["t_start"],
        )
        j_args = (M, (q, p["e_ref"], p["l_ref"]), _DT)
        nt = _median_ms(lambda: get_imr_esigmasur_mode(**np_kw), lambda o: None, n)
        jt = _median_ms(
            lambda: imr(*j_args, t_start=p["t_start"]),
            lambda o: None,
            n,
        )
        times, fn = imr.imr_parameter_space_evaluator(M, _DT, t_start=p["t_start"])
        jfn = jax.jit(fn)
        t0 = time.perf_counter()
        jax.block_until_ready(jfn(q, p["e_ref"], p["l_ref"]))
        compile_s = time.perf_counter() - t0
        et = _median_ms(
            lambda: jfn(q, p["e_ref"], p["l_ref"]),
            lambda o: jax.block_until_ready(o[0]),
            n,
        )
        modes_np = get_imr_esigmasur_mode(**np_kw)
        _, modes_j = imr(*j_args, t_start=p["t_start"])
        mm = _mismatch(np.asarray(modes_np[(2, 2)].data), np.asarray(modes_j[(2, 2)]))
        print(
            f"[{p['tag']}] numpy {nt:8.1f} ms | jax host {jt:8.1f} ms | "
            f"jit fn {et:7.1f} ms (compile {compile_s:.1f} s) | 1-M {mm:.1e}"
        )

    print("\n=== vmap-batched evaluator, ms per waveform ===")
    p = _PARAMS[0]
    M = p["m1"] + p["m2"]
    times, fn = imr.imr_parameter_space_evaluator(M, _DT, t_start=p["t_start"])
    devices = {"cpu": jax.devices("cpu")[0]}
    try:
        devices["gpu"] = jax.devices("gpu")[0]
    except Exception:
        pass
    rng = np.random.default_rng(7)
    for dev_name, dev in devices.items():
        with jax.default_device(dev):
            vfn = jax.jit(jax.vmap(fn))
            for B in (1, 4, 16):
                qs = jnp.asarray(rng.uniform(1.5, 5.0, B))
                es = jnp.asarray(rng.uniform(0.1, 0.4, B))
                ls = jnp.asarray(rng.uniform(0.0, 2 * np.pi, B))
                bt = _median_ms(
                    lambda: vfn(qs, es, ls),
                    lambda o: jax.block_until_ready(o[0]),
                    5,
                )
                print(f"[{dev_name}] B={B:3d}: {bt / B:8.1f} ms/waveform")

    print("\n=== gradient cost (jitted, CPU default device) ===")
    q, e0, l0 = p["m1"] / p["m2"], p["e_ref"], p["l_ref"]

    def loss(qq, ee, ll):
        h, _ = fn(qq, ee, ll)
        return jnp.sum(jnp.abs(h) ** 2) * 1e36

    jloss = jax.jit(loss)
    jgrad = jax.jit(jax.grad(loss, argnums=(0, 1, 2)))
    ft = _median_ms(lambda: jloss(q, e0, l0), jax.block_until_ready, 10)
    gt = _median_ms(lambda: jgrad(q, e0, l0), lambda o: jax.block_until_ready(o[0]), 10)
    print(f"forward {ft:.1f} ms | grad(q,e,l) {gt:.1f} ms | ratio {gt / ft:.1f}x")


if __name__ == "__main__":
    main()
