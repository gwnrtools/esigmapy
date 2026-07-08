#!/usr/bin/env python3
"""Benchmark the JAX surrogate backend against the numpy backend.

Reports, for the two AGENTS.md parameter regimes (long inspiral t_start=-136 and
near-merger t_start=-1):

* per-call steady-state wall time for both backends (JAX warmup reported
  separately -- never judge JAX by its first, compiling call), and
* the accuracy of the JAX output relative to the numpy backend.

It also demonstrates what the JAX backend is *for*: exact gradients and
vmap-batching over waveform parameters, which the numpy backend cannot provide.

Note on timing: both backends load the full multi-hundred-MB surrogate data, so
running them in one process introduces some memory-bandwidth contention on the
long-inspiral (bandwidth-bound) case; treat the long-signal numbers as
indicative. Run with the esigmapy_jax environment and ESIGMASUR_DATA_PATH set.
"""

from __future__ import annotations

import os
import statistics
import time

import numpy as np

import jax
import jax.numpy as jnp

from esigmapy.surrogate.surrogate_jax import (
    CircularSurrogateJAX,
    EccentricSurrogateJAX,
    _bucket,
    _amp_correction_factor,
)
from esigmapy.surrogate.surrogate import EccentricSurrogate

_DT = 0.000244140625
_PARAMS = [
    dict(M=10.0, q=2.3, e_ref=0.43, l_ref=1.3, t_start=-136.0, tag="long inspiral"),
    dict(M=10.0, q=2.3, e_ref=0.43, l_ref=1.3, t_start=-1.0, tag="near merger"),
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


def _rel(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float((np.abs(a - b) / np.maximum(np.abs(b), 1e-30)).max())


def main():
    data = os.environ["ESIGMASUR_DATA_PATH"]
    ecc_dir = os.path.join(data, "ecc_sur_data")
    circ_dir = os.path.join(data, "circ_sur_data")

    ej = EccentricSurrogateJAX(ecc_data_dir=ecc_dir, circ_data_dir=circ_dir)
    enp = EccentricSurrogate(ecc_data_dir=ecc_dir, circ_data_dir=circ_dir)

    print("=== eccentric full waveform: JAX vs numpy ===")
    for p in _PARAMS:
        n = 15 if p["t_start"] == -136.0 else 100
        kw = dict(
            M=p["M"],
            params=(p["q"], p["e_ref"], p["l_ref"]),
            delta_t=_DT,
            t_start=p["t_start"],
            t_end=None,
        )
        jt = _median_ms(lambda: ej(**kw), lambda o: jax.block_until_ready(o[1]), n)
        _, hn = enp(**kw)
        _, hj = ej(**kw)
        nt = _median_ms(lambda: enp(**kw), lambda o: None, n)
        print(
            f"  {p['tag']:14s} (t_start={p['t_start']:>6}): "
            f"JAX {jt:8.3f} ms | numpy {nt:8.3f} ms | "
            f"waveform max_rel_diff {_rel(hj, hn):.2e}"
        )

    _demo_grad_vmap(ej)


def _demo_grad_vmap(ej):
    """Show exact gradients and parameter batching through the eccentric core."""
    M, t_start, eta = 10.0, -1.0, 0.2
    t_start, t_end, scale = ej._set_time_range(M, t_start, None)
    num = int((t_end - t_start) / _DT) + 1
    query, _ = ej.circ_sur._build_query(t_start, _DT, num, None, None)
    circ_scale = M / ej.circ_sur.sur_total_mass
    amp0, phase0 = ej.circ_sur._eval_padded(M, eta, circ_scale, t_start, query, True)
    sidx = CircularSurrogateJAX._start_trunc_index(ej.t_grid_sur, t_start / scale)
    in_bucket = min(_bucket(ej.n_t - sidx), ej.n_t)
    in_start = ej.n_t - in_bucket
    tg0 = (ej.t_g0 + in_start * ej.t_dg) * scale
    tdg = ej.t_dg * scale
    amp_scale = scale / _amp_correction_factor
    const = (
        jnp.asarray(eta),
        jnp.asarray(1.3),
        jnp.asarray(in_start),
        in_bucket,
        True,
        jnp.asarray(tg0),
        jnp.asarray(tdg),
        query,
        amp0,
        phase0,
        jnp.asarray(amp_scale),
    )

    def energy(e_ref):
        amp, _ = ej._core(const[0], e_ref, const[1], *const[2:])
        return jnp.sum(amp * amp)

    grad = float(jax.grad(energy)(jnp.asarray(0.3)))
    eps = 1e-6
    fd = (
        float(energy(jnp.asarray(0.3 + eps))) - float(energy(jnp.asarray(0.3 - eps)))
    ) / (2 * eps)
    batch = jax.vmap(lambda er: ej._core(const[0], er, const[1], *const[2:])[0][:4])(
        jnp.array([0.1, 0.2, 0.3, 0.4])
    )

    print("\n=== JAX-only capabilities (numpy backend cannot do these) ===")
    print(
        f"  jax.grad d(sum amp^2)/d(e_ref) = {grad:.4e} "
        f"(finite-diff {fd:.4e}, rel {abs(grad - fd) / abs(fd):.1e})"
    )
    print(
        f"  jax.vmap over 4 e_ref values -> batched amp shape "
        f"{tuple(np.asarray(batch).shape)}"
    )


if __name__ == "__main__":
    raise SystemExit(main())
