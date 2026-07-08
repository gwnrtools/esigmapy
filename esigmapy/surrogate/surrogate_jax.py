# Copyright (C) 2026 Akash Maurya

"""JAX implementation of the ESIGMASur surrogate.

This is a JAX port of ``surrogate.py`` built on the JAX backend of TPI
(``TPI_jax``). It targets the waveform paths (uniform output grid: ``delta_t``
with ``t_start``/``t_end``; full mode and ``return_amp_phase_only``); a
user-supplied ``times`` array and ``return_orbital_variables`` are supported as
secondary paths.

Design notes
------------
* ``import TPI_jax`` (below) enables JAX 64-bit mode at import time, so this
  module must be imported before any other JAX array is created in the process.
* Dynamic array shapes (truncated native grids, output sample count) are handled
  by *bucketing*: lengths are rounded up (keeping the top few significant bits)
  and the input truncation start index is moved down to include extra real grid
  samples, so a bounded number of jit signatures are compiled and reused across
  calls with different ``M``/``t_start``/``delta_t``.
* Heavy work runs inside a small number of ``jax.jit`` cores; parameter-space
  evaluation uses ``TP_Interpolant_ND_Vector`` (shared-node pieces) and a
  ``vmap`` over stacked per-fit splines (res_amp/res_phase). The cores are pure
  functions of scalars and device constants, so ``vmap``/``grad`` can be layered
  on later.
"""

import os
import sys

import numpy as np
import h5py
import scipy.interpolate as si


def _import_tpi_jax():
    """Import the JAX backend of TPI (``TPI_jax``).

    The compiled ``TPI`` package installs without ``TPI_jax.py``; the JAX
    backend lives in the TPI source tree. Locate it via ``$TPI_JAX_PATH`` or the
    ``TPI`` directory sitting beside the esigmapy repository, then import."""
    try:
        import TPI_jax  # noqa: F401

        return TPI_jax
    except ModuleNotFoundError:
        pass
    candidates = []
    env = os.environ.get("TPI_JAX_PATH")
    if env:
        candidates.append(env)
    repo_root = os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    )
    candidates.append(os.path.join(os.path.dirname(repo_root), "TPI"))
    for cand in candidates:
        if os.path.isfile(os.path.join(cand, "TPI_jax.py")):
            sys.path.insert(0, cand)
            break
    import TPI_jax  # noqa: F401

    return TPI_jax


TPI_jax = _import_tpi_jax()  # enables jax x64 at import; before any jax array
import jax
import jax.numpy as jnp
from functools import partial

from lal import SpinWeightedSphericalHarmonic
from pycbc.conversions import eta_from_q

_amp_correction_factor = float(np.real(SpinWeightedSphericalHarmonic(0, 0, -2, 2, 2)))
_TWO_PI = 2.0 * np.pi


def _bucket(n):
    """Round ``n`` up to a bucket size, keeping the top 3 significant bits.

    This bounds the number of distinct jit signatures (~4 buckets per octave,
    so a few dozen over the full grid range) while capping padding waste at
    ~25% -- far less than the ~2x worst case of power-of-two rounding, which
    matters because the padded length sets the amount of real gemv/interp work.
    """
    n = int(n)
    if n <= 8:
        return 8
    step = 1 << (n.bit_length() - 3)
    return ((n + step - 1) // step) * step


# --- device-side spline / interpolation primitives --------------------------


def _cubic_bspline_vec_eval(knots, coeffs, x):
    """Evaluate a 1D cubic (vector-valued) B-spline at scalar ``x``.

    knots  : (nt,) knot vector (clamped ends, as in a scipy tck).
    coeffs : (ncoef, m) control points (m value components).
    Returns (m,). Cox-de Boor recurrence for the four active cubic basis
    functions; matches ``scipy.interpolate.BSpline`` for in-range ``x`` and
    clamps to the boundary segment outside the knot support.
    """
    nt = knots.shape[0]
    mu = jnp.clip(jnp.searchsorted(knots, x, side="right") - 1, 3, nt - 5)
    k = jax.lax.dynamic_slice(knots, (mu - 2,), (6,))
    l1 = x - k[2]
    l2 = x - k[1]
    l3 = x - k[0]
    r1 = k[3] - x
    r2 = k[4] - x
    r3 = k[5] - x

    def sd(num, den):
        return jnp.where(den != 0.0, num / den, 0.0)

    b0 = sd(r1, r1 + l1)
    b1 = sd(l1, r1 + l1)
    q0 = r1 * sd(b0, r1 + l2)
    q1 = l2 * sd(b0, r1 + l2) + r2 * sd(b1, r2 + l1)
    q2 = l1 * sd(b1, r2 + l1)
    d0 = r1 * sd(q0, r1 + l3)
    d1 = l3 * sd(q0, r1 + l3) + r2 * sd(q1, r2 + l2)
    d2 = l2 * sd(q1, r2 + l2) + r3 * sd(q2, r3 + l1)
    d3 = l1 * sd(q2, r3 + l1)
    basis = jnp.stack((d0, d1, d2, d3))
    zero = jnp.zeros((), mu.dtype)
    block = jax.lax.dynamic_slice(coeffs, (mu - 3, zero), (4, coeffs.shape[1]))
    return basis @ block


def _interp_uniform(values, g0, dg, query):
    """Linear interpolation of ``values`` (on a uniform grid starting at ``g0``
    with spacing ``dg``) at ``query`` points, replicating np.interp clamping."""
    last = values.shape[0] - 1
    pos = (query - g0) / dg
    idx = jnp.clip(jnp.floor(pos).astype(jnp.int64), 0, last)
    idx1 = jnp.minimum(idx + 1, last)
    w = jnp.clip(pos - idx, 0.0, 1.0)
    v0 = jnp.take(values, idx)
    v1 = jnp.take(values, idx1)
    return v0 * (1.0 - w) + v1 * w


# --- circular surrogate ------------------------------------------------------


class CircularSurrogateJAX:
    """JAX port of ``CircularSurrogate`` (waveform paths)."""

    def __init__(self, data_dir):
        self.sur_dir = data_dir
        self.data_piece_names = ["amp", "phase"]
        self._load()

    # -- loading --

    def _get_metadata(self, keys):
        filename = os.path.join(self.sur_dir, "surrogate_metadata.hdf")
        with h5py.File(filename, "r") as f:
            return [f[k][()] for k in keys]

    def _load(self):
        self.sur_total_mass, t_grid_sur = self._get_metadata(["M", "t_grid_sur"])
        self.sur_total_mass = float(self.sur_total_mass)
        t_grid_sur = np.asarray(t_grid_sur, dtype=np.float64)
        self.t_grid_sur = t_grid_sur
        self.n_t = t_grid_sur.shape[0]
        # native uniform grid scalars (best-fit spacing minimizes index drift)
        self.t_g0 = float(t_grid_sur[0])
        self.t_dg = float((t_grid_sur[-1] - t_grid_sur[0]) / (self.n_t - 1))

        norm = np.load(os.path.join(self.sur_dir, "norm_factors.npz"))
        eim = np.load(os.path.join(self.sur_dir, "eim_B.npz"))
        self.eim_B = {}
        self.tck = {}
        for name in self.data_piece_names:
            # fold the scalar norm factor into the EIM matrix rows once
            eim_B = np.asarray(eim[f"eim_B_{name}"], dtype=np.float64)
            eim_B = eim_B * float(norm[f"norm_factor_{name}"])
            self.eim_B[name] = jnp.asarray(eim_B)
            tck = np.load(os.path.join(self.sur_dir, f"{name}_fits.npz"))
            self.tck[name] = (
                jnp.asarray(np.asarray(tck["t"], dtype=np.float64)),
                jnp.asarray(np.asarray(tck["c"], dtype=np.float64)),
            )
            tck.close()
        norm.close()
        eim.close()

    # -- time-range / bucketing helpers (host) --

    def _set_time_range(self, M, t_start, t_end):
        scale = M / self.sur_total_mass
        t_min = self.t_grid_sur[0] * scale
        t_max = self.t_grid_sur[-1] * scale
        if t_start is None:
            t_start = t_min
        if t_end is None:
            t_end = t_max
        if t_start < t_min or t_end > t_max:
            raise ValueError(
                f"Requested time range [{t_start}s, {t_end}s] is out of the "
                f"surrogate's range [{t_min}s, {t_max}s] for M={M:.2f}."
            )
        return t_start, t_end, scale

    @staticmethod
    def _start_trunc_index(grid, val):
        idx = int(np.searchsorted(grid, val, side="right")) - 1 - 5
        return idx if idx > 0 else 0

    # -- evaluation --

    def __call__(
        self,
        M,
        q,
        delta_t=None,
        t_start=None,
        t_end=None,
        times=None,
        remove_initial_phase=False,
        return_amp_phase_only=False,
    ):
        if delta_t is None and times is None:
            raise ValueError("Either delta_t or times must be provided.")

        t_start, t_end, scale = self._set_time_range(M, t_start, t_end)
        start_idx = self._start_trunc_index(self.t_grid_sur, t_start / scale)

        if times is None:
            num_samples = int((t_end - t_start) / delta_t) + 1
            new_t_grid = t_start + np.arange(num_samples) * delta_t
        else:
            new_t_grid = np.asarray(times, dtype=np.float64)
            num_samples = new_t_grid.shape[0]

        q = float(np.asarray(q).reshape(-1)[0])
        eta = float(eta_from_q(q))

        # Bucketed source slice: pull start_idx down to a power-of-two length.
        in_len_true = self.n_t - start_idx
        in_bucket = min(_bucket(in_len_true), self.n_t)
        in_start = self.n_t - in_bucket  # keep the tail, extend with real samples
        g0 = (self.t_g0 + in_start * self.t_dg) * scale
        dg = self.t_dg * scale
        amp_scale = scale / _amp_correction_factor

        if times is None:
            out_bucket = _bucket(num_samples)
            query = jnp.asarray(t_start + np.arange(out_bucket) * delta_t)
        else:
            out_bucket = _bucket(num_samples)
            padded = np.empty(out_bucket, dtype=np.float64)
            padded[:num_samples] = new_t_grid
            padded[num_samples:] = new_t_grid[-1]
            query = jnp.asarray(padded)

        amp, phase = _circ_amp_phase(
            self.tck["amp"][0],
            self.tck["amp"][1],
            self.tck["phase"][0],
            self.tck["phase"][1],
            self.eim_B["amp"],
            self.eim_B["phase"],
            jnp.asarray(eta),
            jnp.asarray(in_start),
            in_bucket,
            jnp.asarray(g0),
            jnp.asarray(dg),
            query,
            jnp.asarray(amp_scale),
            bool(remove_initial_phase),
        )

        amp = np.asarray(amp)[:num_samples]
        phase = np.asarray(phase)[:num_samples]
        if return_amp_phase_only:
            return amp, phase
        return new_t_grid, amp * np.exp(-1j * phase)


@partial(jax.jit, static_argnums=(8, 13))
def _circ_amp_phase(
    amp_t,
    amp_c,
    phase_t,
    phase_c,
    eim_amp,
    eim_phase,
    eta,
    in_start,
    in_bucket,
    g0,
    dg,
    query,
    amp_scale,
    remove_phase0,
):
    """Circular amp/phase core: tck node values -> bucketed EIM gemv -> uniform
    interpolation onto the (padded) output grid. Returns (amp*scale, phase).

    Static args: in_bucket (dynamic_slice size) and remove_phase0. All per-call
    scalars (eta, in_start, g0, dg, amp_scale) are traced so calls with
    different M/t_start/delta_t reuse the same compiled program per bucket."""
    amp_nodes = _cubic_bspline_vec_eval(amp_t, amp_c, eta)
    phase_nodes = _cubic_bspline_vec_eval(phase_t, phase_c, eta)

    in_start = in_start.astype(jnp.int64)
    zero = jnp.zeros((), jnp.int64)
    eim_amp_slice = jax.lax.dynamic_slice(
        eim_amp, (zero, in_start), (eim_amp.shape[0], in_bucket)
    )
    eim_phase_slice = jax.lax.dynamic_slice(
        eim_phase, (zero, in_start), (eim_phase.shape[0], in_bucket)
    )
    amp_native = amp_nodes @ eim_amp_slice
    phase_native = phase_nodes @ eim_phase_slice

    amp = _interp_uniform(amp_native, g0, dg, query) * amp_scale
    phase = _interp_uniform(phase_native, g0, dg, query)
    if remove_phase0:
        phase = phase - phase[0]
    return amp, phase
