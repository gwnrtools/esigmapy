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

    @staticmethod
    def _build_query(t_start, delta_t, num_samples, new_t_grid, times):
        """(padded query array of length out_bucket, out_bucket)."""
        out_bucket = _bucket(num_samples)
        if times is None:
            query = jnp.asarray(t_start + np.arange(out_bucket) * delta_t)
        else:
            padded = np.empty(out_bucket, dtype=np.float64)
            padded[:num_samples] = new_t_grid
            padded[num_samples:] = new_t_grid[-1]
            query = jnp.asarray(padded)
        return query, out_bucket

    def _trunc_setup(self, t_start_over_scale):
        """Host-side native-grid truncation: (in_start, in_bucket) covering the
        native samples from ``t_start_over_scale`` (native units) to the end."""
        start_idx = self._start_trunc_index(self.t_grid_sur, t_start_over_scale)
        in_bucket = min(_bucket(self.n_t - start_idx), self.n_t)
        in_start = self.n_t - in_bucket  # keep the tail, extend with real samples
        return in_start, in_bucket

    def _eval_padded(self, M, eta, scale, t_start, query, remove_initial_phase):
        """Return device amp/phase padded to len(query) (unsliced)."""
        in_start, in_bucket = self._trunc_setup(t_start / scale)
        return self._eval_padded_trunc(
            eta, scale, in_start, in_bucket, query, remove_initial_phase
        )

    def _eval_padded_trunc(
        self, eta, scale, in_start, in_bucket, query, remove_initial_phase
    ):
        """Like ``_eval_padded`` but with the host truncation precomputed, so
        ``scale`` (hence ``M``) may be a traced value."""
        g0 = (self.t_g0 + in_start * self.t_dg) * scale
        dg = self.t_dg * scale
        amp_scale = scale / _amp_correction_factor
        return _circ_amp_phase(
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

        if times is None:
            num_samples = int((t_end - t_start) / delta_t) + 1
            new_t_grid = t_start + np.arange(num_samples) * delta_t
        else:
            new_t_grid = np.asarray(times, dtype=np.float64)
            num_samples = new_t_grid.shape[0]

        q = float(np.asarray(q).reshape(-1)[0])
        eta = float(eta_from_q(q))

        query, out_bucket = self._build_query(
            t_start, delta_t, num_samples, new_t_grid, times
        )
        amp, phase = self._eval_padded(
            M, eta, scale, t_start, query, remove_initial_phase
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


# --- eccentric surrogate -----------------------------------------------------


def _tp3d_point(k0, k1, k2, coeffs, x):
    """Evaluate one 3D cubic tensor-product spline (nodes = knots[3:-3]) at the
    length-3 point ``x``. ``coeffs`` has shape (n0+2, n1+2, n2+2). Vectorized
    over a stack of fits via ``_tp3d_stacked``.

    The de Boor recurrence is unrolled by hand rather than delegated to
    ``TPI_jax._find_active_cubic_bspline_span_jax``: that kernel builds the
    basis with a dynamic-bound ``lax.fori_loop``, which reverse-mode
    differentiation (``jax.grad``) cannot trace through, and grad-ability over
    the waveform parameters is a design goal of this backend. (TPI_jax's own
    vector-interpolant evaluator, used for the shared-node pieces, unrolls the
    recurrence the same way and is grad-safe.)"""
    starts = []
    bases = []
    for axis, knots in enumerate((k0, k1, k2)):
        xa = x[axis]
        nodes = knots[3:-3]
        n_a = nodes.shape[0]
        span = jnp.clip(jnp.searchsorted(nodes, xa, side="right") + 2, 3, n_a + 1)
        k = jax.lax.dynamic_slice(knots, (span - 2,), (6,))
        l1 = xa - k[2]
        l2 = xa - k[1]
        l3 = xa - k[0]
        r1 = k[3] - xa
        r2 = k[4] - xa
        r3 = k[5] - xa

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
        bases.append(jnp.stack((d0, d1, d2, d3)))
        starts.append(span - 3)

    s0, s1, s2 = (s.astype(jnp.int64) for s in starts)
    block = jax.lax.dynamic_slice(coeffs, (s0, s1, s2), (4, 4, 4))
    return jnp.einsum("abc,a,b,c->", block, bases[0], bases[1], bases[2])


# Vectorized over a stack of fits: knots (M, n_a+6) per axis, coeffs (M, ...),
# points (M, 3). Each fit is evaluated at its own point.
_tp3d_stacked = jax.vmap(_tp3d_point, in_axes=(0, 0, 0, 0, 0))


def _construct_knots_np(nodes):
    nodes = np.asarray(nodes, dtype=np.float64)
    return np.concatenate((np.repeat(nodes[:1], 3), nodes, np.repeat(nodes[-1:], 3)))


class EccentricSurrogateJAX:
    """JAX port of ``EccentricSurrogate`` (waveform paths).

    Parameter-space evaluation uses two mechanisms, both set up once in
    ``__init__``:

    * the pieces sharing one grid and evaluated at the reference point
      [eta, e_ref, l_ref] (e, res_circ_phase, shifted_mean_anomaly,
      mean_anomaly_offset, x) are combined into ONE
      ``TP_Interpolant_ND_Vector`` via ``FromComponentSplines``;
    * res_amp and res_phase, whose fits have different nodes and their own
      evaluation points, are stacked and evaluated with a ``vmap`` over
      ``_tp3d_point``.
    """

    # Order of the data pieces stacked into the shared reference-point vector
    # interpolant. Per-piece component counts are NOT hard-coded: they are
    # inferred in _load from the number of fit files each piece provides (and
    # cross-checked against the piece's EIM basis count), so a change to the
    # underlying surrogate data does not require editing this class.
    _SHARED_PIECES = (
        "e",
        "res_circ_phase",
        "shifted_mean_anomaly",
        "mean_anomaly_offset",
        "x",
    )

    def __init__(self, ecc_data_dir, circ_data_dir):
        self.sur_dir = ecc_data_dir
        self.circ_sur = CircularSurrogateJAX(circ_data_dir)
        self._load()
        self.q_min, self.q_max = 1.0, 6.0
        self.e_ref_min, self.e_ref_max = 5.0e-7, 0.431
        self.l_ref_min, self.l_ref_max = 0.0, _TWO_PI
        self._core, self._core_orbital = self._make_core()

    # -- loading --

    def _load(self):
        filename = os.path.join(self.sur_dir, "surrogate_metadata.hdf")
        with h5py.File(filename, "r") as f:
            self.sur_total_mass = float(f["M"][()])
            t_grid_sur = np.asarray(f["t_grid_sur"][()], dtype=np.float64)
            l_grid_sur = np.asarray(f["l_grid_sur"][()], dtype=np.float64)
        self.t_grid_sur = t_grid_sur
        self.l_grid_sur = l_grid_sur
        self.n_t = t_grid_sur.shape[0]
        self.n_l = l_grid_sur.shape[0]
        self.t_g0 = float(t_grid_sur[0])
        self.t_dg = float((t_grid_sur[-1] - t_grid_sur[0]) / (self.n_t - 1))
        self.l_g0 = float(l_grid_sur[0])
        self.l_dg = float((l_grid_sur[-1] - l_grid_sur[0]) / (self.n_l - 1))

        norm = np.load(os.path.join(self.sur_dir, "norm_factors.npz"))
        eim = np.load(os.path.join(self.sur_dir, "eim_B.npz"))
        ei = np.load(os.path.join(self.sur_dir, "ei_indices.npz"))
        piece_names = [
            "res_amp",
            "res_phase",
            "res_circ_phase",
            "shifted_mean_anomaly",
            "e",
            "x",
        ]

        # EIM matrices with norm factor folded in.
        self.eim_B = {}
        for name in piece_names:
            m = np.asarray(eim[f"eim_B_{name}"], dtype=np.float64) * float(
                norm[f"norm_factor_{name}"]
            )
            self.eim_B[name] = jnp.asarray(m)
        self.ei_res_amp = np.asarray(ei["ei_indices_res_amp"])
        self.ei_res_phase = np.asarray(ei["ei_indices_res_phase"])
        norm.close()
        eim.close()
        ei.close()

        # e-EIM columns and l-grid values at the EIM nodes (precomputed).
        self.eim_e_ra = self.eim_B["e"][:, self.ei_res_amp]  # (n_e, n_ra)
        self.eim_e_rp = self.eim_B["e"][:, self.ei_res_phase]  # (n_e, n_rp)
        self.l_grid_ra = jnp.asarray(self.l_grid_sur[self.ei_res_amp])  # (n_ra,)
        self.l_grid_rp = jnp.asarray(self.l_grid_sur[self.ei_res_phase])  # (n_rp,)

        # Parameter-space fits: read each piece's per-fit (nodes, coeffs).
        raw_nodes = {}
        raw_coeffs = {}
        for name in piece_names:
            raw_nodes[name], raw_coeffs[name] = self._load_piece_raw(name)
        # res_amp / res_phase have different nodes per fit -> stacked vmap eval.
        self._res_amp_stack = self._build_stack(
            raw_nodes["res_amp"], raw_coeffs["res_amp"]
        )
        self._res_phase_stack = self._build_stack(
            raw_nodes["res_phase"], raw_coeffs["res_phase"]
        )

        # mean_anomaly_offset (single fit, shares the reference-point grid).
        mao_path = os.path.join(
            self.sur_dir, "fits/mean_anomaly_offset-ref_space-3D-fit_spline.h5"
        )
        with h5py.File(mao_path, "r") as f:
            shared_nodes_1d = [np.asarray(n, dtype=np.float64) for n in f["nodes"][()]]
            mao_coeffs = np.asarray(f["coefficients"][()], dtype=np.float64)

        # Shared-node vector interpolant via FromComponentSplines (one-time
        # setup): the _SHARED_PIECES, all on the common reference grid. Each
        # piece contributes as many components as it has fits on disk; for the
        # pieces with an EIM basis the two counts must agree.
        raw_coeffs["mean_anomaly_offset"] = [mao_coeffs]
        components = []
        self._shared_slices = {}
        offset = 0
        for name in self._SHARED_PIECES:
            count = len(raw_coeffs[name])
            if name in self.eim_B and self.eim_B[name].shape[0] != count:
                raise ValueError(
                    f"Data piece '{name}': {count} parameter-space fits on "
                    f"disk but the EIM basis has {self.eim_B[name].shape[0]} "
                    "rows; the surrogate data is inconsistent."
                )
            self._shared_slices[name] = slice(offset, offset + count)
            components.extend(raw_coeffs[name])
            offset += count
        self._vec_interp = TPI_jax.TP_Interpolant_ND_Vector.FromComponentSplines(
            shared_nodes_1d, components
        )

    def _load_piece_raw(self, name):
        """Return (nodes_per_fit, coeffs_list) for a data piece, fit-index order."""
        import re

        load_dir = os.path.join(self.sur_dir, f"fits/{name}_fits")
        pattern = re.compile(r"-(\d+)_spline\.h5$")
        indexed = {}
        for fn in os.listdir(load_dir):
            m = pattern.search(fn) if fn.endswith(".h5") else None
            if m:
                indexed[int(m.group(1))] = fn
        sorted_files = [indexed[i] for i in range(max(indexed) + 1)]
        nodes_per_fit = []
        coeffs_list = []
        for fn in sorted_files:
            with h5py.File(os.path.join(load_dir, fn), "r") as f:
                nodes_per_fit.append(
                    [np.asarray(n, dtype=np.float64) for n in f["nodes"][()]]
                )
                coeffs_list.append(np.asarray(f["coefficients"][()], dtype=np.float64))
        return nodes_per_fit, coeffs_list

    @staticmethod
    def _build_stack(nodes_per_fit, coeffs_list):
        knots = [
            jnp.asarray(np.stack([_construct_knots_np(nl[a]) for nl in nodes_per_fit]))
            for a in range(3)
        ]
        coeffs = jnp.asarray(np.stack(coeffs_list))
        return {"knots": knots, "coeffs": coeffs}

    # -- host helpers --

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

    def check_param_range(self, q, e_ref, l_ref):
        if not (self.q_min <= q <= self.q_max):
            raise ValueError(f"q={q} out of range [{self.q_min}, {self.q_max}].")
        if not (0.0 <= e_ref <= self.e_ref_max):
            raise ValueError(f"e_ref={e_ref} out of range [0, {self.e_ref_max}].")
        if not (self.l_ref_min <= l_ref <= self.l_ref_max):
            raise ValueError(
                f"l_ref={l_ref} out of range [{self.l_ref_min}, {self.l_ref_max}]."
            )

    # -- jitted core --

    def _make_core(self):
        vec = self._vec_interp
        ra_k0, ra_k1, ra_k2 = self._res_amp_stack["knots"]
        ra_c = self._res_amp_stack["coeffs"]
        rp_k0, rp_k1, rp_k2 = self._res_phase_stack["knots"]
        rp_c = self._res_phase_stack["coeffs"]
        eim_e_ra, eim_e_rp = self.eim_e_ra, self.eim_e_rp
        l_grid_ra, l_grid_rp = self.l_grid_ra, self.l_grid_rp
        eim_ra, eim_rp = self.eim_B["res_amp"], self.eim_B["res_phase"]
        eim_rcp = self.eim_B["res_circ_phase"]
        eim_sma = self.eim_B["shifted_mean_anomaly"]
        l_g0, l_dg = self.l_g0, self.l_dg
        eim_e_full, eim_x_full = self.eim_B["e"], self.eim_B["x"]
        s_amp, s_rcp, s_sma = (
            self._shared_slices["e"],
            self._shared_slices["res_circ_phase"],
            self._shared_slices["shifted_mean_anomaly"],
        )
        s_mao = self._shared_slices["mean_anomaly_offset"]
        s_x = self._shared_slices["x"]

        def body(
            eta,
            e_ref,
            l_ref,
            in_start_t,
            in_bucket_t,
            remove_phase0,
            tg0,
            tdg,
            query,
            amp0,
            phase0,
            amp_scale,
            with_orbital,
        ):
            shared = vec.TPInterpolationND(jnp.array([eta, e_ref, l_ref]))
            e_nodes = shared[s_amp]
            rcp_nodes = shared[s_rcp]
            sma_nodes = shared[s_sma]
            mao = shared[s_mao][0]

            e_eim_ra = e_nodes @ eim_e_ra  # (n_ra,)
            e_eim_rp = e_nodes @ eim_e_rp  # (n_rp,)
            l_eim_ra = jnp.mod(l_grid_ra + mao, _TWO_PI)
            l_eim_rp = jnp.mod(l_grid_rp + mao, _TWO_PI)
            eta_ra = jnp.full_like(e_eim_ra, eta)
            eta_rp = jnp.full_like(e_eim_rp, eta)
            ra_pts = jnp.stack((eta_ra, e_eim_ra, l_eim_ra), axis=1)  # (n_ra, 3)
            rp_pts = jnp.stack((eta_rp, e_eim_rp, l_eim_rp), axis=1)  # (n_rp, 3)
            ra_nodes = _tp3d_stacked(ra_k0, ra_k1, ra_k2, ra_c, ra_pts)
            rp_nodes = _tp3d_stacked(rp_k0, rp_k1, rp_k2, rp_c, rp_pts)

            zero = jnp.zeros((), jnp.int64)
            eim_sma_slice = jax.lax.dynamic_slice(
                eim_sma,
                (zero, in_start_t.astype(jnp.int64)),
                (eim_sma.shape[0], in_bucket_t),
            )
            lt_relation = sma_nodes @ eim_sma_slice  # (in_bucket_t,)
            delta_amp = ra_nodes @ eim_ra  # (n_l,)
            delta_phase = rp_nodes @ eim_rp
            res_circ = rcp_nodes @ eim_rcp

            l_s = _interp_uniform(lt_relation, tg0, tdg, query)
            dA = _interp_uniform(delta_amp, l_g0, l_dg, l_s)
            dP = _interp_uniform(delta_phase, l_g0, l_dg, l_s)
            rC = _interp_uniform(res_circ, l_g0, l_dg, l_s)

            amp = amp0 + amp_scale * dA
            phase = phase0 + rC + dP
            if remove_phase0:
                phase = phase - phase[0]
            if with_orbital:
                e_orb = _interp_uniform(e_nodes @ eim_e_full, l_g0, l_dg, l_s)
                x_orb = _interp_uniform(shared[s_x] @ eim_x_full, l_g0, l_dg, l_s)
                return amp, phase, e_orb, x_orb, l_s, mao
            return amp, phase

        # in_bucket_t (arg 4), remove_phase0 (5), with_orbital (6) are static.
        core = jax.jit(partial(body, with_orbital=False), static_argnums=(4, 5))
        core_orbital = jax.jit(partial(body, with_orbital=True), static_argnums=(4, 5))
        # Un-jitted body, kept so parameter_space_evaluator can trace it inside
        # user-supplied jit/vmap/grad transformations.
        self._body = body
        return core, core_orbital

    # -- evaluation --

    def _grid_setup(self, M, delta_t, t_start, t_end, times):
        """Host-side prologue shared by ``__call__`` and
        ``parameter_space_evaluator``: resolve the time range, build the output
        grid and padded device query, and compute the bucketed native-grid
        truncation scalars."""
        if delta_t is None and times is None:
            raise ValueError("Either delta_t or times must be provided.")

        t_start, t_end, scale = self._set_time_range(M, t_start, t_end)

        if times is None:
            num_samples = int((t_end - t_start) / delta_t) + 1
            new_t_grid = t_start + np.arange(num_samples) * delta_t
        else:
            new_t_grid = np.asarray(times, dtype=np.float64)
            num_samples = new_t_grid.shape[0]

        query, _ = self.circ_sur._build_query(
            t_start, delta_t, num_samples, new_t_grid, times
        )

        start_idx_t = CircularSurrogateJAX._start_trunc_index(
            self.t_grid_sur, t_start / scale
        )
        in_bucket_t = min(_bucket(self.n_t - start_idx_t), self.n_t)
        in_start_t = self.n_t - in_bucket_t
        tg0 = (self.t_g0 + in_start_t * self.t_dg) * scale
        tdg = self.t_dg * scale
        amp_scale = scale / _amp_correction_factor

        return (
            t_start,
            t_end,
            new_t_grid,
            num_samples,
            query,
            in_start_t,
            in_bucket_t,
            tg0,
            tdg,
            amp_scale,
        )

    def __call__(
        self,
        M,
        params,
        delta_t=None,
        t_start=None,
        t_end=None,
        times=None,
        remove_start_phase=True,
        return_orbital_variables=False,
    ):
        q, e_ref, l_ref = params
        self.check_param_range(q=q, e_ref=e_ref, l_ref=l_ref)

        (
            t_start,
            t_end,
            new_t_grid,
            num_samples,
            query,
            in_start_t,
            in_bucket_t,
            tg0,
            tdg,
            amp_scale,
        ) = self._grid_setup(M, delta_t, t_start, t_end, times)

        # Below the smallest supported eccentricity, fall back to the circular
        # surrogate (matching the numpy backend).
        if e_ref <= self.e_ref_min:
            if return_orbital_variables:
                raise NotImplementedError(
                    "return_orbital_variables with the circular fallback "
                    "(e_ref < e_ref_min) is not supported by the JAX backend."
                )
            if e_ref != 0.0:
                print(
                    f"Warning: e_ref={e_ref} < {self.e_ref_min}. Using circular "
                    "surrogate instead."
                )
            return self.circ_sur(
                M=M,
                q=q,
                delta_t=delta_t,
                t_start=t_start,
                t_end=t_end,
                times=(new_t_grid if times is not None else None),
                remove_initial_phase=True,
            )

        eta = float(eta_from_q(q))

        # Circular baseline amp/phase (initial phase removed), padded on device.
        circ_scale = M / self.circ_sur.sur_total_mass
        amp0, phase0 = self.circ_sur._eval_padded(
            M, eta, circ_scale, t_start, query, remove_initial_phase=True
        )

        core_args = (
            jnp.asarray(eta),
            jnp.asarray(float(e_ref)),
            jnp.asarray(float(l_ref)),
            jnp.asarray(in_start_t),
            in_bucket_t,
            bool(remove_start_phase),
            jnp.asarray(tg0),
            jnp.asarray(tdg),
            query,
            amp0,
            phase0,
            jnp.asarray(amp_scale),
        )

        if return_orbital_variables:
            amp, phase, e_orb, x_orb, l_s, mao = self._core_orbital(*core_args)
            amp = np.asarray(amp)[:num_samples]
            phase = np.asarray(phase)[:num_samples]
            e = np.asarray(e_orb)[:num_samples]
            x = np.asarray(x_orb)[:num_samples]
            l = np.asarray(l_s)[:num_samples] + float(mao)
            l -= 2 * np.pi * np.floor(l[0] / (2 * np.pi))
            orb_vars = {"e": e, "l": l, "x": x}
            return new_t_grid, orb_vars, amp * np.exp(-1j * phase)

        amp, phase = self._core(*core_args)
        amp = np.asarray(amp)[:num_samples]
        phase = np.asarray(phase)[:num_samples]
        return new_t_grid, amp * np.exp(-1j * phase)

    def parameter_space_evaluator(
        self,
        M=None,
        delta_t=None,
        t_start=None,
        t_end=None,
        times=None,
        remove_start_phase=True,
        output="mode",
        M_min=None,
    ):
        """Return ``(new_t_grid, fn)`` with ``fn`` a pure, JAX-traceable
        function of the waveform parameters.

        With ``output="mode"`` (the default) ``fn`` returns the complex (2,2)
        mode ``h = amp * exp(-1j * phase)`` on ``new_t_grid`` -- together with
        the returned time grid this is exactly the ``(t, h)`` output of
        ``__call__``. ``output="amp_phase"`` returns the ``(amp, phase)`` pair
        instead: two real arrays, the natural inputs for ``jax.grad`` of
        phase-dependent functionals. ``fn`` composes with ``jax.jit``,
        ``jax.vmap``, and ``jax.grad``/``jax.jacfwd``; for fast repeated calls
        wrap it in ``jax.jit`` (called eagerly it computes the same thing,
        just slower).

        Two calling conventions, selected by ``M``:

        * ``M`` given (a float): the full time-grid configuration (``M``,
          ``delta_t``/``times``, ``t_start``/``t_end``) is fixed host-side and
          ``fn(q, e_ref, l_ref)`` evaluates at one parameter-space point.
        * ``M=None``: only the output grid is fixed host-side (``t_start``,
          ``delta_t``, and optionally ``t_end`` <= 0, default 0), and
          ``fn(q, e_ref, l_ref, M, t_start=None)`` takes the total mass -- and
          optionally a per-waveform start time -- as *traced* arguments, so one
          ``jax.vmap`` batch may mix masses and durations on the common grid
          (the parameter-estimation workload: one data-segment grid, varying
          sources). Samples of ``new_t_grid`` earlier than the waveform's
          start, ``max(t_start, t_min(M))`` where ``t_min(M)`` is the
          surrogate's earliest time for that mass, are exactly zero (mode) or
          zeroed amp/phase; with ``remove_start_phase=True`` the phase is
          referenced at that start. Pass ``M_min``, a lower bound on every
          ``M`` the returned ``fn`` will see, to truncate the native grids
          host-side (worst case is the smallest mass); without it the full
          native grids are used -- always correct, but slower on CPU.
          A user-supplied ``times`` array is not supported in this mode.

        Unlike ``__call__``, ``fn`` performs NO parameter-range checks (traced
        values cannot be inspected) and has no circular fallback at very small
        ``e_ref``: out-of-range parameters are silently extrapolated with the
        boundary spline segments, exactly as in TPI's jit-traced evaluation.
        Validate inputs with ``check_param_range`` beforehand if that matters.
        """
        if output not in ("mode", "amp_phase"):
            raise ValueError(f"output must be 'mode' or 'amp_phase', got {output!r}.")
        if M is None:
            return self._traced_mass_evaluator(
                delta_t, t_start, t_end, times, remove_start_phase, output, M_min
            )
        if M_min is not None:
            raise ValueError("M_min only applies to the traced-mass mode (M=None).")

        (
            t_start,
            _t_end,
            new_t_grid,
            num_samples,
            query,
            in_start_t,
            in_bucket_t,
            tg0,
            tdg,
            amp_scale,
        ) = self._grid_setup(M, delta_t, t_start, t_end, times)

        body = self._body
        circ = self.circ_sur
        circ_scale = M / circ.sur_total_mass
        in_start_j = jnp.asarray(in_start_t)
        tg0_j = jnp.asarray(tg0)
        tdg_j = jnp.asarray(tdg)
        amp_scale_j = jnp.asarray(amp_scale)
        remove = bool(remove_start_phase)

        def fn(q, e_ref, l_ref):
            eta = q / (1.0 + q) ** 2  # eta_from_q, in traceable form
            amp0, phase0 = circ._eval_padded(
                M, eta, circ_scale, t_start, query, remove_initial_phase=True
            )
            amp, phase = body(
                eta,
                e_ref,
                l_ref,
                in_start_j,
                in_bucket_t,
                remove,
                tg0_j,
                tdg_j,
                query,
                amp0,
                phase0,
                amp_scale_j,
                with_orbital=False,
            )
            amp = amp[:num_samples]
            phase = phase[:num_samples]
            if output == "mode":
                return amp * jnp.exp(-1j * phase)
            return amp, phase

        return new_t_grid, fn

    def _traced_mass_evaluator(
        self, delta_t, t_start, t_end, times, remove_start_phase, output, M_min
    ):
        """Build the ``M=None`` evaluator: fixed output grid, traced mass and
        per-waveform start time. See ``parameter_space_evaluator``."""
        if times is not None:
            raise NotImplementedError(
                "A user-supplied `times` grid is not supported with M=None; "
                "provide t_start and delta_t (uniform common grid)."
            )
        if delta_t is None or t_start is None:
            raise ValueError("With M=None, both t_start and delta_t are required.")
        if t_end is None:
            t_end = 0.0
        if t_end > 0.0 or t_start >= t_end:
            raise ValueError(
                f"Need t_start < t_end <= 0 (t=0 is the end of the inspiral); "
                f"got t_start={t_start}, t_end={t_end}."
            )

        num_samples = int((t_end - t_start) / delta_t) + 1
        new_t_grid = t_start + np.arange(num_samples) * delta_t
        query, _ = self.circ_sur._build_query(
            t_start, delta_t, num_samples, new_t_grid, None
        )

        # Host-side native-grid truncation. The needed native span is widest
        # for the smallest mass (t_start/scale most negative), so a bound
        # computed at M_min is safe for every M >= M_min; without a bound,
        # keep the full grids.
        if M_min is not None:
            start_idx_t = CircularSurrogateJAX._start_trunc_index(
                self.t_grid_sur, t_start / (M_min / self.sur_total_mass)
            )
            in_bucket_t = min(_bucket(self.n_t - start_idx_t), self.n_t)
            in_start_t = self.n_t - in_bucket_t
            in_start_c, in_bucket_c = self.circ_sur._trunc_setup(
                t_start / (M_min / self.circ_sur.sur_total_mass)
            )
        else:
            in_start_t, in_bucket_t = 0, self.n_t
            in_start_c, in_bucket_c = 0, self.circ_sur.n_t

        body = self._body
        circ = self.circ_sur
        sur_mass = self.sur_total_mass
        circ_mass = circ.sur_total_mass
        t_g0, t_dg = self.t_g0, self.t_dg
        in_start_j = jnp.asarray(in_start_t)
        t_grid_j = jnp.asarray(new_t_grid)
        seg_start = float(new_t_grid[0])
        remove = bool(remove_start_phase)

        def fn(q, e_ref, l_ref, M, t_start=None):
            eta = q / (1.0 + q) ** 2  # eta_from_q, in traceable form
            scale = M / sur_mass
            amp0, phase0 = circ._eval_padded_trunc(
                eta,
                M / circ_mass,
                in_start_c,
                in_bucket_c,
                query,
                remove_initial_phase=True,
            )
            amp, phase = body(
                eta,
                e_ref,
                l_ref,
                in_start_j,
                in_bucket_t,
                False,  # phase is re-referenced at the waveform start below
                (t_g0 + in_start_t * t_dg) * scale,
                t_dg * scale,
                query,
                amp0,
                phase0,
                scale / _amp_correction_factor,
                with_orbital=False,
            )
            amp = amp[:num_samples]
            phase = phase[:num_samples]
            # Effective waveform start: the caller's t_start, clipped to the
            # surrogate's earliest time for this mass.
            t_min = t_g0 * scale
            start_eff = t_min if t_start is None else jnp.maximum(t_start, t_min)
            if remove:
                phase = phase - _interp_uniform(phase, seg_start, delta_t, start_eff)
            keep = t_grid_j >= start_eff
            if output == "mode":
                return jnp.where(keep, amp * jnp.exp(-1j * phase), 0.0)
            return jnp.where(keep, amp, 0.0), jnp.where(keep, phase, 0.0)

        return new_t_grid, fn
