# Copyright (C) 2026 Akash Maurya

"""ESIGMASur surrogate base evaluation code"""

import os

os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from numba import njit
import h5py
import scipy.interpolate as si
import TPI
from lal import SpinWeightedSphericalHarmonic
from pycbc.conversions import eta_from_q
import re

_amp_correction_factor = np.real(
    SpinWeightedSphericalHarmonic(0, 0, -2, 2, 2)
)  # The amplitude normalization factor saved in the surrogate files is off by this amount. So we correct for it during surrogate evaluation by using this factor.


@njit(fastmath=True)
def mode_from_amp_phase(amp, phase):
    return amp * np.exp(-1j * phase)


@njit(fastmath=True, cache=True)
def _fused_interp_mode_uniform(
    g0, dg, amp_table, phase_table, q0, dq, n_out, amp_scale, remove_phase0
):
    """Single-pass amp/phase interpolation onto a uniform output grid, followed
    by amplitude scaling, optional initial-phase removal, and complex-mode
    assembly ``amp * exp(-1j*phase)``.

    Both the native (source) grid and the output/query grid are uniform, so the
    source index is found by arithmetic instead of a binary search. ``np.interp``
    boundary clamping is reproduced exactly: queries at or beyond the grid ends
    take the boundary table value. Fusing all steps writes the large output
    array once and reads the source tables once, avoiding the several full-size
    temporaries of the ``np.interp`` -> ``*=`` -> ``exp`` pipeline.
    """
    ntab = amp_table.shape[0]
    last = ntab - 1
    inv_dg = 1.0 / dg
    out = np.empty(n_out, dtype=np.complex128)

    if remove_phase0:
        pos0 = (q0 - g0) * inv_dg
        if pos0 <= 0.0:
            phase0 = phase_table[0]
        elif pos0 >= last:
            phase0 = phase_table[last]
        else:
            i0 = int(pos0)
            w0 = pos0 - i0
            phase0 = phase_table[i0] + w0 * (phase_table[i0 + 1] - phase_table[i0])
    else:
        phase0 = 0.0

    for i in range(n_out):
        pos = (q0 + i * dq - g0) * inv_dg
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
def _fused_interp_amp_phase_uniform(
    g0, dg, amp_table, phase_table, q0, dq, n_out, amp_scale, remove_phase0
):
    """As ``_fused_interp_mode_uniform`` but returns the scaled ``(amp, phase)``
    arrays instead of the complex mode (the ``return_amp_phase_only`` path)."""
    ntab = amp_table.shape[0]
    last = ntab - 1
    inv_dg = 1.0 / dg
    amp = np.empty(n_out, dtype=np.float64)
    phase = np.empty(n_out, dtype=np.float64)

    if remove_phase0:
        pos0 = (q0 - g0) * inv_dg
        if pos0 <= 0.0:
            phase0 = phase_table[0]
        elif pos0 >= last:
            phase0 = phase_table[last]
        else:
            i0 = int(pos0)
            w0 = pos0 - i0
            phase0 = phase_table[i0] + w0 * (phase_table[i0 + 1] - phase_table[i0])
    else:
        phase0 = 0.0

    for i in range(n_out):
        pos = (q0 + i * dq - g0) * inv_dg
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


def _uniform_grid_scalars(grid, start_idx):
    """First value and best-fit uniform spacing of ``grid[start_idx:]``.

    The spacing is taken as the average over the retained span (endpoints /
    count) rather than a single local difference, which minimizes drift of the
    arithmetic index ``(query - g0) / dg`` accumulated over long output grids.
    """
    last = grid.shape[0] - 1
    g0 = grid[start_idx]
    dg = (grid[last] - g0) / (last - start_idx)
    return g0, dg


@njit(fastmath=True, cache=True, inline="always")
def _lerp_uniform(x, g0, inv_dg, table, last):
    """Linear interpolation of ``table`` (on a uniform grid ``g0`` with spacing
    ``1/inv_dg``) at query ``x``, with np.interp boundary clamping."""
    pos = (x - g0) * inv_dg
    if pos <= 0.0:
        return table[0]
    if pos >= last:
        return table[last]
    idx = int(pos)
    w = pos - idx
    return table[idx] + w * (table[idx + 1] - table[idx])


@njit(fastmath=True, cache=True)
def _fused_ecc_mode_uniform(
    tg0,
    tdg,
    lt_relation,
    q0,
    dq,
    n_out,
    lg0,
    ldg,
    damp,
    dphase,
    rcp,
    amp0,
    phase0,
    amp_scale,
    remove_phase0,
):
    """Single-pass eccentric final stage.

    For each uniform output sample the mean-anomaly ``l_s`` is interpolated from
    ``lt_relation`` on the uniform (scaled) time grid, then the three residual
    data pieces (delta amplitude, delta phase, residual circular phase) are
    interpolated from that ``l_s`` on the uniform mean-anomaly grid -- the
    index/weight into the l-grid is computed once and shared by all three. The
    circular baseline (``amp0``/``phase0``) is combined in, the initial phase is
    optionally removed, and the complex mode is assembled, all in one pass over
    the output grid instead of four np.interp calls plus several full-size
    temporaries.
    """
    last_t = lt_relation.shape[0] - 1
    inv_tdg = 1.0 / tdg
    last_l = damp.shape[0] - 1
    inv_ldg = 1.0 / ldg
    out = np.empty(n_out, dtype=np.complex128)

    if remove_phase0:
        ls0 = _lerp_uniform(q0, tg0, inv_tdg, lt_relation, last_t)
        phase0_ref = (
            phase0[0]
            + _lerp_uniform(ls0, lg0, inv_ldg, rcp, last_l)
            + _lerp_uniform(ls0, lg0, inv_ldg, dphase, last_l)
        )
    else:
        phase0_ref = 0.0

    for i in range(n_out):
        ls = _lerp_uniform(q0 + i * dq, tg0, inv_tdg, lt_relation, last_t)
        # Shared l-grid position for the three residual pieces.
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
        ph = phase0[i] + rcpi + dP - phase0_ref
        out[i] = amp * (np.cos(ph) - 1j * np.sin(ph))
    return out


def _unwrap_single_float(val):
    if isinstance(val, (float, int, np.floating, np.integer)):
        return float(val)
    elif isinstance(val, np.ndarray) and val.size == 1:
        return float(val.item())
    else:
        raise ValueError(
            f"Expected a float or a numpy array of size 1, got: {val} ({type(val)})"
        )


class Surrogate:
    def __init__(self, data_dir):
        # The directory where the surrogate data is stored
        self.sur_dir = data_dir
        # Surrogate data pieces. This will be set in the child classes, as the circular and eccentric surrogates are built of different data pieces.
        self.data_piece_names = None
        # Normalization factors for the surrogate data pieces, which are used to un-normalize the surrogate data pieces before reconstructing the waveform; read from the surrogate files via self.load_norm_factors().
        self.norm_factor = {}
        # The EIM B-matrices for the surrogate data pieces; read from the surrogate files via self.load_eim_B_matrices().
        self.eim_B = {}
        # The parametric fits for the surrogate data pieces at EI nodes; read from the surrogate files via self.load_param_space_fits(), defined in child classes.
        self.fit = {}

    def get_metadata(self, key):
        """
        Returns metadata of the surrogate

        Parameters:
        -----------
        key -- Metadata key (string) or list of keys to return.

        Returns:
        --------
        Metadata value(s) corresponding to the key(s).
        If a single key is provided, returns the corresponding value,
        otherwise returns a list of values corresponding to the list of keys.
        """
        filename = os.path.join(self.sur_dir, "surrogate_metadata.hdf")
        with h5py.File(filename, "r") as f:
            if isinstance(key, str):
                return f[key][()]
            return [f[k][()] for k in key]

    def load_norm_factors(self):
        filename_norm_factors = os.path.join(self.sur_dir, f"norm_factors.npz")
        norm_factor_dataset = np.load(filename_norm_factors)
        for data_piece_name in self.data_piece_names:
            self.norm_factor[data_piece_name] = norm_factor_dataset[
                f"norm_factor_{data_piece_name}"
            ]
        norm_factor_dataset.close()

    def load_eim_B_matrices(self, filename=None, data_piece_names=None):
        if filename is None:
            filename = "eim_B.npz"
        if data_piece_names is None:
            data_piece_names = self.data_piece_names
        if isinstance(data_piece_names, str):
            data_piece_names = [data_piece_names]
        filename_eim = os.path.join(self.sur_dir, filename)
        eim_B_dataset = np.load(filename_eim)
        for data_piece_name in data_piece_names:
            self.eim_B[data_piece_name] = eim_B_dataset[f"eim_B_{data_piece_name}"]
        eim_B_dataset.close()

    def _set_time_range(self, M, times, t_start, t_end):
        mass_scaling_factor = M / self.sur_total_mass

        t_min_sur = self.t_grid_sur[0] * mass_scaling_factor
        t_max_sur = self.t_grid_sur[-1] * mass_scaling_factor

        if times is None:
            if t_start is None:
                t_start = t_min_sur
            if t_end is None:
                t_end = t_max_sur
        else:
            t_start = times[0]
            t_end = times[-1]

        if t_start < t_min_sur or t_end > t_max_sur:
            raise ValueError(
                f"""Requested time range [{t_start}s, {t_end}s] is out of the surrogate's time range of [{t_min_sur}s, {t_max_sur}s] for the given total mass of {M:.2f} M_sun. 
Please choose a time interval within these bounds."""
            )

        return t_start, t_end, mass_scaling_factor

    @staticmethod
    def _find_conservative_starting_truncation_index(grid, val):
        idx = np.searchsorted(grid, val, side="right") - 1
        idx -= 5  # Leaving some buffer to avoid edge effects in spline interpolation
        # Checking if we even have data this deep in the starting
        if idx < 0:
            idx = 0
        return idx


class CircularSurrogate(Surrogate):
    def __init__(self, data_dir):
        super().__init__(data_dir=data_dir)
        self.data_piece_names = [
            "amp",
            "phase",
            "x",
            "l",
        ]  # The surrogate data pieces for the circular surrogate. The "x" and "l" data pieces are used only for reconstructing the orbital variables, and are not needed if one only wants the waveform.
        self.load_sur_metadata()
        self.load_norm_factors()
        self.load_eim_B_matrices(
            filename="eim_B.npz", data_piece_names=["amp", "phase"]
        )
        self.load_eim_B_matrices(
            filename="eim_B-orb_vars.npz", data_piece_names=["x", "l"]
        )
        self.load_param_space_fits()

    def load_sur_metadata(self):
        self.sur_total_mass, self.t_grid_sur = self.get_metadata(["M", "t_grid_sur"])

    @staticmethod
    def load_interpolant(filename):
        # scipy BSplines
        tck = np.load(filename)
        vec_interpolant = si.BSpline(tck["t"], tck["c"], tck["k"], extrapolate=False)
        return vec_interpolant

    def load_param_space_fits(self):
        for data_piece_name in self.data_piece_names:
            filename_interp = os.path.join(self.sur_dir, f"{data_piece_name}_fits.npz")
            self.fit[data_piece_name] = self.load_interpolant(filename_interp)

    def __call__(
        self,
        M,
        q,
        reference_mean_anomaly=0.0,  # This is only used for returning orbital variables
        delta_t=None,
        t_start=None,
        t_end=None,
        times=None,
        remove_initial_phase=False,
        return_amp_phase_only=False,
        return_orbital_variables=False,
    ):
        if delta_t is None and times is None:
            raise ValueError("Either delta_t or times must be provided.")

        t_grid_sur = self.t_grid_sur

        t_start, t_end, mass_scaling_factor = self._set_time_range(
            M=M, times=times, t_start=t_start, t_end=t_end
        )

        start_idx = self._find_conservative_starting_truncation_index(
            grid=t_grid_sur, val=t_start / mass_scaling_factor
        )

        if times is None:
            num_samples = int((t_end - t_start) / delta_t) + 1
            # The output grid array itself is built lazily: the amp/phase-only
            # fast path never returns it, and for long waveforms materializing
            # the multi-megasample array is a measurable cost.
            new_t_grid = None
        else:
            new_t_grid = times

        q = _unwrap_single_float(
            q
        )  # This is to ensure that q is a single value and not an array
        eta = eta_from_q(q)

        amp_node_vals = self.fit["amp"](eta)
        phase_node_vals = self.fit["phase"](eta)

        amp_native = self.norm_factor["amp"] * np.dot(
            amp_node_vals, self.eim_B["amp"][:, start_idx:]
        )
        phase_native = self.norm_factor["phase"] * np.dot(
            phase_node_vals, self.eim_B["phase"][:, start_idx:]
        )

        amp_scale = (
            mass_scaling_factor / _amp_correction_factor
        )  # Correcting for the amplitude normalization factor that was off in the surrogate files.

        # Fast path: for a uniform output grid (times is None) the amp/phase
        # interpolation, amplitude scaling, optional initial-phase removal, and
        # complex-mode assembly are fused into a single pass over the output
        # grid. Both the native and output grids are uniform, so interpolation
        # uses index arithmetic. When the caller supplies arbitrary `times`, the
        # generic np.interp path below is used instead.
        if times is None and not return_orbital_variables:
            g0, dg = _uniform_grid_scalars(t_grid_sur, start_idx)
            g0 *= mass_scaling_factor
            dg *= mass_scaling_factor
            if return_amp_phase_only:
                return _fused_interp_amp_phase_uniform(
                    g0,
                    dg,
                    amp_native,
                    phase_native,
                    t_start,
                    delta_t,
                    num_samples,
                    amp_scale,
                    remove_initial_phase,
                )
            new_t_grid = t_start + np.arange(num_samples) * delta_t
            return new_t_grid, _fused_interp_mode_uniform(
                g0,
                dg,
                amp_native,
                phase_native,
                t_start,
                delta_t,
                num_samples,
                amp_scale,
                remove_initial_phase,
            )

        if new_t_grid is None:
            new_t_grid = t_start + np.arange(num_samples) * delta_t
        t_grid_sur = t_grid_sur[start_idx:] * mass_scaling_factor
        amp = np.interp(new_t_grid, t_grid_sur, amp_native)
        phase = np.interp(new_t_grid, t_grid_sur, phase_native)

        amp *= amp_scale

        if remove_initial_phase:
            phase -= phase[0]
        if return_amp_phase_only:
            return amp, phase

        if return_orbital_variables:
            e = np.zeros(len(new_t_grid))

            l_node_vals = self.fit["l"](eta)
            l_native = self.norm_factor["l"] * np.dot(
                l_node_vals, self.eim_B["l"][:, start_idx:]
            )
            l = (
                np.interp(new_t_grid, t_grid_sur, l_native) + reference_mean_anomaly
            )  # The circular surrogate mean anomaly data was prepared such that the mean anomaly at the reference time is 0, so we add the reference_mean_anomaly here to get the correct mean anomaly values at the reference time. This is physically correct because circular waveforms just differing in mean anomalies are simply constant phase-shifted versions of each other
            # Caution: This also means that one should make sure during surrogate construction that the reference time of the circular surrogate is the same as that of the eccentric surrogate.
            l -= (
                2 * np.pi * np.floor(l[0] / (2 * np.pi))
            )  # Bringing starting value of mean anomaly in [0, 2pi)

            x_node_vals = self.fit["x"](eta)
            x_native = self.norm_factor["x"] * np.dot(
                x_node_vals, self.eim_B["x"][:, start_idx:]
            )
            x = np.interp(new_t_grid, t_grid_sur, x_native)

            orb_vars = {"e": e, "l": l, "x": x}
            return new_t_grid, orb_vars, mode_from_amp_phase(amp, phase)

        return new_t_grid, mode_from_amp_phase(amp, phase)


class EccentricSurrogate(Surrogate):
    def __init__(self, ecc_data_dir, circ_data_dir, verbose=False):
        super().__init__(data_dir=ecc_data_dir)
        self.circ_sur_dir = circ_data_dir
        self.circ_sur = CircularSurrogate(data_dir=self.circ_sur_dir)

        self.data_piece_names = [
            "res_amp",
            "res_phase",
            "res_circ_phase",
            "shifted_mean_anomaly",
            "e",
            "x",
        ]
        self.ei_indices = {}

        self.load_sur_metadata()
        self.load_norm_factors()
        self.load_eim_B_matrices()
        self.load_ei_indices()
        self.load_param_space_fits(verbose=verbose)

        self.q_min = 1.0
        self.q_max = 6.0
        self.e_ref_min = 5.0e-7  # The TPI fits below this value fail to evaluate
        self.e_ref_max = 0.431
        self.l_ref_min = 0.0
        self.l_ref_max = 2 * np.pi

    def check_param_range(self, q, e_ref, l_ref, override=False):
        if not override:
            if not (self.q_min <= q <= self.q_max):
                raise ValueError(
                    f"Mass ratio q={q} is out of range [{self.q_min}, {self.q_max}]. Please choose a value within this range."
                )
            if not (0.0 <= e_ref <= self.e_ref_max):
                raise ValueError(
                    f"Reference eccentricity e_ref={e_ref} is out of range [{0.}, {self.e_ref_max}]. Please choose a value within this range."
                )
            if not (self.l_ref_min <= l_ref <= self.l_ref_max):
                raise ValueError(
                    f"Reference mean anomaly l_ref={l_ref} is out of range [{self.l_ref_min}, {self.l_ref_max}]. Please choose a value within this range."
                )

    def load_sur_metadata(self):
        (
            self.sur_total_mass,
            self.t_ref,
            self.t_grid_sur,
            self.l_grid_sur,
        ) = self.get_metadata(["M", "t_ref", "t_grid_sur", "l_grid_sur"])

    def load_ei_indices(self):
        filename_ei_indices = os.path.join(self.sur_dir, f"ei_indices.npz")
        ei_indices_dataset = np.load(filename_ei_indices)
        for data_piece_name in ["res_amp", "res_phase"]:
            self.ei_indices[data_piece_name] = ei_indices_dataset[
                f"ei_indices_{data_piece_name}"
            ]
        ei_indices_dataset.close()

    @staticmethod
    def _get_sorted_fit_filenames(filenames):
        # Use a regex to extract the number before .h5, regardless of the prefix
        pattern = re.compile(r"-(\d+)_spline\.h5$")
        indexed_files = {}
        for fname in filenames:
            match = pattern.search(
                fname
            )  # search instead of match allows for variable prefix
            if match:
                idx = int(match.group(1))
                indexed_files[idx] = fname

        # Create array of appropriate size
        max_index = max(indexed_files.keys())
        sorted_filenames = [None] * (max_index + 1)
        # Fill the array
        for idx, fname in indexed_files.items():
            sorted_filenames[idx] = fname
        return sorted_filenames

    @staticmethod
    def _read_fit_raw(filepath):
        with h5py.File(filepath, "r") as f:
            nodes = [np.asarray(n, dtype=np.float64) for n in f["nodes"][()]]
            coeffs = np.asarray(f["coefficients"][()], dtype=np.float64)
        return nodes, coeffs

    @staticmethod
    def load_fit(filepath):
        nodes, coeffs = EccentricSurrogate._read_fit_raw(filepath)
        return TPI.TP_Interpolant_ND(nodes, coeffs=coeffs)

    def load_param_space_fits(self, verbose=False):
        raw_coeffs = {}

        for data_piece_name in self.data_piece_names:
            load_dir = os.path.join(self.sur_dir, f"fits/{data_piece_name}_fits")
            # List of all the files which end in .h5
            filenames = [f for f in os.listdir(load_dir) if f.endswith(".h5")]
            sorted_filenames = self._get_sorted_fit_filenames(filenames=filenames)
            self.fit[data_piece_name] = []
            coeffs_list = []

            for filename in sorted_filenames:
                filepath = os.path.join(load_dir, filename)
                if verbose:
                    print(f"Loading spline interpolant of GPR from {filename}...")
                nodes, coeffs = self._read_fit_raw(filepath)
                self.fit[data_piece_name].append(
                    TPI.TP_Interpolant_ND(nodes, coeffs=coeffs)
                )
                coeffs_list.append(coeffs)

            raw_coeffs[data_piece_name] = coeffs_list

        filepath = os.path.join(
            self.sur_dir, f"fits/mean_anomaly_offset-ref_space-3D-fit_spline.h5"
        )
        mao_nodes, mao_coeffs = self._read_fit_raw(filepath)
        self.fit["mean_anomaly_offset_fit"] = [
            TPI.TP_Interpolant_ND(mao_nodes, coeffs=mao_coeffs)
        ]
        raw_coeffs["mean_anomaly_offset"] = [mao_coeffs]

        # Every piece evaluated at the shared reference point [eta, e_ref, l_ref]
        # lives on the same parameter-space grid (the mean_anomaly_offset grid),
        # so all their fits are combined into ONE vector-valued TPI interpolant:
        # one evaluation returns every component, sharing the span search and
        # basis computation across all of them.
        shared_piece_names = [
            "e",
            "res_circ_phase",
            "shifted_mean_anomaly",
            "mean_anomaly_offset",
            "x",
        ]
        self._shared_ref_slices = {}
        components = []
        offset = 0
        for name in shared_piece_names:
            count = len(raw_coeffs[name])
            self._shared_ref_slices[name] = slice(offset, offset + count)
            components.extend(raw_coeffs[name])
            offset += count
        self._shared_ref_interp = TPI.TP_Interpolant_ND_Vector.FromComponentSplines(
            mao_nodes, components
        )

    def _eval_shared_ref(self, eta, e_ref, l_ref):
        """Evaluate all shared-node reference-point pieces at once; returns a dict
        of per-piece node-value arrays keyed by data-piece name."""
        vals = self._shared_ref_interp.TPInterpolationND(np.array([eta, e_ref, l_ref]))
        return {name: vals[sl] for name, sl in self._shared_ref_slices.items()}

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
        if delta_t is None and times is None:
            raise ValueError("Either delta_t or times must be provided.")

        t_grid_sur = self.t_grid_sur  # The native time grid of the surrogate
        l_grid_sur = (
            self.l_grid_sur
        )  # The native shifted mean anomaly grid of the surrogate

        t_start, t_end, mass_scaling_factor = self._set_time_range(
            M=M, times=times, t_start=t_start, t_end=t_end
        )

        start_idx_t = self._find_conservative_starting_truncation_index(
            grid=t_grid_sur, val=t_start / mass_scaling_factor
        )
        # t_grid_sur stays the full native grid here; the scaled/truncated array
        # is only materialized on the generic (non-fused) path below.

        if times is None:
            num_samples = int((t_end - t_start) / delta_t) + 1
            new_t_grid = t_start + np.arange(num_samples) * delta_t
        else:
            new_t_grid = times

        q, e_ref, l_ref = params
        self.check_param_range(q=q, e_ref=e_ref, l_ref=l_ref)

        if e_ref > self.e_ref_min:
            # Pass times=None (with delta_t/t_start/t_end) when this call built
            # the output grid itself, so the circular surrogate rebuilds the
            # bitwise-identical uniform grid and takes its fused fast path. Only
            # forward an explicit `times` array when the caller supplied one.
            amp0_, phase0_ = self.circ_sur(
                M=M,
                q=q,
                delta_t=delta_t,
                t_start=t_start,
                t_end=t_end,
                times=(None if times is None else new_t_grid),
                remove_initial_phase=True,
                return_amp_phase_only=True,
                return_orbital_variables=False,
            )
        else:
            if e_ref != 0.0:
                print(
                    f"Warning: e_ref={e_ref} < {self.e_ref_min}. Setting e_ref to 0 and using circular surrogate instead."
                )
            return self.circ_sur(
                M=M,
                q=q,
                reference_mean_anomaly=l_ref,
                delta_t=delta_t,
                t_start=t_start,
                t_end=t_end,
                times=new_t_grid,
                remove_initial_phase=True,
                return_orbital_variables=return_orbital_variables,
            )
        eta = eta_from_q(q)

        # All parametric fits evaluated at the shared reference point
        # [eta, e_ref, l_ref] (e, res_circ_phase, shifted_mean_anomaly,
        # mean_anomaly_offset, x) are computed in one shared-basis contraction.
        shared_ref_vals = self._eval_shared_ref(eta, e_ref, l_ref)
        e_node_vals = shared_ref_vals["e"]
        res_circ_phase_node_vals = shared_ref_vals["res_circ_phase"]
        shifted_mean_anomaly_node_vals = shared_ref_vals["shifted_mean_anomaly"]

        e_eim_res_amp = self.norm_factor["e"] * np.dot(
            e_node_vals, self.eim_B["e"][:, self.ei_indices["res_amp"]]
        )
        e_eim_res_phase = self.norm_factor["e"] * np.dot(
            e_node_vals, self.eim_B["e"][:, self.ei_indices["res_phase"]]
        )

        mean_anomaly_offset_of_shifted_mean_anomaly = float(
            shared_ref_vals["mean_anomaly_offset"][0]
        )  # This is the mean anomaly offset of the shifted mean anomaly
        # Caution: It's important here that l_grid_sur is the un-truncated, full sur.l_grid_sur
        l_eim_res_amp = (
            l_grid_sur[self.ei_indices["res_amp"]]
            + mean_anomaly_offset_of_shifted_mean_anomaly
        ) % (2 * np.pi)
        l_eim_res_phase = (
            l_grid_sur[self.ei_indices["res_phase"]]
            + mean_anomaly_offset_of_shifted_mean_anomaly
        ) % (2 * np.pi)

        # Each res_amp/res_phase fit has its own parameter-space nodes and its
        # own evaluation point, so they cannot share one vector-valued
        # interpolant; evaluate the per-fit TPI splines one by one.
        res_amp_fits = self.fit["res_amp"]
        res_amp_points = np.empty((e_eim_res_amp.shape[0], 3))
        res_amp_points[:, 0] = eta
        res_amp_points[:, 1] = e_eim_res_amp
        res_amp_points[:, 2] = l_eim_res_amp
        res_amp_node_vals = np.fromiter(
            (
                fit.TPInterpolationND(point)
                for fit, point in zip(res_amp_fits, res_amp_points)
            ),
            dtype=np.float64,
            count=len(res_amp_fits),
        )

        res_phase_fits = self.fit["res_phase"]
        res_phase_points = np.empty((e_eim_res_phase.shape[0], 3))
        res_phase_points[:, 0] = eta
        res_phase_points[:, 1] = e_eim_res_phase
        res_phase_points[:, 2] = l_eim_res_phase
        res_phase_node_vals = np.fromiter(
            (
                fit.TPInterpolationND(point)
                for fit, point in zip(res_phase_fits, res_phase_points)
            ),
            dtype=np.float64,
            count=len(res_phase_fits),
        )

        lt_relation = self.norm_factor["shifted_mean_anomaly"] * np.dot(
            shifted_mean_anomaly_node_vals,
            self.eim_B["shifted_mean_anomaly"][:, start_idx_t:],
        )

        l_start = lt_relation[0]
        start_idx_l = self._find_conservative_starting_truncation_index(
            grid=l_grid_sur, val=l_start
        )

        delta_amp_native = self.norm_factor["res_amp"] * np.dot(
            res_amp_node_vals, self.eim_B["res_amp"][:, start_idx_l:]
        )
        delta_phase_native = self.norm_factor["res_phase"] * np.dot(
            res_phase_node_vals, self.eim_B["res_phase"][:, start_idx_l:]
        )
        res_circ_phase_native = self.norm_factor["res_circ_phase"] * np.dot(
            res_circ_phase_node_vals, self.eim_B["res_circ_phase"][:, start_idx_l:]
        )

        amp_scale = mass_scaling_factor / _amp_correction_factor

        # Fast path: for a uniform output grid the l_s interpolation, the three
        # residual-piece interpolations (sharing one l-grid index/weight), the
        # circular-baseline combination, initial-phase removal, and complex-mode
        # assembly are fused into a single pass. See _fused_ecc_mode_uniform.
        if times is None and not return_orbital_variables:
            tg0, tdg = _uniform_grid_scalars(t_grid_sur, start_idx_t)
            tg0 *= mass_scaling_factor
            tdg *= mass_scaling_factor
            lg0, ldg = _uniform_grid_scalars(l_grid_sur, start_idx_l)
            return new_t_grid, _fused_ecc_mode_uniform(
                tg0,
                tdg,
                lt_relation,
                t_start,
                delta_t,
                num_samples,
                lg0,
                ldg,
                delta_amp_native,
                delta_phase_native,
                res_circ_phase_native,
                amp0_,
                phase0_,
                amp_scale,
                remove_start_phase,
            )

        t_grid_sur = t_grid_sur[start_idx_t:] * mass_scaling_factor
        l_grid_sur = l_grid_sur[start_idx_l:]

        l_s = np.interp(new_t_grid, t_grid_sur, lt_relation)
        delta_A = np.interp(l_s, l_grid_sur, delta_amp_native)
        delta_phi = np.interp(l_s, l_grid_sur, delta_phase_native)
        res_circ_phi = np.interp(l_s, l_grid_sur, res_circ_phase_native)

        delta_A *= amp_scale
        amp = amp0_ + delta_A
        phase = phase0_ + res_circ_phi + delta_phi
        if remove_start_phase:
            phase -= phase[0]

        if return_orbital_variables:
            e_native = self.norm_factor["e"] * np.dot(
                e_node_vals, self.eim_B["e"][:, start_idx_l:]
            )
            e = np.interp(l_s, l_grid_sur, e_native)

            l = l_s + mean_anomaly_offset_of_shifted_mean_anomaly
            l -= (
                2 * np.pi * np.floor(l[0] / (2 * np.pi))
            )  # Bringing starting value of mean anomaly in [0, 2pi)

            x_node_vals = shared_ref_vals["x"]
            x_native = self.norm_factor["x"] * np.dot(
                x_node_vals, self.eim_B["x"][:, start_idx_l:]
            )
            x = np.interp(l_s, l_grid_sur, x_native)

            orb_vars = {"e": e, "l": l, "x": x}
            return new_t_grid, orb_vars, mode_from_amp_phase(amp, phase)

        return new_t_grid, mode_from_amp_phase(amp, phase)


_surrogate_instance = None


def _get_surrogate():
    global _surrogate_instance
    if _surrogate_instance is None:
        sur_data_dir = os.environ.get("ESIGMASUR_DATA_PATH", None)
        if sur_data_dir is None:
            raise RuntimeError(
                "Surrogate data not found. Please set the ESIGMASUR_DATA_PATH "
                "environment variable to the path of the surrogate data directory."
            )
        ecc_sur_data_dir = os.path.join(sur_data_dir, "ecc_sur_data")
        circ_sur_data_dir = os.path.join(sur_data_dir, "circ_sur_data")
        if not os.path.isdir(ecc_sur_data_dir) or not os.path.isdir(circ_sur_data_dir):
            raise RuntimeError(
                f"Surrogate data not found. Please ensure that the environment variable ESIGMASUR_DATA_PATH points to the surrogate data directory."
            )
        _surrogate_instance = EccentricSurrogate(
            ecc_data_dir=ecc_sur_data_dir, circ_data_dir=circ_sur_data_dir
        )
        print("Surrogate data loaded successfully.")
    return _surrogate_instance
