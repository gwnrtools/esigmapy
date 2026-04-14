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
        t_grid_sur = t_grid_sur[start_idx:] * mass_scaling_factor

        if times is None:
            num_samples = int((t_end - t_start) / delta_t) + 1
            new_t_grid = t_start + np.arange(num_samples) * delta_t
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

        amp = np.interp(new_t_grid, t_grid_sur, amp_native)
        phase = np.interp(new_t_grid, t_grid_sur, phase_native)

        amp *= (
            mass_scaling_factor / _amp_correction_factor
        )  # Correcting for the amplitude normalization factor that was off in the surrogate files.

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
        self.sur_total_mass, self.t_ref, self.t_grid_sur, self.l_grid_sur = (
            self.get_metadata(["M", "t_ref", "t_grid_sur", "l_grid_sur"])
        )

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
    def load_fit(filepath):
        with h5py.File(filepath, "r") as f:
            nodes = f["nodes"][()]
            coeffs = f["coefficients"][()]
        fit = TPI.TP_Interpolant_ND(list(nodes), coeffs=coeffs)
        return fit

    def load_param_space_fits(self, verbose=False):
        for data_piece_name in self.data_piece_names:
            load_dir = os.path.join(self.sur_dir, f"fits/{data_piece_name}_fits")
            # List of all the files which end in .h5
            filenames = [f for f in os.listdir(load_dir) if f.endswith(".h5")]
            sorted_filenames = self._get_sorted_fit_filenames(filenames=filenames)
            self.fit[data_piece_name] = []

            for filename in sorted_filenames:
                filepath = os.path.join(load_dir, filename)
                if verbose:
                    print(f"Loading spline interpolant of GPR from {filename}...")
                fit = self.load_fit(filepath)
                self.fit[data_piece_name].append(fit)

        filepath = os.path.join(
            self.sur_dir, f"fits/mean_anomaly_offset-ref_space-3D-fit_spline.h5"
        )
        fit = self.load_fit(filepath)
        self.fit["mean_anomaly_offset_fit"] = [fit]

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
        t_grid_sur = t_grid_sur[start_idx_t:] * mass_scaling_factor

        if times is None:
            num_samples = int((t_end - t_start) / delta_t) + 1
            new_t_grid = t_start + np.arange(num_samples) * delta_t
        else:
            new_t_grid = times

        q, e_ref, l_ref = params
        self.check_param_range(q=q, e_ref=e_ref, l_ref=l_ref)

        if e_ref > self.e_ref_min:
            amp0_, phase0_ = self.circ_sur(
                M=M,
                q=q,
                delta_t=delta_t,
                t_start=t_start,
                t_end=t_end,
                times=new_t_grid,
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

        e_node_vals = np.asarray(
            [self.fit["e"][i]([eta, e_ref, l_ref]) for i in range(len(self.fit["e"]))]
        )
        e_eim_res_amp = self.norm_factor["e"] * np.dot(
            e_node_vals, self.eim_B["e"][:, self.ei_indices["res_amp"]]
        )
        e_eim_res_phase = self.norm_factor["e"] * np.dot(
            e_node_vals, self.eim_B["e"][:, self.ei_indices["res_phase"]]
        )

        mean_anomaly_offset_of_shifted_mean_anomaly = self.fit[
            "mean_anomaly_offset_fit"
        ][0](
            [eta, e_ref, l_ref]
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

        res_amp_node_vals = np.asarray(
            [
                self.fit["res_amp"][i]([eta, e_eim_res_amp[i], l_eim_res_amp[i]])
                for i in range(len(self.fit["res_amp"]))
            ]
        )
        res_phase_node_vals = np.asarray(
            [
                self.fit["res_phase"][i]([eta, e_eim_res_phase[i], l_eim_res_phase[i]])
                for i in range(len(self.fit["res_phase"]))
            ]
        )
        res_circ_phase_node_vals = np.asarray(
            [
                self.fit["res_circ_phase"][i]([eta, e_ref, l_ref])
                for i in range(len(self.fit["res_circ_phase"]))
            ]
        )
        shifted_mean_anomaly_node_vals = np.asarray(
            [
                self.fit["shifted_mean_anomaly"][i]([eta, e_ref, l_ref])
                for i in range(len(self.fit["shifted_mean_anomaly"]))
            ]
        )

        lt_relation = self.norm_factor["shifted_mean_anomaly"] * np.dot(
            shifted_mean_anomaly_node_vals,
            self.eim_B["shifted_mean_anomaly"][:, start_idx_t:],
        )

        l_start = lt_relation[0]
        start_idx_l = self._find_conservative_starting_truncation_index(
            grid=l_grid_sur, val=l_start
        )
        l_grid_sur = l_grid_sur[start_idx_l:]

        delta_amp_native = self.norm_factor["res_amp"] * np.dot(
            res_amp_node_vals, self.eim_B["res_amp"][:, start_idx_l:]
        )
        delta_phase_native = self.norm_factor["res_phase"] * np.dot(
            res_phase_node_vals, self.eim_B["res_phase"][:, start_idx_l:]
        )
        res_circ_phase_native = self.norm_factor["res_circ_phase"] * np.dot(
            res_circ_phase_node_vals, self.eim_B["res_circ_phase"][:, start_idx_l:]
        )

        l_s = np.interp(new_t_grid, t_grid_sur, lt_relation)
        delta_A = np.interp(l_s, l_grid_sur, delta_amp_native)
        delta_phi = np.interp(l_s, l_grid_sur, delta_phase_native)
        res_circ_phi = np.interp(l_s, l_grid_sur, res_circ_phase_native)

        delta_A *= mass_scaling_factor / _amp_correction_factor
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

            x_node_vals = np.asarray(
                [
                    self.fit["x"][i]([eta, e_ref, l_ref])
                    for i in range(len(self.fit["x"]))
                ]
            )
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
