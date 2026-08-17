import numpy as np
import pytest
from scipy import integrate
from unittest.mock import MagicMock, patch

import lal
import pycbc.types as pt

from esigmapy.generator import (
    _get_window_start,
    _get_transition_frequency_window,
    get_imr_esigma_modes,
    get_inspiral_esigma_modes,
    eccentricity_at_reference_frequency,
)

# ---------------------------------------------------------------------------
# Shared helpers for Tier 2 (LAL-mocked) tests
# ---------------------------------------------------------------------------


def _lal_ts(data):
    """Minimal mock of a LAL REAL8TimeSeries."""
    ts = MagicMock()
    ts.data.data = np.asarray(data, dtype=float)
    return ts


def _lal_mode(n, freq_hz=50.0, dt=1.0 / 4096):
    """Minimal mock of a LAL COMPLEX16TimeSeries for a GW mode."""
    m = MagicMock()
    t = np.arange(n) * dt
    m.data.data = np.exp(-2j * np.pi * freq_hz * t)
    return m


def _fake_dynamics(n=64, dt=1.0 / 4096):
    """Return a tuple of 9 mock LAL time series mimicking SimInspiralESIGMADynamics output."""
    t_arr = np.arange(n, dtype=float) * dt
    return (
        _lal_ts(t_arr.copy()),  # t  (copied so in-place *= doesn't affect original)
        _lal_ts(np.ones(n) * 0.1),  # x
        _lal_ts(np.zeros(n)),  # e
        _lal_ts(np.zeros(n)),  # l
        _lal_ts(np.zeros(n)),  # u
        _lal_ts(np.linspace(0, 5, n)),  # phi
        _lal_ts(np.ones(n) * 100.0),  # phidot
        _lal_ts(np.ones(n) * 10.0),  # r
        _lal_ts(np.zeros(n)),  # rdot
    )


# ---------------------------------------------------------------------------
# _get_window_start
# ---------------------------------------------------------------------------


class TestGetWindowStart:
    """Tests for the private helper that integrates a frequency series and
    returns the first index where the cumulative phase exceeds a threshold."""

    _DT = 1.0 / 4096
    _F0 = 100.0  # Hz — constant frequency used in several tests

    def test_forward_meets_threshold(self):
        freq = np.ones(500) * self._F0
        idx = _get_window_start(freq, self._DT, 1.0, direction="forward")
        assert idx is not None
        assert abs(integrate.trapezoid(freq[: idx + 1], dx=self._DT)) >= 1.0

    def test_forward_is_first_crossing(self):
        # The index immediately before should NOT yet meet the threshold
        freq = np.ones(500) * self._F0
        idx = _get_window_start(freq, self._DT, 1.0, direction="forward")
        assert abs(integrate.trapezoid(freq[:idx], dx=self._DT)) < 1.0

    def test_forward_returns_none_when_unreachable(self):
        # Total integral of [1 Hz × 5 samples × dt] is far below 1000 rad
        freq = np.ones(5) * 1.0
        assert _get_window_start(freq, 0.001, 1000.0, direction="forward") is None

    def test_backward_meets_threshold(self):
        freq = np.ones(500) * self._F0
        idx = _get_window_start(freq, self._DT, 1.0, direction="backward")
        assert idx is not None
        assert abs(integrate.trapezoid(freq[idx:], dx=self._DT)) >= 1.0

    def test_backward_is_last_crossing(self):
        # The index one step to the right should NOT meet the threshold alone
        freq = np.ones(500) * self._F0
        idx = _get_window_start(freq, self._DT, 1.0, direction="backward")
        assert abs(integrate.trapezoid(freq[idx + 1 :], dx=self._DT)) < 1.0

    def test_backward_returns_none_when_unreachable(self):
        freq = np.ones(5) * 1.0
        assert _get_window_start(freq, 0.001, 1000.0, direction="backward") is None

    def test_forward_index_in_valid_range(self):
        freq = np.linspace(10.0, 200.0, 500)
        idx = _get_window_start(freq, self._DT, 1.0, direction="forward")
        assert idx is not None
        assert 0 < idx < len(freq)


# ---------------------------------------------------------------------------
# _get_transition_frequency_window
# ---------------------------------------------------------------------------


class TestGetTransitionFrequencyWindow:
    """Tests for the private helper that converts num_hyb_orbits into a
    hybridization frequency window width."""

    @staticmethod
    def _setup(n=2000, dt=1.0 / 4096, f_start=10.0, f_end=200.0):
        freq = np.linspace(f_start, f_end, n)
        phase = np.cumsum(freq) * dt  # monotonically increasing
        return freq, phase, dt

    def test_end_mode_returns_positive_width(self):
        freq, phase, dt = self._setup()
        result = _get_transition_frequency_window(
            phase,
            freq,
            dt,
            f_mr_transition=freq[1000],
            num_hyb_orbits=0.1,
            keep_f_mr_transition_at_center=False,
            blend_using_avg_orbital_frequency=False,
            failsafe=True,
        )
        assert result > 0

    def test_center_mode_returns_positive_width(self):
        freq, phase, dt = self._setup()
        result = _get_transition_frequency_window(
            phase,
            freq,
            dt,
            f_mr_transition=freq[1000],
            num_hyb_orbits=0.1,
            keep_f_mr_transition_at_center=True,
            blend_using_avg_orbital_frequency=False,
            failsafe=True,
        )
        assert result > 0

    def test_avg_orbital_frequency_mode_returns_positive_width(self):
        freq, phase, dt = self._setup()
        result = _get_transition_frequency_window(
            phase,
            freq,
            dt,
            f_mr_transition=freq[1000],
            num_hyb_orbits=0.1,
            keep_f_mr_transition_at_center=False,
            blend_using_avg_orbital_frequency=True,
            failsafe=True,
        )
        assert result > 0

    def test_more_orbits_wider_or_equal_window(self):
        freq, phase, dt = self._setup()
        f_tr = freq[1000]
        w_narrow = _get_transition_frequency_window(
            phase,
            freq,
            dt,
            f_tr,
            num_hyb_orbits=0.1,
            keep_f_mr_transition_at_center=False,
            blend_using_avg_orbital_frequency=False,
            failsafe=True,
        )
        w_wide = _get_transition_frequency_window(
            phase,
            freq,
            dt,
            f_tr,
            num_hyb_orbits=0.5,
            keep_f_mr_transition_at_center=False,
            blend_using_avg_orbital_frequency=False,
            failsafe=True,
        )
        assert w_wide >= w_narrow


# ---------------------------------------------------------------------------
# get_imr_esigma_modes — input validation guards (fire before any LAL call)
# ---------------------------------------------------------------------------


class TestGetImrEsigmaModesValidation:
    _BASE = dict(
        mass1=20.0,
        mass2=20.0,
        f_lower=20.0,
        delta_t=1.0 / 2048,
        merger_ringdown_approximant="NRSur7dq4",
    )

    def test_invalid_approximant_raises_before_lal(self):
        kwargs = {**self._BASE, "merger_ringdown_approximant": "IMRPhenomD"}
        with pytest.raises(IOError):
            get_imr_esigma_modes(**kwargs, mean_anomaly=0.0)

    def test_both_phase_angles_none_raises(self):
        with pytest.raises(IOError):
            get_imr_esigma_modes(**self._BASE, mean_anomaly=None, coa_phase=None)

    def test_align_merger_without_mean_anomaly_raises(self):
        # blend_aligning_merger_to_inspiral=True (default) requires mean_anomaly
        with pytest.raises(IOError):
            get_imr_esigma_modes(
                **self._BASE,
                blend_aligning_merger_to_inspiral=True,
                mean_anomaly=None,
                coa_phase=0.0,
            )

    def test_align_inspiral_without_coa_phase_raises(self):
        # blend_aligning_merger_to_inspiral=False requires coa_phase
        with pytest.raises(IOError):
            get_imr_esigma_modes(
                **self._BASE,
                blend_aligning_merger_to_inspiral=False,
                mean_anomaly=0.0,
                coa_phase=None,
            )


# ---------------------------------------------------------------------------
# get_inspiral_esigma_modes — Tier 2 (LAL mocked)
# ---------------------------------------------------------------------------


class TestGetInspiralEsigmaModesWithMock:
    """Tests for the Python orchestration logic in get_inspiral_esigma_modes.
    lalsimulation is replaced by a MagicMock so no C library is exercised."""

    @pytest.fixture
    def ctx(self):
        n, dt = 64, 1.0 / 4096
        dyn = _fake_dynamics(n, dt)
        mode = _lal_mode(n, 50.0, dt)
        return n, dt, dyn, mode

    def _call(self, mock_ls, dyn, mode, dt, **kwargs):
        mock_ls.SimInspiralESIGMADynamics.return_value = dyn
        mock_ls.SimInspiralESIGMAModeFromDynamics.return_value = mode
        defaults = dict(
            mass1=10.0,
            mass2=10.0,
            f_lower=20.0,
            delta_t=dt,
            modes_to_use=[(2, 2)],  # always explicit — avoids mutating the default
            include_conjugate_modes=False,
        )
        return get_inspiral_esigma_modes(**{**defaults, **kwargs})

    def test_returned_modes_contain_requested_key(self, ctx):
        n, dt, dyn, mode = ctx
        with patch("esigmapy.generator.ls") as mock_ls:
            modes = self._call(mock_ls, dyn, mode, dt)
        assert (2, 2) in modes

    def test_conjugate_modes_added_when_flag_true(self, ctx):
        n, dt, dyn, mode = ctx
        with patch("esigmapy.generator.ls") as mock_ls:
            modes = self._call(mock_ls, dyn, mode, dt, include_conjugate_modes=True)
        assert (2, -2) in modes

    def test_conjugate_modes_absent_when_flag_false(self, ctx):
        n, dt, dyn, mode = ctx
        with patch("esigmapy.generator.ls") as mock_ls:
            modes = self._call(mock_ls, dyn, mode, dt, include_conjugate_modes=False)
        assert (2, -2) not in modes

    def test_returns_pycbc_timeseries_by_default(self, ctx):
        n, dt, dyn, mode = ctx
        with patch("esigmapy.generator.ls") as mock_ls:
            modes = self._call(mock_ls, dyn, mode, dt)
        assert isinstance(modes[(2, 2)], pt.TimeSeries)

    def test_returns_numpy_array_when_not_pycbc(self, ctx):
        n, dt, dyn, mode = ctx
        with patch("esigmapy.generator.ls") as mock_ls:
            t_arr, modes = self._call(
                mock_ls, dyn, mode, dt, return_pycbc_timeseries=False
            )
        assert isinstance(modes[(2, 2)], np.ndarray)

    def test_distance_converted_from_mpc_to_si(self, ctx):
        n, dt, dyn, mode = ctx
        distance_mpc = 200.0
        with patch("esigmapy.generator.ls") as mock_ls:
            self._call(mock_ls, dyn, mode, dt, distance=distance_mpc)
            passed_distance = mock_ls.SimInspiralESIGMAModeFromDynamics.call_args.args[
                -1
            ]
        assert passed_distance == pytest.approx(
            distance_mpc * 1e6 * lal.PC_SI, rel=1e-12
        )

    def test_mode_generation_called_once_per_mode(self, ctx):
        n, dt, dyn, mode = ctx
        with patch("esigmapy.generator.ls") as mock_ls:
            self._call(
                mock_ls,
                dyn,
                mode,
                dt,
                modes_to_use=[(2, 2), (3, 3)],
                include_conjugate_modes=False,
            )
        assert mock_ls.SimInspiralESIGMAModeFromDynamics.call_count == 2


# ---------------------------------------------------------------------------
# eccentricity_at_reference_frequency — Tier 2 (LAL mocked)
# ---------------------------------------------------------------------------


class TestEccentricityAtReferenceFrequencyWithMock:
    """The function finds the x index closest to x_reference and returns e there.
    Tests verify that index-finding logic and the return value."""

    @staticmethod
    def _setup(n=100, target_idx=60, M=20.0, dt=1.0 / 4096):
        x_arr = np.linspace(0.05, 0.5, n)
        e_arr = np.linspace(0.20, 0.01, n)
        # Choose f_reference so x_reference == x_arr[target_idx] exactly
        f_ref = x_arr[target_idx] ** 1.5 / (np.pi * M * lal.MTSUN_SI)
        dyn = (
            _lal_ts(np.arange(n, dtype=float) * dt),
            _lal_ts(x_arr),
            _lal_ts(e_arr),
            _lal_ts(np.zeros(n)),
            _lal_ts(np.zeros(n)),
            _lal_ts(np.zeros(n)),
            _lal_ts(np.ones(n) * 100.0),
            _lal_ts(np.ones(n) * 10.0),
            _lal_ts(np.zeros(n)),
        )
        return dyn, x_arr, e_arr, f_ref

    def test_returns_eccentricity_at_x_reference(self):
        M, target_idx = 20.0, 60
        dyn, x_arr, e_arr, f_ref = self._setup(target_idx=target_idx, M=M)
        with patch("esigmapy.generator.ls") as mock_ls:
            mock_ls.SimInspiralESIGMADynamics.return_value = dyn
            result = eccentricity_at_reference_frequency(
                M / 2, M / 2, 0, 0, 0.1, 0, 20.0, 4096, f_ref
            )
        assert result == pytest.approx(e_arr[target_idx], rel=1e-10)

    def test_returns_float(self):
        M = 20.0
        dyn, _, _, f_ref = self._setup(M=M)
        with patch("esigmapy.generator.ls") as mock_ls:
            mock_ls.SimInspiralESIGMADynamics.return_value = dyn
            result = eccentricity_at_reference_frequency(
                M / 2, M / 2, 0, 0, 0.1, 0, 20.0, 4096, f_ref
            )
        assert np.isscalar(result) or result.ndim == 0
