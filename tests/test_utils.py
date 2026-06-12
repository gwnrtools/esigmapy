import numpy as np
import pytest
import lal
import pycbc.types as pt

from esigmapy.utils import (
    f22_from_x,
    x_from_f22,
    f_ISCO_spin,
    get_peak_freqs,
    get_polarizations_from_multipoles,
)

_MTSUN = lal.MTSUN_SI
_PI = np.pi


# ---------------------------------------------------------------------------
# f22_from_x / x_from_f22
# ---------------------------------------------------------------------------


class TestF22X:
    @pytest.mark.parametrize(
        "x,M",
        [
            (0.05, 20.0),
            (0.1, 50.0),
            (0.2, 100.0),
            (0.3, 200.0),
        ],
    )
    def test_roundtrip_x_to_f_to_x(self, x, M):
        assert x_from_f22(f22_from_x(x, M), M) == pytest.approx(x, rel=1e-12)

    @pytest.mark.parametrize(
        "f,M",
        [
            (20.0, 30.0),
            (50.0, 50.0),
            (100.0, 100.0),
        ],
    )
    def test_roundtrip_f_to_x_to_f(self, f, M):
        assert f22_from_x(x_from_f22(f, M), M) == pytest.approx(f, rel=1e-12)

    def test_f22_mass_scaling(self):
        # Doubling total mass at fixed x halves the GW frequency
        x = 0.1
        assert f22_from_x(x, 40.0) == pytest.approx(f22_from_x(x, 20.0) / 2, rel=1e-12)

    def test_x_mass_scaling(self):
        # At fixed f, doubling M scales x by 2^(2/3)
        f = 100.0
        ratio = x_from_f22(f, 40.0) / x_from_f22(f, 20.0)
        assert ratio == pytest.approx(2.0 ** (2.0 / 3.0), rel=1e-12)

    def test_f22_positive(self):
        assert f22_from_x(0.1, 50.0) > 0

    def test_x_positive(self):
        assert x_from_f22(100.0, 50.0) > 0

    def test_f22_known_value(self):
        x, M = 0.1, 50.0
        expected = x**1.5 / (M * _MTSUN * _PI)
        assert f22_from_x(x, M) == pytest.approx(expected, rel=1e-14)

    def test_x_known_value(self):
        f, M = 100.0, 50.0
        expected = (M * _MTSUN * _PI * f) ** (2.0 / 3.0)
        assert x_from_f22(f, M) == pytest.approx(expected, rel=1e-14)


# ---------------------------------------------------------------------------
# f_ISCO_spin
# ---------------------------------------------------------------------------


class TestFISCOSpin:
    def test_positive_for_zero_spin(self):
        assert f_ISCO_spin(20.0, 20.0, 0.0, 0.0) > 0

    def test_positive_for_nonzero_spin(self):
        assert f_ISCO_spin(10.0, 20.0, 0.3, -0.2) > 0

    def test_mass_scaling_equal_mass_zero_spin(self):
        # For equal-mass, zero-spin: frequency ∝ 1/M (eta fixed, Scap=0, SS=1)
        f1 = f_ISCO_spin(10.0, 10.0, 0.0, 0.0)
        f2 = f_ISCO_spin(20.0, 20.0, 0.0, 0.0)
        assert f1 / f2 == pytest.approx(2.0, rel=1e-6)

    def test_prograde_spin_raises_isco_frequency(self):
        f_zero = f_ISCO_spin(10.0, 10.0, 0.0, 0.0)
        f_pro = f_ISCO_spin(10.0, 10.0, 0.5, 0.5)
        assert f_pro > f_zero

    def test_retrograde_spin_lowers_isco_frequency(self):
        f_zero = f_ISCO_spin(10.0, 10.0, 0.0, 0.0)
        f_retro = f_ISCO_spin(10.0, 10.0, -0.5, -0.5)
        assert f_retro < f_zero

    def test_spin_monotone_in_spin_magnitude(self):
        # Larger prograde spin → higher ISCO frequency
        f_lo = f_ISCO_spin(10.0, 10.0, 0.3, 0.3)
        f_hi = f_ISCO_spin(10.0, 10.0, 0.7, 0.7)
        assert f_hi > f_lo


# ---------------------------------------------------------------------------
# get_peak_freqs
# ---------------------------------------------------------------------------


class TestGetPeakFreqs:
    @staticmethod
    def _ts(data, dt=1.0 / 4096):
        return pt.TimeSeries(data.astype(float), delta_t=dt)

    def test_no_peaks_for_monotone_increasing(self):
        ts = self._ts(np.linspace(10.0, 200.0, 500))
        times, peaks = get_peak_freqs(ts)
        assert len(peaks) == 0
        assert len(times) == 0

    def test_no_peaks_for_constant_series(self):
        ts = self._ts(np.full(100, 42.0))
        times, peaks = get_peak_freqs(ts)
        assert len(peaks) == 0

    def test_single_isolated_peak(self):
        data = np.zeros(101)
        data[50] = 1.0
        ts = self._ts(data)
        times, peaks = get_peak_freqs(ts)
        assert len(peaks) == 1
        assert peaks[0] == pytest.approx(1.0)

    def test_three_isolated_peaks_count(self):
        data = np.zeros(301)
        data[50] = 1.0
        data[150] = 2.0
        data[250] = 3.0
        ts = self._ts(data)
        _, peaks = get_peak_freqs(ts)
        assert len(peaks) == 3

    def test_three_isolated_peaks_values(self):
        data = np.zeros(301)
        data[50] = 1.0
        data[150] = 2.0
        data[250] = 3.0
        ts = self._ts(data)
        _, peaks = get_peak_freqs(ts)
        assert sorted(peaks) == pytest.approx([1.0, 2.0, 3.0])

    def test_boundary_indices_ignored(self):
        # Peaks at index 0 and n-1 are never reported
        data = np.zeros(51)
        data[0] = 10.0
        data[-1] = 10.0
        data[25] = 5.0
        ts = self._ts(data)
        _, peaks = get_peak_freqs(ts)
        assert len(peaks) == 1
        assert peaks[0] == pytest.approx(5.0)

    def test_peak_times_correspond_to_peak_values(self):
        dt = 1.0 / 1024
        data = np.zeros(201)
        data[100] = 7.0
        ts = self._ts(data, dt=dt)
        times, peaks = get_peak_freqs(ts)
        assert len(times) == 1
        assert times[0] == pytest.approx(100 * dt, rel=1e-6)


# ---------------------------------------------------------------------------
# get_polarizations_from_multipoles
# ---------------------------------------------------------------------------


class TestGetPolarizationsFromMultipoles:
    def test_zero_mode_gives_zero_polarizations(self):
        modes = {(2, 2): np.array([0.0 + 0j, 0.0 + 0j])}
        hp, hc = get_polarizations_from_multipoles(
            modes, inclination=0.3, coa_phase=0.0
        )
        assert np.allclose(hp, 0.0)
        assert np.allclose(hc, 0.0)

    def test_linear_amplitude_scaling(self):
        modes1 = {(2, 2): np.array([1.0 + 0j])}
        modes2 = {(2, 2): np.array([3.0 + 0j])}
        hp1, hc1 = get_polarizations_from_multipoles(
            modes1, inclination=0.5, coa_phase=0.0
        )
        hp2, hc2 = get_polarizations_from_multipoles(
            modes2, inclination=0.5, coa_phase=0.0
        )
        assert np.allclose(hp2, 3.0 * hp1)
        assert np.allclose(hc2, 3.0 * hc1)

    def test_polarization_norm_single_mode(self):
        # hp^2 + hc^2 == |h22|^2 * |Y_{-2,22}|^2  for a single mode
        inc, coa = 0.5, 0.3
        h22 = 2.0 + 1.5j
        modes = {(2, 2): np.array([h22])}
        hp, hc = get_polarizations_from_multipoles(
            modes, inclination=inc, coa_phase=coa
        )
        ylm = lal.SpinWeightedSphericalHarmonic(inc, coa, -2, 2, 2)
        expected = abs(h22) ** 2 * abs(ylm) ** 2
        actual = float(hp[0]) ** 2 + float(hc[0]) ** 2
        assert actual == pytest.approx(expected, rel=1e-12)

    def test_two_modes_superpose_linearly(self):
        inc, coa = 0.4, 0.2
        h22 = np.array([1.0 + 0j])
        h21 = np.array([0.5 + 0.3j])
        hp22, hc22 = get_polarizations_from_multipoles({(2, 2): h22}, inc, coa)
        hp21, hc21 = get_polarizations_from_multipoles({(2, 1): h21}, inc, coa)
        hp_both, hc_both = get_polarizations_from_multipoles(
            {(2, 2): h22, (2, 1): h21}, inc, coa
        )
        assert np.allclose(hp_both, hp22 + hp21)
        assert np.allclose(hc_both, hc22 + hc21)

    def test_output_is_real_valued(self):
        modes = {(2, 2): np.array([1.0 + 2.0j, 3.0 - 1.0j])}
        hp, hc = get_polarizations_from_multipoles(
            modes, inclination=0.3, coa_phase=0.0
        )
        assert np.all(np.imag(hp) == 0.0)
        assert np.all(np.imag(hc) == 0.0)

    def test_coa_phase_shifts_polarizations(self):
        # Different coalescence phases must give different results (not trivially zero-effect)
        modes = {(2, 2): np.array([1.0 + 0j])}
        hp0, hc0 = get_polarizations_from_multipoles(
            modes, inclination=0.5, coa_phase=0.0
        )
        hp1, hc1 = get_polarizations_from_multipoles(
            modes, inclination=0.5, coa_phase=0.5
        )
        assert not np.allclose(hp0, hp1) or not np.allclose(hc0, hc1)
