import numpy as np
import pytest

from esigmapy.blend import (
    align_in_phase,
    blend_modes,
    blend_series,
    compute_amplitude,
    compute_frequency,
    compute_phase,
    find_first_value_location_in_series,
    find_last_value_location_in_series,
    mismatch_discrete,
)

# ---------------------------------------------------------------------------
# find_first_value_location_in_series
# ---------------------------------------------------------------------------


class TestFindFirstValue:
    def test_below_min_raises(self):
        with pytest.raises(Exception, match="[Ll]ower"):
            find_first_value_location_in_series(np.array([1.0, 2.0, 3.0]), 0.5)

    def test_above_max_raises(self):
        with pytest.raises(Exception, match="[Hh]igher"):
            find_first_value_location_in_series(np.array([1.0, 2.0, 3.0]), 3.5)

    def test_returns_first_bracket_index(self):
        # 1.5 lies between arr[1]=1 and arr[2]=2; nearest is equidistant → either valid
        arr = np.array([0.0, 1.0, 2.0, 3.0])
        idx = find_first_value_location_in_series(arr, 1.5)
        assert idx in (1, 2)

    def test_nearer_left_element_chosen(self):
        # 1.2 is closer to arr[1]=1.0 than arr[2]=2.0
        arr = np.array([0.0, 1.0, 2.0, 3.0])
        idx = find_first_value_location_in_series(arr, 1.2)
        assert idx == 1

    def test_nearer_right_element_chosen(self):
        # 1.8 is closer to arr[2]=2.0 than arr[1]=1.0
        arr = np.array([0.0, 1.0, 2.0, 3.0])
        idx = find_first_value_location_in_series(arr, 1.8)
        assert idx == 2

    def test_exact_match_at_start(self):
        arr = np.array([1.0, 2.0, 3.0, 4.0])
        idx = find_first_value_location_in_series(arr, 1.0)
        assert idx == 0

    def test_first_occurrence_returned_for_non_monotone(self):
        # arr crosses 0.7 between idx 0→1 (first) and again between 2→3
        arr = np.array([0.0, 1.0, 0.5, 1.0, 2.0])
        idx = find_first_value_location_in_series(arr, 0.7)
        # First bracket is [0→1]; nearest of 0.0 and 1.0 to 0.7 is 1.0 → idx 1
        assert idx <= 1

    def test_result_within_valid_range(self):
        arr = np.linspace(1.0, 5.0, 50)
        idx = find_first_value_location_in_series(arr, 3.0)
        assert 0 <= idx < len(arr)


class TestFindLastValue:
    def test_below_min_raises(self):
        with pytest.raises(Exception, match="[Ll]ower"):
            find_last_value_location_in_series(np.array([1.0, 2.0, 3.0]), 0.5)

    def test_above_max_raises(self):
        with pytest.raises(Exception, match="[Hh]igher"):
            find_last_value_location_in_series(np.array([1.0, 2.0, 3.0]), 3.5)

    def test_monotone_returns_index_at_least_as_large_as_first(self):
        arr = np.linspace(0.0, 10.0, 100)
        target = 5.0
        idx_first = find_first_value_location_in_series(arr, target)
        idx_last = find_last_value_location_in_series(arr, target)
        assert idx_last >= idx_first

    def test_returns_later_occurrence_for_non_monotone(self):
        # Two crossings of 1.5 in arr: between idx 1-2 and between idx 3-4
        arr = np.array([0.0, 1.0, 2.0, 1.0, 2.0, 3.0])
        idx = find_last_value_location_in_series(arr, 1.5)
        # The last crossing is in the 3-4 region; result should be >= 3
        assert idx >= 3

    def test_result_within_valid_range(self):
        arr = np.linspace(1.0, 5.0, 50)
        idx = find_last_value_location_in_series(arr, 3.0)
        assert 0 <= idx < len(arr)

    def test_returns_last_bracket_index(self):
        # 1.5 lies between arr[1]=1.0 and arr[2]=2.0; nearest is equidistant → returns later index
        arr = np.array([0.0, 1.0, 2.0, 3.0])
        idx = find_last_value_location_in_series(arr, 1.5)
        assert idx == 2

    def test_nearer_left_element_chosen(self):
        # 1.2 is closer to arr[1]=1.0 than arr[2]=2.0
        arr = np.array([0.0, 1.0, 2.0, 3.0])
        idx = find_last_value_location_in_series(arr, 1.2)
        assert idx == 1

    def test_nearer_right_element_chosen(self):
        # 1.8 is closer to arr[2]=2.0 than arr[1]=1.0
        arr = np.array([0.0, 1.0, 2.0, 3.0])
        idx = find_last_value_location_in_series(arr, 1.8)
        assert idx == 2

    def test_exact_match_at_start(self):
        arr = np.array([1.0, 2.0, 3.0, 4.0])
        idx = find_last_value_location_in_series(arr, 1.0)
        assert idx == 0

    def test_exact_match_at_end(self):
        arr = np.array([1.0, 2.0, 3.0, 4.0])
        idx = find_last_value_location_in_series(arr, 4.0)
        assert idx == 3


# ---------------------------------------------------------------------------
# compute_amplitude
# ---------------------------------------------------------------------------


class TestComputeAmplitude:
    def test_known_magnitude(self):
        # |3 + 4j| == 5
        wave = np.full(10, 3.0 + 4.0j)
        amp = compute_amplitude(wave)
        assert np.allclose(amp, 5.0, rtol=1e-15, atol=1e-15)

    def test_real_valued_input(self):
        wave = np.array([-3.0, -1.0, 0.0, 1.0, 3.0])
        assert np.allclose(
            compute_amplitude(wave), np.abs(wave), rtol=1e-15, atol=1e-15
        )

    def test_zero_amplitude(self):
        wave = np.zeros(5, dtype=complex)
        assert np.allclose(compute_amplitude(wave), 0.0, rtol=1e-15, atol=1e-15)

    def test_output_nonnegative(self):
        rng = np.random.default_rng(0)
        wave = rng.standard_normal(50) + 1j * rng.standard_normal(50)
        assert np.all(compute_amplitude(wave) >= 0)

    def test_unit_circle(self):
        # Points on the unit circle all have amplitude 1
        phases = np.linspace(0, 2 * np.pi, 100, endpoint=False)
        wave = np.exp(1j * phases)
        assert np.allclose(compute_amplitude(wave), 1.0, rtol=1e-15, atol=1e-15)


# ---------------------------------------------------------------------------
# compute_phase
# ---------------------------------------------------------------------------


class TestComputePhase:
    def test_linear_phase_recovered(self):
        # waveform = exp(-i * phi) → compute_phase returns phi (GW convention)
        N, dt, f = 512, 1.0 / 4096, 80.0
        t = np.arange(N) * dt
        phi = 2 * np.pi * f * t
        wave = np.exp(-1j * phi)
        recovered = compute_phase(wave)
        # Allow for unwrapping edge at boundaries; check interior
        assert np.allclose(recovered[1:-1], phi[1:-1], rtol=1e-15, atol=1e-14)

    def test_constant_phase_is_flat(self):
        # All samples at the same phase → output is a constant array
        wave = np.exp(-1j * 0.7) * np.ones(50)
        phase = compute_phase(wave)
        assert np.allclose(phase, phase[0], rtol=1e-15, atol=1e-15)

    def test_opposite_sign_convention(self):
        # exp(+i*phi) → compute_phase returns -phi (negative sign in formula)
        N, dt, f = 256, 1.0 / 4096, 50.0
        t = np.arange(N) * dt
        phi = 2 * np.pi * f * t
        wave = np.exp(1j * phi)
        recovered = compute_phase(wave)
        assert np.allclose(recovered[1:-1], -phi[1:-1], rtol=1e-15, atol=1e-14)

    def test_output_unwrapped(self):
        # Phase must not have jumps larger than π between consecutive samples
        N, dt, f = 1024, 1.0 / 4096, 200.0
        t = np.arange(N) * dt
        wave = np.exp(-1j * 2 * np.pi * f * t)
        phase = compute_phase(wave)
        assert np.all(np.abs(np.diff(phase)) < np.pi)


# ---------------------------------------------------------------------------
# compute_frequency
# ---------------------------------------------------------------------------


class TestComputeFrequency:
    def test_constant_frequency_recovered(self):
        N, dt, f0 = 512, 1.0 / 4096, 100.0
        t = np.arange(N) * dt
        phase = 2 * np.pi * f0 * t  # linear phase
        freq = compute_frequency(phase, dt)
        # np.gradient is exact for linear functions at all points
        assert np.allclose(freq, f0, rtol=1e-13, atol=1e-15)

    def test_negative_frequency_for_decreasing_phase(self):
        N, dt, f0 = 256, 1.0 / 4096, 50.0
        t = np.arange(N) * dt
        phase = -2 * np.pi * f0 * t  # decreasing phase
        freq = compute_frequency(phase, dt)
        assert np.allclose(freq, -f0, rtol=1e-13, atol=1e-15)

    def test_output_length_matches_input(self):
        phase = np.linspace(0, 10 * np.pi, 200)
        freq = compute_frequency(phase, 1.0 / 4096)
        assert len(freq) == len(phase)

    def test_roundtrip_phase_to_frequency(self):
        # compute_frequency(compute_phase(exp(-i*2*pi*f*t)), dt) ≈ f
        N, dt, f0 = 1024, 1.0 / 4096, 120.0
        t = np.arange(N) * dt
        wave = np.exp(-1j * 2 * np.pi * f0 * t)
        freq = compute_frequency(compute_phase(wave), dt)
        # Ignore first and last points where gradient uses one-sided differences
        assert np.allclose(freq[1:-1], f0, rtol=1e-12, atol=1e-15)


# ---------------------------------------------------------------------------
# mismatch_discrete
# ---------------------------------------------------------------------------


class TestMismatchDiscrete:
    @staticmethod
    def _all_idx(n):
        return np.arange(n)

    def test_identical_waveforms_zero_mismatch(self):
        w = np.exp(-1j * np.linspace(0, 4 * np.pi, 50))
        idx = self._all_idx(len(w))
        assert mismatch_discrete(w, w, idx, idx) == pytest.approx(0.0, abs=1e-14)

    def test_quadrature_phase_gives_mismatch_one(self):
        # w2 = i * w1  →  |w1 - w2|^2 = 2|w1|^2  →  mm = 0.5 * 2 = 1.0
        w1 = np.ones(20, dtype=complex)
        w2 = 1j * w1
        idx = self._all_idx(len(w1))
        assert mismatch_discrete(w1, w2, idx, idx) == pytest.approx(1.0, rel=1e-12)

    def test_antiphase_gives_mismatch_two(self):
        # w2 = -w1  →  |w1 - w2|^2 = 4|w1|^2  →  mm = 0.5 * 4 = 2.0
        w1 = np.ones(20, dtype=complex)
        w2 = -w1
        idx = self._all_idx(len(w1))
        assert mismatch_discrete(w1, w2, idx, idx) == pytest.approx(2.0, rel=1e-12)

    def test_nonnegative(self):
        rng = np.random.default_rng(1)
        w1 = rng.standard_normal(30) + 1j * rng.standard_normal(30)
        w2 = rng.standard_normal(30) + 1j * rng.standard_normal(30)
        idx = self._all_idx(len(w1))
        assert mismatch_discrete(w1, w2, idx, idx) >= 0

    def test_subsample_indices(self):
        # Using every other sample should still give 0 for identical waveforms
        w = np.exp(-1j * np.linspace(0, 2 * np.pi, 40))
        idx = np.arange(0, 40, 2)
        assert mismatch_discrete(w, w, idx, idx) == pytest.approx(0.0, abs=1e-14)


# ---------------------------------------------------------------------------
# blend_series
# ---------------------------------------------------------------------------


class TestBlendSeries:
    def test_left_edge_equals_x1(self):
        x1 = np.ones(100)
        x2 = np.zeros(100)
        result = blend_series(x1, x2, 20, 40, 20, 40)
        assert result[0] == pytest.approx(1.0)  # tau=0 at left edge

    def test_right_edge_approaches_x2(self):
        x1 = np.zeros(100)
        x2 = np.ones(100)
        result = blend_series(x1, x2, 20, 40, 20, 40)
        # tau at last point = sin^2(pi/2 * (n-1)/n), very close to 1
        assert result[-1] > 0.97

    def test_output_length(self):
        x1 = np.ones(100)
        x2 = np.zeros(100)
        t1, t2 = 10, 30
        result = blend_series(x1, x2, t1, t2, t1, t2)
        assert len(result) == t2 - t1

    def test_monotone_blend_between_constants(self):
        # Blending from 0→1 over a window must be monotone non-decreasing
        x1 = np.zeros(100)
        x2 = np.ones(100)
        result = blend_series(x1, x2, 20, 60, 20, 60)
        assert np.all(np.diff(result) >= -1e-14)

    def test_unequal_window_lengths_raise(self):
        x1 = np.ones(100)
        x2 = np.ones(100)
        with pytest.raises(AssertionError):
            blend_series(x1, x2, 10, 30, 10, 25)  # insp window=20, mr window=15

    def test_different_offsets_same_length(self):
        # Windows can be at different positions but must have the same length
        x1 = np.arange(100, dtype=float)
        x2 = np.arange(100, dtype=float) * 2
        result = blend_series(x1, x2, 10, 30, 50, 70)  # both windows length 20
        assert len(result) == 20
        assert result[0] == pytest.approx(x1[10])


# ---------------------------------------------------------------------------
# align_in_phase
# ---------------------------------------------------------------------------


class TestAlignInPhase:
    @staticmethod
    def _make_chirp(n=200, f0=50.0, dt=1.0 / 4096):
        t = np.arange(n) * dt
        return np.exp(-1j * 2 * np.pi * f0 * t)

    @staticmethod
    def _sample_idx(t1, t2, no_sp=8):
        return np.linspace(t1, t2, no_sp, dtype=int) - t1

    def test_empty_inspiral_raises(self):
        with pytest.raises(IOError, match="[Zz]ero length"):
            align_in_phase(
                np.array([]),
                np.ones(10, dtype=complex),
                np.arange(5),
                np.arange(5),
                0,
                4,
                0,
                4,
            )

    def test_narrow_window_raises(self):
        with pytest.raises(IOError):
            align_in_phase(
                np.ones(10, dtype=complex),
                np.ones(10, dtype=complex),
                np.arange(1),
                np.arange(1),
                5,
                4,
                5,
                4,  # t2 < t1: zero-width window
            )

    def test_identical_waveforms_zero_shift(self):
        wave = self._make_chirp()
        t1, t2 = 50, 100
        idx = self._sample_idx(t1, t2)
        _, _, shift = align_in_phase(wave, wave, idx, idx, t1, t2, t1, t2, m_mode=2)
        # Optimal shift for identical waveforms is 0 (mod 2π/m)
        assert np.cos(2 * shift) == pytest.approx(1.0, abs=1e-3)

    def test_reduces_mismatch(self):
        wave = self._make_chirp()
        delta = np.pi / 4  # deliberate phase offset
        mr = wave * np.exp(1j * delta)
        t1, t2 = 50, 100
        idx = self._sample_idx(t1, t2)

        mm_before = mismatch_discrete(wave[t1 : t2 + 1], mr[t1 : t2 + 1], idx, idx)
        insp_a, mr_a, _ = align_in_phase(
            wave,
            mr,
            idx,
            idx,
            t1,
            t2,
            t1,
            t2,
            m_mode=2,
            align_merger_to_inspiral=True,
        )
        mm_after = mismatch_discrete(insp_a[t1 : t2 + 1], mr_a[t1 : t2 + 1], idx, idx)
        assert mm_after < mm_before

    def test_align_merger_to_inspiral_order(self):
        wave = self._make_chirp()
        mr = wave * np.exp(1j * 0.3)
        t1, t2 = 30, 80
        idx = self._sample_idx(t1, t2)
        insp_out, mr_out, _ = align_in_phase(
            wave,
            mr,
            idx,
            idx,
            t1,
            t2,
            t1,
            t2,
            align_merger_to_inspiral=True,
        )
        # Inspiral is unchanged; merger-ringdown is shifted toward inspiral
        assert np.allclose(insp_out, wave, rtol=1e-15, atol=1e-15)
        assert not np.allclose(mr_out, mr, rtol=1e-15, atol=1e-15)  # mr was modified
        assert np.allclose(mr_out, wave, rtol=1e-15, atol=1e-7)  # and now matches wave

    def test_align_inspiral_to_merger_order(self):
        wave = self._make_chirp()
        mr = wave * np.exp(1j * 0.3)
        t1, t2 = 30, 80
        idx = self._sample_idx(t1, t2)
        insp_out, mr_out, _ = align_in_phase(
            wave,
            mr,
            idx,
            idx,
            t1,
            t2,
            t1,
            t2,
            align_merger_to_inspiral=False,
        )
        # Merger-ringdown is unchanged; inspiral is shifted
        assert np.allclose(mr_out, mr, rtol=1e-15, atol=1e-15)
        assert not np.allclose(insp_out, wave, rtol=1e-15, atol=1e-15)


# ---------------------------------------------------------------------------
# blend_modes — input validation guards (pure Python, no LAL required)
# ---------------------------------------------------------------------------


class TestBlendModesValidation:
    @staticmethod
    def _modes(n=500, freq_hz=50.0, dt=1.0 / 4096):
        t = np.arange(n) * dt
        mode = np.exp(-2j * np.pi * freq_hz * t)
        return {(2, 2): mode, (2, -2): np.conj(mode)}

    def test_negative_frq_width_raises(self):
        modes = self._modes()
        with pytest.raises(IOError, match="negative"):
            blend_modes(modes, modes, np.ones(500), 50.0, frq_width=-5.0)

    def test_zero_frq_width_raises(self):
        modes = self._modes()
        with pytest.raises(IOError):
            blend_modes(modes, modes, np.ones(500), 50.0, frq_width=0.0)

    def test_mismatched_orbital_freq_length_raises(self):
        modes = self._modes(n=500)
        with pytest.raises(IOError):
            blend_modes(
                modes,
                modes,
                inspiral_orbital_frequency=np.ones(100),  # wrong length
                frq_attach=50.0,
                frq_width=5.0,
                blend_using_avg_orbital_frequency=True,
            )

    def test_mode_missing_from_inspiral_raises(self):
        insp = {(2, 2): np.ones(500, dtype=complex)}
        mr = {(2, 2): np.ones(500, dtype=complex), (3, 3): np.ones(500, dtype=complex)}
        with pytest.raises(IOError):
            blend_modes(
                insp,
                mr,
                np.ones(500),
                50.0,
                frq_width=5.0,
                modes_to_blend=[(2, 2), (3, 3)],
                include_conjugate_modes=False,
            )

    def test_mode_missing_from_mr_raises(self):
        insp = {
            (2, 2): np.ones(500, dtype=complex),
            (3, 3): np.ones(500, dtype=complex),
        }
        mr = {(2, 2): np.ones(500, dtype=complex)}
        with pytest.raises(IOError):
            blend_modes(
                insp,
                mr,
                np.ones(500),
                50.0,
                frq_width=5.0,
                modes_to_blend=[(2, 2), (3, 3)],
                include_conjugate_modes=False,
            )
