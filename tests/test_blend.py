import numpy as np
import pytest

from esigmapy.blend import (
    blend_modes,
    blend_series,
    compute_amplitude,
    compute_frequency,
    compute_phase,
    find_first_value_location_in_series,
    find_last_value_location_in_series,
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
        # 1.5 lies between arr[1]=1 and arr[2]=2; nearest is equidistant → returns earlier index
        arr = np.array([0.0, 1.0, 2.0, 3.0])
        idx = find_first_value_location_in_series(arr, 1.5)
        assert idx == 1

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
        assert result[-1] == pytest.approx(1.0)  # tau=1 at right edge

    def test_output_length(self):
        x1 = np.ones(100)
        x2 = np.zeros(100)
        t1, t2 = 10, 30
        result = blend_series(x1, x2, t1, t2, t1, t2)
        # blend_series assumes t1 and t2 are inclusive indices;
        # output length should be t2 - t1 + 1
        assert len(result) == t2 - t1 + 1

    def test_monotone_blend_between_constants(self):
        # Blending from 0→1 over a window must be monotone non-decreasing
        x1 = np.zeros(100)
        x2 = np.ones(100)
        result = blend_series(x1, x2, 20, 60, 20, 60)
        assert np.all(np.diff(result) >= -1e-14)

    def test_unequal_window_lengths_raise(self):
        x1 = np.ones(100)
        x2 = np.ones(100)
        with pytest.raises(ValueError):
            blend_series(x1, x2, 10, 30, 10, 25)  # insp window=21, mr window=16

    def test_different_offsets_same_length(self):
        # Windows can be at different positions but must have the same length
        x1 = np.arange(100, dtype=float)
        x2 = np.arange(100, dtype=float) * 2
        result = blend_series(x1, x2, 10, 30, 50, 70)  # both windows length 20+1=21
        assert len(result) == 21
        assert result[0] == pytest.approx(x1[10])
        assert result[-1] == pytest.approx(x2[70])


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


# ---------------------------------------------------------------------------
# blend_modes — output correctness
# ---------------------------------------------------------------------------


@pytest.fixture
def blend_setup():
    """Synthetic linear-chirp inspiral + MR modes with a known blending window.

    (2,2) inspiral: 40→160 Hz GW over 500 samples  (orbital: 20→80 Hz)
    (2,2) MR      : 60→300 Hz GW over 200 samples
    Attach at 120 Hz ± 10 Hz, using avg orbital frequency.
    """
    dt = 1.0 / 4096
    N_insp, N_mr = 500, 200

    def chirp(N, f0, f1):
        t = np.arange(N) * dt
        T = (N - 1) * dt
        phi = 2 * np.pi * (f0 * t + 0.5 * (f1 - f0) / T * t**2)
        return np.exp(-1j * phi)

    insp_22 = chirp(N_insp, 40, 160)
    mr_22 = chirp(N_mr, 60, 300)
    return {
        "inspiral_modes": {(2, 2): insp_22, (2, -2): np.conj(insp_22)},
        "mr_modes": {(2, 2): mr_22, (2, -2): np.conj(mr_22)},
        "orbital_freq": np.linspace(20, 80, N_insp),
        "frq_attach": 120.0,
        "frq_width": 20.0,
        "dt": dt,
    }


def _run_blend(setup, **overrides):
    kwargs = dict(
        modes_to_blend=[(2, 2), (2, -2)],
        mode_to_align_by=(2, 2),
        blend_using_avg_orbital_frequency=True,
        blend_aligning_merger_to_inspiral=True,
        include_conjugate_modes=False,
        delta_t=setup["dt"],
        frq_width=setup["frq_width"],
    )
    kwargs.update(overrides)
    return blend_modes(
        setup["inspiral_modes"],
        setup["mr_modes"],
        setup["orbital_freq"],
        setup["frq_attach"],
        **kwargs,
    )


class TestBlendModesOutput:
    """Correctness tests for blend_modes output values.

    Return tuple layout:
        [0] hybrid_modes dict
        [1] t1_index_insp,  [2] t1_index_mr
        [3] t2_index_insp,  [4] t2_index_mr
        [5] frq_insp_window, [6] frq_hyb_window, [7] frq_mr_window
        [8] amp_insp_window, [9] amp_hyb_window,  [10] amp_mr_window
    """

    def test_return_tuple_length(self, blend_setup):
        assert len(_run_blend(blend_setup)) == 11

    def test_hybrid_keys_match_modes_to_blend(self, blend_setup):
        result = _run_blend(blend_setup)
        assert set(result[0].keys()) == {(2, 2), (2, -2)}

    def test_hybrid_mode_is_ndarray(self, blend_setup):
        result = _run_blend(blend_setup)
        assert isinstance(result[0][(2, 2)], np.ndarray)

    # ── inspiral head / MR tail preservation ─────────────────────────────

    def test_inspiral_head_unchanged_when_aligning_mr(self, blend_setup):
        """Before t1_index_insp the hybrid must equal the inspiral exactly."""
        result = _run_blend(blend_setup, blend_aligning_merger_to_inspiral=True)
        hybrid = result[0][(2, 2)]
        t1_i = result[1]
        np.testing.assert_array_equal(
            hybrid[:t1_i], blend_setup["inspiral_modes"][(2, 2)][:t1_i]
        )

    def test_mr_tail_unchanged_when_aligning_inspiral(self, blend_setup):
        """After t2_index_mr the MR tail must be unchanged (phase shift lands on inspiral head)."""
        result = _run_blend(blend_setup, blend_aligning_merger_to_inspiral=False)
        hybrid = result[0][(2, 2)]
        t2_mr = result[4]
        mr = blend_setup["mr_modes"][(2, 2)]
        np.testing.assert_array_equal(hybrid[-(len(mr) - t2_mr - 1) :], mr[t2_mr + 1 :])

    # ── amplitude correctness ─────────────────────────────────────────────

    def test_window_amplitude_equals_returned_amp_hyb_window(self, blend_setup):
        """Amplitude of hybrid in the window must equal the returned amp_hyb_window."""
        result = _run_blend(blend_setup)
        hybrid = result[0][(2, 2)]
        t1_i, t2_i = result[1], result[3]
        amp_hyb_w = result[9][(2, 2)]
        np.testing.assert_allclose(
            compute_amplitude(hybrid[t1_i : t2_i + 1]), amp_hyb_w, rtol=1e-12
        )

    def test_mr_tail_amplitude_unaffected_by_phasor_shift(self, blend_setup):
        """The constant phasor applied to the MR tail must not change its amplitude."""
        result = _run_blend(blend_setup, blend_aligning_merger_to_inspiral=True)
        hybrid = result[0][(2, 2)]
        t1_i, t2_i, t2_mr = result[1], result[3], result[4]
        mr = blend_setup["mr_modes"][(2, 2)]
        hybrid_tail = hybrid[t1_i + (t2_i - t1_i) + 1 :]
        np.testing.assert_allclose(
            compute_amplitude(hybrid_tail),
            compute_amplitude(mr[t2_mr + 1 :]),
            rtol=1e-12,
        )

    def test_window_amplitude_starts_at_inspiral_ends_at_mr(self, blend_setup):
        """blend_series guarantees τ=0 at left → amp_hyb[0]=amp_insp[t1_i]
        and τ=1 at right → amp_hyb[-1]=amp_mr[t2_mr]."""
        result = _run_blend(blend_setup)
        amp_hyb_w = result[9][(2, 2)]
        amp_insp_w = result[8][(2, 2)]
        amp_mr_w = result[10][(2, 2)]
        assert amp_hyb_w[0] == pytest.approx(amp_insp_w[0], rel=1e-12)
        assert amp_hyb_w[-1] == pytest.approx(amp_mr_w[-1], rel=1e-12)

    # ── phase / frequency correctness ─────────────────────────────────────

    def test_window_phase_is_integral_of_blended_frequency(self, blend_setup):
        """Core test: phase in the blending window must equal the cumulative
        trapezoid of frq_hyb_window plus the inspiral angle at window start.

        This directly validates the numpy-cumsum trapezoid integration that
        replaced scipy.integrate.cumulative_trapezoid.
        """
        result = _run_blend(blend_setup)
        hybrid = result[0][(2, 2)]
        t1_i, t2_i = result[1], result[3]
        frq_hyb_w = result[6][(2, 2)]
        dt = blend_setup["dt"]

        expected = np.empty(len(frq_hyb_w))
        expected[0] = -np.angle(blend_setup["inspiral_modes"][(2, 2)][t1_i])
        expected[1:] = expected[0] + 2 * np.pi * np.cumsum(
            0.5 * (frq_hyb_w[:-1] + frq_hyb_w[1:]) * dt
        )
        np.testing.assert_allclose(
            compute_phase(hybrid[t1_i : t2_i + 1]), expected, rtol=1e-10
        )

    def test_phase_left_boundary_continuous(self, blend_setup):
        """Phase of hybrid at t1_i must match inspiral phase at t1_i (no jump)."""
        result = _run_blend(blend_setup)
        hybrid = result[0][(2, 2)]
        t1_i = result[1]
        insp = blend_setup["inspiral_modes"][(2, 2)]
        # The window phase is initialised from inspiral_angle_at_window_start
        assert -np.angle(hybrid[t1_i]) == pytest.approx(
            -np.angle(insp[t1_i]), abs=1e-12
        )

    def test_mr_tail_is_constant_phasor_shifted(self, blend_setup):
        """The MR tail in the hybrid is the original MR tail times a constant
        phasor — so the ratio hybrid_tail / mr_tail must be constant."""
        result = _run_blend(blend_setup, blend_aligning_merger_to_inspiral=True)
        hybrid = result[0][(2, 2)]
        t1_i, t2_i, t2_mr = result[1], result[3], result[4]
        mr = blend_setup["mr_modes"][(2, 2)]
        hybrid_tail = hybrid[t1_i + (t2_i - t1_i) + 1 :]
        mr_tail = mr[t2_mr + 1 :]
        if len(mr_tail) == 0:
            pytest.skip("no MR tail samples")
        ratio = hybrid_tail / mr_tail
        np.testing.assert_allclose(np.angle(ratio), np.angle(ratio[0]), atol=1e-10)

    def test_inspiral_head_is_constant_phasor_shifted_when_aligning_inspiral(
        self, blend_setup
    ):
        """When the inspiral is phase-shifted to match MR, the head is multiplied
        by a constant phasor — ratio hybrid_head / inspiral_head must be constant."""
        result = _run_blend(blend_setup, blend_aligning_merger_to_inspiral=False)
        hybrid = result[0][(2, 2)]
        t1_i = result[1]
        insp = blend_setup["inspiral_modes"][(2, 2)]
        if t1_i == 0:
            pytest.skip("no inspiral head samples")
        ratio = hybrid[:t1_i] / insp[:t1_i]
        np.testing.assert_allclose(np.angle(ratio), np.angle(ratio[0]), atol=1e-10)

    def test_blended_frequency_starts_at_inspiral_ends_at_mr(self, blend_setup):
        """At τ=0 frq_hyb == frq_insp; at τ=1 frq_hyb == frq_mr."""
        result = _run_blend(blend_setup)
        frq_i_w = result[5][(2, 2)]
        frq_h_w = result[6][(2, 2)]
        frq_mr_w = result[7][(2, 2)]
        assert frq_h_w[0] == pytest.approx(frq_i_w[0], rel=1e-12)
        assert frq_h_w[-1] == pytest.approx(frq_mr_w[-1], rel=1e-12)

    # ── both alignment directions ─────────────────────────────────────────

    def test_align_inspiral_to_mr_produces_valid_output(self, blend_setup):
        """blend_aligning_merger_to_inspiral=False must complete and return modes."""
        result = _run_blend(blend_setup, blend_aligning_merger_to_inspiral=False)
        assert (2, 2) in result[0]
        assert len(result[0][(2, 2)]) > 0

    def test_window_phase_integral_correct_when_aligning_inspiral(self, blend_setup):
        """Same phase-integral check for the align-inspiral-to-MR branch."""
        result = _run_blend(blend_setup, blend_aligning_merger_to_inspiral=False)
        hybrid = result[0][(2, 2)]
        t1_i, t2_i = result[1], result[3]
        frq_hyb_w = result[6][(2, 2)]
        dt = blend_setup["dt"]

        # In this branch, phase_h is offset so the window *ends* at mr_angle_at_window_end
        mr_angle_end = -np.angle(blend_setup["mr_modes"][(2, 2)][result[4]])
        raw_integral = (
            2
            * np.pi
            * np.cumsum(
                np.concatenate([[0], 0.5 * (frq_hyb_w[:-1] + frq_hyb_w[1:]) * dt])
            )
        )
        offset = mr_angle_end - raw_integral[-1]
        expected = raw_integral + offset
        np.testing.assert_allclose(
            compute_phase(hybrid[t1_i : t2_i + 1]), expected, rtol=1e-10
        )
