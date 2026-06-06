# Copyright (C) 2023 Kartikey Sharma, Prayush Kumar, Akash Maurya
#
"""Master function to hybridise any complex timeseries using the 'frequency' as a user input,
specifically to be used for gravitational waveform hybridisation, fine-tuned for a single mode.
"""

import numpy as np
from scipy.integrate import cumulative_trapezoid


# Akash: Can use @njit here, but numba is not a dependency of default ESIGMAPy currently
def find_first_value_location_in_series(frq_timeseries, frq_desired):
    if frq_desired < np.min(frq_timeseries):
        raise Exception("Desired frequency out of bounds, lower than min frequency")

    if frq_desired > np.max(frq_timeseries):
        raise Exception("Desired frequency out of bounds, higher than max frequency")

    """ 
    We traverse the array to find the location where the i_th value is less than
    the desired value while the i+1_th value is more, hence locating the desired 
    value somewhere between those two points. We then choose the value closer to 
    the value desired (among i and i+1) and call it the location of the desired value. 
    """
    for idx in range(len(frq_timeseries) - 1):
        if frq_timeseries[idx] <= frq_desired <= frq_timeseries[idx + 1]:
            fr1 = frq_timeseries[idx]
            fr2 = frq_timeseries[idx + 1]
            # If equidistant, <= ensures we select the earlier index (idx)
            return idx if abs(frq_desired - fr1) <= abs(frq_desired - fr2) else idx + 1
    raise ValueError("frq_desired not found in frq_timeseries")


# Akash: Can use @njit here, but numba is not a dependency of default ESIGMAPy currently
def find_last_value_location_in_series(frq_timeseries, frq_desired):
    if frq_desired < np.min(frq_timeseries):
        raise Exception(
            f"""Desired value {frq_desired} out of bounds, lower than min value {np.min(frq_timeseries)}"""
        )

    if frq_desired > np.max(frq_timeseries):
        raise Exception(
            f"""Desired value {frq_desired} out of bounds, higher than max value {np.max(frq_timeseries)}"""
        )

    """ 
    We reverse the array and traverse it to find the location where the i_th value is more than
    the desired value while the i+1_th value is less, hence locating the desired value somewhere
    between those two points. We then choose the value closer to the value desired (among i and i+1) 
    and call it the location of the desired value. 
    """
    reversed_freq_timeseries = frq_timeseries[::-1]
    final_idx = None

    for idx in range(len(reversed_freq_timeseries) - 1):
        # This assumes that the frequency is monotonically increasing with time,
        # which is true only for orbit-averaged frequency for eccentric waveforms
        if (
            reversed_freq_timeseries[idx]
            >= frq_desired
            >= reversed_freq_timeseries[idx + 1]
        ):
            fr1 = reversed_freq_timeseries[idx]
            fr2 = reversed_freq_timeseries[idx + 1]
            final_idx = (
                idx if abs(frq_desired - fr1) <= abs(frq_desired - fr2) else idx + 1
            )
            break
    if final_idx is None:
        raise ValueError("frq_desired not found in frq_timeseries")
    return len(frq_timeseries) - 1 - final_idx


def blend_series(x1, x2, t1_index_insp, t2_index_insp, t1_index_mr, t2_index_mr):
    """
    Function to blend two series x1 and x2 over the window defined by t1_index_insp, t2_index_insp for x1 and t1_index_mr, t2_index_mr for x2.

    Both indices are assumed to be inclusive, i.e. the blending window includes the values at
    t1_index_insp and t2_index_insp for x1 and t1_index_mr and t2_index_mr for x2. Therefore,
    x1 should have data in x1[t1_index_insp:t2_index_insp+1] and x2 should have data in x2[t1_index_mr:t2_index_mr+1].
    The blended series is then given by (1 - tau) * x1 + tau * x2, where tau is the blending function,
    and tau = 0 at t1_index_insp and tau = 1 at t2_index_insp, with a smooth transition in between.
    """

    insp_index_range = t2_index_insp - t1_index_insp
    mr_index_range = t2_index_mr - t1_index_mr

    if insp_index_range <= 0 or mr_index_range <= 0:
        raise ValueError(
            "Invalid index range for blending. Ensure that t2_index > t1_index for both inspiral and merger-ringdown."
        )

    if insp_index_range != mr_index_range:
        raise ValueError("Inconsistent indices passed to blending function")

    # blending fn is an array
    blfn_var = np.arange(
        t1_index_insp, t2_index_insp + 1
    )  # +1 since the window is inclusive of t2_index_insp
    tau = np.square(
        (np.sin((np.pi / 2) * (blfn_var - t1_index_insp) / insp_index_range))
    )

    x_hyb = (1 - tau) * x1[t1_index_insp : t2_index_insp + 1] + tau * x2[
        t1_index_mr : t2_index_mr + 1
    ]
    return x_hyb


def compute_amplitude(waveform):
    amplitude = np.abs(waveform)
    return amplitude


def compute_phase(waveform):
    phase = np.unwrap(-np.angle(waveform))
    return phase


def compute_frequency(phase, delta_t):
    frequency = np.gradient(phase, delta_t) / (2 * np.pi)
    return frequency


def blend_modes(
    inspiral_modes,
    merger_ringdown_modes,
    inspiral_orbital_frequency,
    frq_attach,
    frq_width=10.0,
    delta_t=1.0 / 4096,
    modes_to_blend=[(2, 2), (3, 3), (4, 4)],
    mode_to_align_by=(2, 2),
    blend_using_avg_orbital_frequency=True,
    blend_aligning_merger_to_inspiral=True,
    include_conjugate_modes=True,
    verbose=False,
):
    """blend inspiral and merger-ringdown modes

    Inputs
    ------
    inspiral_modes: dict
        Dictionary indexed by (l, m) containing numpy-like arrays of
        complex-valued mode timeseries.
    merger_ringdown_modes: dict
        Dictionary indexed by (l, m) containing numpy-like arrays of
        complex-valued mode timeseries.

    inspiral_orbital_frequency: 1D NumPy array
        Orbit-averaged orbital frequency timeseries corresponding to
        the inspiral. To be used for finding hybridization window in
        time if `blend_using_avg_orbital_frequency` is set to True.
        Should be a NumPy array of the same length as the inspiral modes.
    frq_attach: float
        Frequency (Hz) at which to align the inspiral and merger-ringdown modes.
        This is understood to be the frequency corresponding to `mode_to_align_by`
        (see below).
    frq_width: {10.0, float}
        Frequency (Hz) window around the central attachment frequency over which
        hybridization of modes is performed.
        This is understood to be the frequency width in terms of the frequency of
        `mode_to_align_by` (see below).
    delta_t: {1/4096, float}
        time grid-spacing for timeseries (in seconds)

    modes_to_blend: {[(2, 2), (3, 3), (4, 4)], list}
        List of modes as tuples of (l, m) values to blend
    mode_to_align_by: {(2, 2), tuple}
        One specific mode (l, m) value that is to be treated as baseline for
        time/phase alignment. Should be a +ve m mode.
        We recommend using only the (2, 2) mode for this.
    blend_using_avg_orbital_frequency: (True, bool)
        Flag that enables use of orbit averaged orbital frequency for all
        calculations of inspiral to merger-ringdown attachment
    blend_aligning_merger_to_inspiral: (True, bool)
        Flag that ensures any phase difference between inspiral and merger-rd
        modes is removed by phase shifting the merger-rd portion. If set to
        False, the inspiral portion is phase shifted instead
    include_conjugate_modes: {True, bool}
        When set to True, we also consider (l, -m) modes in addition to (l, m) ones.
    verbose: {False, bool}
        Set this to True to enable logging output.
    """
    if frq_width <= 0:
        raise IOError(f"""You are trying to blend over a frequency window of
            negative length (= {frq_width}Hz). Fix this.""")
    if blend_using_avg_orbital_frequency:
        if len(inspiral_orbital_frequency) != len(inspiral_modes[mode_to_align_by]):
            raise IOError(
                f"""You asked for hybridization using orbital frequency, but
                the orbital frequency array and inspiral modes array have
                different lengths: {len(inspiral_orbital_frequency)}, {len(inspiral_modes[mode_to_align_by])}"""
            )
    modes_not_aligned_by = modes_to_blend.copy()
    if mode_to_align_by in modes_not_aligned_by:
        modes_not_aligned_by.remove(mode_to_align_by)
    if include_conjugate_modes:
        for el, em in modes_to_blend.copy():
            if (el, -em) not in modes_to_blend:
                modes_to_blend.append((el, -em))
        for el, em in modes_not_aligned_by.copy():
            if (el, -em) not in modes_not_aligned_by:
                modes_not_aligned_by.append((el, -em))

    # Input checks
    for lm in modes_to_blend:
        if lm not in inspiral_modes or lm not in merger_ringdown_modes:
            raise IOError(
                "We cannot blend {} mode as its missing in the input inspiral modes"
                " ({}) or the merger ringdown modes ({})".format(
                    lm, lm in inspiral_modes, lm in merger_ringdown_modes
                )
            )
    if verbose:
        print("Hybridizing the following modes: {}".format(modes_to_blend))
        print("By aligning {} mode".format(mode_to_align_by))
        print(
            "..and inheriting the phase/time shifts for alignment of {} modes".format(
                modes_not_aligned_by
            )
        )

    """ first we need to find the attachment region, based on the frequency """

    """ 
        We search left to right in merger-ringdown to avoid frequency fluctuations 
        after the merger, and right to left in inspiral to avoid frequency degeneracy
        caused by eccentricity 
    """
    el, em = mode_to_align_by

    if em <= 0:
        raise ValueError(
            f"mode_to_align_by must have positive m, got {mode_to_align_by}"
        )

    phase_mr_align = compute_phase(merger_ringdown_modes[(el, em)])
    frq_mr_align = compute_frequency(phase_mr_align, delta_t)

    t1_index_mr = find_first_value_location_in_series(
        frq_mr_align, frq_attach - frq_width / 2
    )
    t2_index_mr = find_first_value_location_in_series(
        frq_mr_align, frq_attach + frq_width / 2
    )

    """ 
    For eccentric inspiral, there will be multiple instances of the 
    same frequency. Pick the one having the highest index value (i.e. 
    the one at the rightmost occurrence in time) 

    """
    if blend_using_avg_orbital_frequency:
        t2_index_insp = find_last_value_location_in_series(
            inspiral_orbital_frequency, (frq_attach + frq_width / 2) / em
        )
        if verbose > 1:
            print(f"""Hybridizing using orbital frequency. Frequency
                {frq_attach + frq_width / 2}Hz found at {t2_index_insp}.
                """)
    else:
        phase_insp_align = compute_phase(inspiral_modes[(el, em)])
        frq_insp_align = compute_frequency(phase_insp_align, delta_t)
        t2_index_insp = find_last_value_location_in_series(
            frq_insp_align, frq_attach + frq_width / 2
        )

    # another way to define t2_index_mr is through number of points in the inspiral window
    t1_index_insp = t2_index_insp - (t2_index_mr - t1_index_mr)

    if verbose > 4:
        print(f"""In the merger-ringdown waveform, the hybridization frequency window
              [{frq_attach - frq_width / 2}, {frq_attach + frq_width / 2}]
              was found at indices [{t1_index_mr}, {t2_index_mr}]
              """)
        print(f"""In the inspiral waveform, the same window is to be found at
              indices [{t1_index_insp}, {t2_index_insp}.]
              """)
    """ 
        Theoretically, we NEED a timeshift to align the waveforms in frequency. 
        Instead of shifting one of the two waveforms for alignment, we are defining
        the time such that the frequencies are pre-aligned to the best of the 
        discrete interval errors. That is: 
            deltaT (timeshift) = t1_index_insp - t1_index_mr
        The mathematical way is to optimise the difference in frequencies over the matching 
        region and using that to determine deltaT, hence arriving at t1_index_mr. 
    """
    # For frequency computation below, we need at least one point before t1_index_insp/t1_index_mr
    # and at least one point after t2_index_insp/t2_index_mr to compute frequency at these endpoints
    # of the window using np.gradient with central differences. Making that sure here.
    if t1_index_insp < 1:
        raise ValueError("""Inspiral too short to accommodate the hybridization window. 
Try reducing frq_width or moving frq_attach to a higher value.""")
    if t2_index_insp + 2 > len(inspiral_modes[(el, em)]):
        raise ValueError(
            """Inspiral too short to accommodate the high attachment frequency. 
Try moving frq_attach to a lower value."""
        )

    if t1_index_mr < 1:
        raise ValueError(
            """Merger-ringdown too short to accommodate the hybridization window. 
Try reducing frq_width or moving frq_attach to a higher value."""
        )
    if t2_index_mr + 2 > len(merger_ringdown_modes[(el, em)]):
        raise ValueError(
            """Merger-ringdown too short to accommodate the high attachment frequency. 
Try moving frq_attach to a lower value."""
        )

    frq_insp_window = {}
    frq_mr_window = {}
    # t2_index_mr + 1 since the window is inclusive of t2_index_mr
    # Also, t2_index_mr + 1 should be <= len(merger_ringdown_modes[(el, em)]),
    # which is guaranteed because find_first_value_location_in_series
    # returns an index in the range [0, len(frq_mr_align)-1]
    frq_mr_window[(el, em)] = frq_mr_align[t1_index_mr : t2_index_mr + 1]

    if blend_using_avg_orbital_frequency:
        # We compute frequency only within the window (inclusive of endpoints).
        # Since it's computed by np.gradient which uses centered differences
        # to compute the phase derivative, we need to compute the phase at one
        # extra point on either side of the window to get accurate frequency at
        # the endpoints of the window.
        phase_insp_temp = compute_phase(
            inspiral_modes[(el, em)][t1_index_insp - 1 : t2_index_insp + 1 + 1]
        )
        frq_insp_window[(el, em)] = compute_frequency(phase_insp_temp, delta_t)[1:-1]
    else:
        frq_insp_window[(el, em)] = frq_insp_align[t1_index_insp : t2_index_insp + 1]

    # Non-alignment modes: window-only phase/frequency
    for el_i, em_i in modes_to_blend:
        if (el_i, em_i) == mode_to_align_by:
            continue
        phase_insp_temp = compute_phase(
            inspiral_modes[(el_i, em_i)][t1_index_insp - 1 : t2_index_insp + 1 + 1]
        )
        frq_insp_window[(el_i, em_i)] = compute_frequency(phase_insp_temp, delta_t)[
            1:-1
        ]

        phase_mr_temp = compute_phase(
            merger_ringdown_modes[(el_i, em_i)][t1_index_mr - 1 : t2_index_mr + 1 + 1]
        )
        frq_mr_window[(el_i, em_i)] = compute_frequency(phase_mr_temp, delta_t)[1:-1]

    amp_insp_window = {}
    amp_mr_window = {}

    for el_i, em_i in modes_to_blend:
        amp_insp_window[(el_i, em_i)] = compute_amplitude(
            inspiral_modes[(el_i, em_i)][t1_index_insp : t2_index_insp + 1]
        )
        amp_mr_window[(el_i, em_i)] = compute_amplitude(
            merger_ringdown_modes[(el_i, em_i)][t1_index_mr : t2_index_mr + 1]
        )

    """ Blending the mode amplitudes and frequencies using the blending function """

    amp_hyb_window = {}
    frq_hyb_window = {}
    phase_hyb_window = {}
    hybrid_modes = {}

    len_insp_window = t2_index_insp - t1_index_insp
    len_mr_window = t2_index_mr - t1_index_mr

    for el, em in modes_to_blend:
        amp_hyb_window[(el, em)] = blend_series(
            amp_insp_window[(el, em)],  # This should have length len_insp_window + 1
            amp_mr_window[(el, em)],  # This should have length len_mr_window + 1
            0,
            len_insp_window,
            0,
            len_mr_window,
        )
        frq_hyb_window[(el, em)] = blend_series(
            frq_insp_window[(el, em)],  # This should have length len_insp_window + 1
            frq_mr_window[(el, em)],  # This should have length len_mr_window + 1
            0,
            len_insp_window,
            0,
            len_mr_window,
        )
        """ Integrating frq_hyb to obtain phase_hyb. """

        phase_hyb_window[(el, em)] = (2 * np.pi) * cumulative_trapezoid(
            frq_hyb_window[(el, em)], dx=delta_t, initial=0
        )

    """ Right now the phase is integrated only inside the hybridization window, 
    need to add constants to preserve phase continuity and compile full IMR modes. """

    for el, em in modes_to_blend:
        inspiral_angle_at_window_start = -np.angle(
            inspiral_modes[(el, em)][t1_index_insp]
        )
        mr_angle_at_window_end = -np.angle(merger_ringdown_modes[(el, em)][t2_index_mr])

        if blend_aligning_merger_to_inspiral:
            # Because phase_hyb_window[(el, em)] starts from 0 phase
            # this will ensure that phase_hyb_window starts from
            # inspiral_angle_at_window_start
            phase_hyb_window[(el, em)] += inspiral_angle_at_window_start
            mode_within_window = amp_hyb_window[(el, em)] * np.exp(
                -1j * phase_hyb_window[(el, em)]
            )
            mr_phasor_shift = np.exp(
                -1j * (phase_hyb_window[(el, em)][-1] - mr_angle_at_window_end)
            )
            hybrid_modes[(el, em)] = np.concatenate(
                (
                    inspiral_modes[(el, em)][:t1_index_insp],
                    mode_within_window,
                    merger_ringdown_modes[(el, em)][t2_index_mr + 1 :]
                    * mr_phasor_shift,
                )
            )
        else:
            # Because phase_hyb_window[(el, em)] starts from 0 phase
            # this will ensure that phase_hyb_window ends at
            # mr_angle_at_window_end
            phase_hyb_window[(el, em)] += (
                mr_angle_at_window_end - phase_hyb_window[(el, em)][-1]
            )
            mode_within_window = amp_hyb_window[(el, em)] * np.exp(
                -1j * phase_hyb_window[(el, em)]
            )
            inspiral_phasor_shift = np.exp(
                -1j * (phase_hyb_window[(el, em)][0] - inspiral_angle_at_window_start)
            )
            hybrid_modes[(el, em)] = np.concatenate(
                (
                    inspiral_modes[(el, em)][:t1_index_insp] * inspiral_phasor_shift,
                    mode_within_window,
                    merger_ringdown_modes[(el, em)][t2_index_mr + 1 :],
                )
            )

    return (
        hybrid_modes,
        t1_index_insp,
        t1_index_mr,
        t2_index_insp,
        t2_index_mr,
        frq_insp_window,
        frq_hyb_window,
        frq_mr_window,
        amp_insp_window,
        amp_hyb_window,
        amp_mr_window,
    )
