# Copyright (C) 2023 Kartikey Sharma, Prayush Kumar
#
"""Master function to hybridise any complex timeseries using the 'frequency' as a user input, 
specifically to be used for gravitational waveform hybridisation, fine-tuned for a single mode.
"""

import numpy as np
import scipy.optimize
from scipy.integrate import cumulative_trapezoid


def find_first_value_location_in_series(frq_timeseries, frq_desired):
    if frq_desired < np.min(frq_timeseries):
        raise Exception("Desired frequency out of bounds, lower than min frequency")

    if frq_desired > np.max(frq_timeseries):
        raise Exception("Desired frequency out of bounds, higher than max frequency")
    """ 
        We reverse the array and traverse it to find the location where the i_th value is more than
        the desired value while the i+1_th value is less, hence locating the desired value somewhere
        between those two points. We then choose the value closer to the value desired (among i and i+1) 
        and call it the location of the desired value. 
    """

    for idx, f_value in enumerate(frq_timeseries):
        if idx != len(frq_timeseries) - 1:
            if (
                frq_timeseries[idx] <= frq_desired
                and frq_timeseries[idx + 1] >= frq_desired
            ):
                fr1 = frq_timeseries[idx]
                fr2 = frq_timeseries[idx + 1]

                if abs(frq_desired - fr1) <= abs(frq_desired - fr2):
                    final_idx = idx
                else:
                    final_idx = idx + 1
                break
    return final_idx


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
    final_idx = len(reversed_freq_timeseries) - 1

    for idx, f_value in enumerate(reversed_freq_timeseries):
        if idx != len(reversed_freq_timeseries) - 1:
            if (
                reversed_freq_timeseries[idx] >= frq_desired
                and reversed_freq_timeseries[idx + 1] <= frq_desired
            ):
                fr1 = reversed_freq_timeseries[idx]
                fr2 = reversed_freq_timeseries[idx + 1]

                if abs(frq_desired - fr1) <= abs(frq_desired - fr2):
                    final_idx = idx
                else:
                    final_idx = idx + 1
                break
    return len(frq_timeseries) - 1 - final_idx


def mismatch_discrete(w1, w2, sample_indices_insp, sample_indices_mr):
    w1_d = w1[sample_indices_insp]
    w2_d = w2[sample_indices_mr]  # can't give the same comb to w2
    w2sq = np.square(np.abs(w2_d))
    # w1sq = np.square(np.abs(w2_d)) # another normalising factor can be (w1sq + w2sq) / 2
    diff = np.abs(w1_d - w2_d)
    diffsq = np.square(diff)
    mm = 0.5 * (np.sum(diffsq) / np.sum(w2sq))
    return mm


def align_in_phase(
    inspiral,
    merger_ringdown,
    sample_indices_insp,
    sample_indices_mr,
    t1_index_insp,
    t2_index_insp,
    t1_index_mr,
    t2_index_mr,
    m_mode=2,
    align_merger_to_inspiral=True,
):
    if len(inspiral) == 0:
        raise IOError(
            f"""You passed an inspiral waveform of zero length, to align
                      with merger-ringdown!"""
        )
    if (t2_index_insp + 1 - t1_index_insp) < 1:
        raise IOError(
            f"""You have passed a very narrow window for the inspiral
                      waveform's hybridization. As per your input, the inspiral
                      waveform from index {t1_index_insp} to {t2_index_insp + 1}
                      should be used to attach merger-ringdown"""
        )

    # Function alignes the two waveforms using the phase, optimised over the attachment region
    # m from l,m mode
    def optfn_ph(phaseshift_correction):
        if align_merger_to_inspiral:
            phase_corrected_inspiral = inspiral
            phase_corrected_merger_ringdown = merger_ringdown * np.exp(
                1j * m_mode * phaseshift_correction
            )
        else:
            phase_corrected_inspiral = inspiral * np.exp(
                1j * m_mode * phaseshift_correction
            )
            phase_corrected_merger_ringdown = merger_ringdown
        m_d = mismatch_discrete(
            phase_corrected_inspiral[t1_index_insp : t2_index_insp + 1],
            phase_corrected_merger_ringdown[t1_index_mr : t2_index_mr + 1],
            sample_indices_insp,
            sample_indices_mr,
        )
        return m_d

    phase_optimizer = scipy.optimize.minimize(optfn_ph, 0)
    phaseshift_required_for_alignment = phase_optimizer.x

    if align_merger_to_inspiral:
        aligned = merger_ringdown * np.exp(
            1j * m_mode * phaseshift_required_for_alignment
        )
        return (inspiral, aligned, phaseshift_required_for_alignment)
    else:
        aligned = inspiral * np.exp(1j * m_mode * phaseshift_required_for_alignment)
        return (aligned, merger_ringdown, phaseshift_required_for_alignment)


def blend_series(x1, x2, t1_index_insp, t2_index_insp, t1_index_mr, t2_index_mr):
    assert (
        t1_index_mr - t2_index_mr == t1_index_insp - t2_index_insp
    ), "Inconsistent indices passed to blending function"

    # blending fn is an array
    blfn_var = np.arange(t1_index_insp, t2_index_insp)
    tau = np.square(
        (
            np.sin(
                (np.pi / 2)
                * (blfn_var - t1_index_insp)
                / (t2_index_insp - t1_index_insp)
            )
        )
    )

    x_hyb = (1 - tau) * x1[t1_index_insp:t2_index_insp] + tau * x2[
        t1_index_mr:t2_index_mr
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
    no_sp=8,
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

    frq_attach: float
        Frequency (Hz) at which to align the inspiral and merger-ringdown modes.
    frq_width: {10.0, float}
        Frequency (Hz) window around the central attachment frequency over which
        hybridization of modes is performed.
    delta_t: {1/4096, float}
        Sample rate for timeseries (Hz)
    no_sp: {8, int}

    modes_to_blend: {[(2, 2), (3, 3), (4, 4)], list}
        List of modes as tuples of (l, m) values to blend
    mode_to_align_by: {(2, 2), tuple}
        One specific mode (l, m) value that is to be treated as baseline for
        time/phase alignment. We recommend using only the (2, 2) mode for this.
    blend_using_avg_orbital_frequency: (False, bool)
        Flag that enables use of orbit averaged orbital frequency for all
        calculations of inspiral to merger-ringdown attachment
    blend_aligning_merger_to_inspiral: (True, bool)
        Flag that ensures any phase difference between inspiral and merger-rd
        modes is removed by phase shifting the merger-rd portion. If set to
        False, the inspiral portion is phase shifted instead
    include_conjugate_modes: {True, bool}
        When set to True, we also consider (l, -m) modes in addition to (l, m) ones.
    verbose: {True, bool}
        Set this to True to enable logging output.
    """
    if frq_width <= 0:
        raise IOError(
            f"""You are trying to blend over a frequency window of
            negative length (= {frq_width}Hz). Fix this."""
        )
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

    # Get amplitude and phase for all modes
    phase_insp = {}
    frq_insp = {}
    phase_mr = {}
    frq_mr = {}

    for el, em in modes_to_blend:
        phase_insp[(el, em)] = compute_phase(inspiral_modes[(el, em)])
        frq_insp[(el, em)] = compute_frequency(phase_insp[(el, em)], delta_t)

        phase_mr[(el, em)] = compute_phase(merger_ringdown_modes[(el, em)])
        frq_mr[(el, em)] = compute_frequency(phase_mr[(el, em)], delta_t)

        if verbose:
            print(
                f"INSPIRAL mode ({el}, {em}) goes from {np.min(frq_insp[(el, em)])}Hz to"
                f" {np.max(frq_insp[(el, em)])}Hz"
            )
            print(
                f"MERGER mode ({el}, {em}) goes from {np.min(frq_mr[(el, em)])}Hz to"
                f" {np.max(frq_mr[(el, em)])}Hz"
            )

    """ first we need to find the attachment region, based on the frequency """

    """ 
        We search left to right in merger-ringdown to avoid frequency fluctuations 
        after the merger, and right to left in inspiral to avoid frequency degeneracy
        caused by eccentricity 
    """
    el, em = mode_to_align_by

    t1_index_mr = find_first_value_location_in_series(
        frq_mr[(el, em)], frq_attach - frq_width / 2
    )

    t2_index_mr = find_first_value_location_in_series(
        frq_mr[(el, em)], frq_attach + frq_width / 2
    )
    """ 
    For eccentric inspiral, there will be multiple instances of the 
    same frequency. Pick the one having the highest index value (i.e. 
    the one at the rightmost occurance in time) 

    """
    if blend_using_avg_orbital_frequency:
        t2_index_insp = find_last_value_location_in_series(
            inspiral_orbital_frequency, (frq_attach + frq_width / 2) / em
        )
        if verbose > 1:
            print(
                f"""Hybridizing using orbital frequency. Frequency
                  {frq_attach + frq_width / 2}Hz found at {t2_index_insp}.
                  The same frequency would have been found at index
                  {find_last_value_location_in_series(frq_insp[(el, em)],
                  frq_attach + frq_width / 2)} of mode frequency evolution.
                  """
            )
    else:
        t2_index_insp = find_last_value_location_in_series(
            frq_insp[(el, em)], frq_attach + frq_width / 2
        )

    # another way to define t2_index_mr is through number of points in the inspiral window
    t1_index_insp = t2_index_insp - (t2_index_mr - t1_index_mr)

    if verbose > 4:
        print(
            f"""In the merger-ringdown waveform, the hybridization frequency window
              [{frq_attach - frq_width / 2}, {frq_attach + frq_width / 2}]
              was found at indices [{t1_index_mr}, {t2_index_mr}]
              """
        )
        print(
            f"""In the inspiral waveform, the same window is to be found at
              indices [{t1_index_insp}, {t2_index_insp}.]
              """
        )
    """ 
        Theoretically, we NEED a timeshift to align the waveforms in frequency. 
        Instead of shifting one of the two waveforms for alignment, we are defining
        the time such that the frequencies are pre-aligned to the best of the 
        discrete interval errors. That is: 
            deltaT (timeshift) = t1_index_insp - t1_index_mr
        The mathematical way is to optimise the difference in frequencies over the matching 
        region and using that to determine deltaT, hence arriving at t1_index_mr. 
    """

    sample_indices_insp = (
        np.linspace(t1_index_insp, t2_index_insp, no_sp).astype(int) - t1_index_insp
    )
    sample_indices_mr = (
        sample_indices_insp  # since the attachment region in both has the same length
    )
    """ alignment using corrective phase addition """

    inspiral_modes_aligned = {}
    amp_insp_aligned = {}
    phase_insp_aligned = {}
    frq_insp_aligned = {}

    merger_ringdown_modes_aligned = {}
    amp_mr_aligned = {}
    phase_mr_aligned = {}
    frq_mr_aligned = {}

    (
        inspiral_modes_aligned[(el, em)],
        merger_ringdown_modes_aligned[(el, em)],
        phase_correction,
    ) = align_in_phase(
        inspiral_modes[(el, em)],
        merger_ringdown_modes[(el, em)],
        sample_indices_insp,
        sample_indices_mr,
        t1_index_insp,
        t2_index_insp,
        t1_index_mr,
        t2_index_mr,
        m_mode=em,
        align_merger_to_inspiral=blend_aligning_merger_to_inspiral,
    )

    amp_insp_aligned[(el, em)] = compute_amplitude(inspiral_modes_aligned[(el, em)])
    phase_insp_aligned[(el, em)] = compute_phase(inspiral_modes_aligned[(el, em)])

    amp_mr_aligned[(el, em)] = compute_amplitude(
        merger_ringdown_modes_aligned[(el, em)]
    )
    phase_mr_aligned[(el, em)] = compute_phase(merger_ringdown_modes_aligned[(el, em)])

    for el, em in modes_not_aligned_by:
        if blend_aligning_merger_to_inspiral:
            inspiral_modes_aligned[(el, em)] = inspiral_modes[(el, em)]
            merger_ringdown_modes_aligned[(el, em)] = merger_ringdown_modes[
                (el, em)
            ] * np.exp(1j * em * phase_correction)
        else:
            inspiral_modes_aligned[(el, em)] = inspiral_modes[(el, em)] * np.exp(
                1j * em * phase_correction
            )
            merger_ringdown_modes_aligned[(el, em)] = merger_ringdown_modes[(el, em)]
        amp_insp_aligned[(el, em)] = compute_amplitude(inspiral_modes_aligned[(el, em)])
        phase_insp_aligned[(el, em)] = compute_phase(inspiral_modes_aligned[(el, em)])
        amp_mr_aligned[(el, em)] = compute_amplitude(
            merger_ringdown_modes_aligned[(el, em)]
        )
        phase_mr_aligned[(el, em)] = compute_phase(
            merger_ringdown_modes_aligned[(el, em)]
        )

    """
        It would be same as frq_mr as the corrected phase factor will be canceled in the derivative, 
        defining frq_insp_aligned just for consistency 
    """
    frq_insp_aligned = frq_insp
    frq_mr_aligned = frq_mr

    """ Performing attachment using the blending function """

    amp_hyb_window = {}
    amp_hyb_full = {}
    frq_hyb_window = {}
    phase_hyb_window = {}
    phase_hyb_full = {}
    hybrid_modes = {}

    for el, em in modes_to_blend:
        amp_hyb_window[(el, em)] = blend_series(
            amp_insp_aligned[(el, em)],
            amp_mr_aligned[(el, em)],
            t1_index_insp,
            t2_index_insp,
            t1_index_mr,
            t2_index_mr,
        )
        frq_hyb_window[(el, em)] = blend_series(
            frq_insp_aligned[(el, em)],
            frq_mr_aligned[(el, em)],
            t1_index_insp,
            t2_index_insp,
            t1_index_mr,
            t2_index_mr,
        )
        """ Integrating frq_hyb to obtain phase_hyb and removing discontinuities, 
            compiling amp_hyb and phase_hyb to obtain the hybrid waveform.     """

        phase_hyb_window[(el, em)] = (2 * np.pi) * cumulative_trapezoid(
            frq_hyb_window[(el, em)], dx=delta_t, initial=0
        )

    """ Right now the phase is integrated only inside the hybrid window, 
    need to add constants to preserve phase continuity and compile full IMR phase """

    def remove_phase_discontinuity(phase_insp_, phase_hyb_window_, phase_mr_):
        if blend_aligning_merger_to_inspiral:
            delta1 = phase_insp_[t1_index_insp] - phase_hyb_window_[0]
            phase_hyb_1 = np.append(
                phase_insp_[:t1_index_insp], phase_hyb_window_ + delta1
            )
            delta2 = phase_hyb_1[t2_index_insp - 1] - phase_mr_[t2_index_mr - 1]
            phase_hyb_2 = np.append(
                phase_hyb_1[: t2_index_insp - 1], phase_mr_[t2_index_mr - 1 :] + delta2
            )
        else:
            delta1 = phase_hyb_window_[0] - phase_insp_[t1_index_insp]
            phase_hyb_1 = np.append(
                phase_insp_[:t1_index_insp] + delta1, phase_hyb_window_
            )
            delta2 = phase_mr_[t2_index_mr - 1] - phase_hyb_1[t2_index_insp - 1]
            phase_hyb_2 = np.append(
                phase_hyb_1[: t2_index_insp - 1] + delta2, phase_mr_[t2_index_mr - 1 :]
            )
        return phase_hyb_2

    for el, em in modes_to_blend:
        phase_hyb_full[(el, em)] = remove_phase_discontinuity(
            phase_insp_aligned[(el, em)],
            phase_hyb_window[(el, em)],
            phase_mr_aligned[(el, em)],
        )

        amp_hyb_full[(el, em)] = np.append(
            np.concatenate(
                [amp_insp_aligned[(el, em)][:t1_index_insp], amp_hyb_window[(el, em)]]
            )[: t2_index_insp - 1],
            amp_mr_aligned[(el, em)][t2_index_mr - 1 :],
        )

        hybrid_modes[(el, em)] = amp_hyb_full[(el, em)] * np.exp(
            -1j * phase_hyb_full[(el, em)]
        )

    return (
        hybrid_modes,
        t1_index_insp,
        t1_index_mr,
        t2_index_insp,
        t2_index_mr,
        frq_insp,
        frq_mr,
        frq_insp_aligned,
        frq_hyb_window,
        inspiral_modes_aligned,
        merger_ringdown_modes_aligned,
        sample_indices_insp,
        sample_indices_mr,
        amp_insp_aligned,
        amp_hyb_window,
        amp_mr_aligned,
        amp_hyb_full,
        phase_insp,
        phase_insp_aligned,
        phase_hyb_window,
        phase_mr_aligned,
        phase_hyb_full,
        phase_correction,
    )
