# Copyright (C) 2026 Kaushik Paul, Akash Maurya
#
import numpy as np
from numba import njit
from scipy.signal import find_peaks
from pycbc.waveform.utils import td_taper

# AM: This code is basically the Python version of the Planck tapering C code in LAL:
# https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/_l_a_l_sim_inspiral_waveform_taper_8c_source.html
# The only new thing here is that while the LAL C code is restricted to only 2 extrema wide tapering, the
# following Python code can taper a user-specified number of extrema of the signal

LALSIMULATION_RINGING_EXTENT = 19


@njit
def Planck_window_LAL(data, taper_method, num_extrema_start=2, num_extrema_end=2):
    """
    Parameters:
    -----------
    data: 1D numpy array of reals
        data to taper
    taper_method: string
        Tapering method. Available methods are:
        "LAL_SIM_INSPIRAL_TAPER_START"
        "LAL_SIM_INSPIRAL_TAPER_END"
        "LAL_SIM_INSPIRAL_TAPER_STARTEND"
    num_extrema_start: int
        number of extrema till which to taper from the start
    num_extrema_end: int
        number of extrema till which to taper from the end

    Returns:
    --------
    window: 1D numpy array
        Planck tapering window
    """
    start = 0
    end = 0
    n = 0
    length = len(data)

    # Search for start and end of signal
    flag = 0
    i = 0
    while flag == 0 and i < length:
        if data[i] != 0.0:
            start = i
            flag = 1
        i += 1
    if flag == 0:
        raise ValueError("No signal found in the vector. Cannot taper.\n")

    flag = 0
    i = length - 1
    while flag == 0:
        if data[i] != 0.0:
            end = i
            flag = 1
        i -= 1

    # Check we have more than 2 data points
    if (end - start) <= 1:
        raise RuntimeError("Data less than 3 points, cannot taper!\n")

    # Calculate middle point in case of short waveform
    mid = int((start + end) / 2)

    window = np.ones(length)
    # If requested search for num_extrema_start-th peak from start and taper
    if taper_method != "LAL_SIM_INSPIRAL_TAPER_END":
        flag = 0
        i = start + 1
        while flag < num_extrema_start and i != mid:
            if abs(data[i]) >= abs(data[i - 1]) and abs(data[i]) >= abs(data[i + 1]):

                if abs(data[i]) == abs(data[i + 1]):
                    i += 1
                # only count local extrema more than 19 samples in
                if i - start > LALSIMULATION_RINGING_EXTENT:
                    flag += 1
                n = i - start
            i += 1

        # Have we reached the middle without finding `num_extrema_start` peaks?
        if flag < num_extrema_start:
            n = mid - start
            print(
                f"""WARNING: Reached the middle of waveform without finding {num_extrema_start} extrema.
Tapering only till the middle from the beginning."""
            )

        # Taper to that point
        realN = n
        window[: start + 1] = 0.0
        realI = np.arange(1, n - 1)
        z = (realN - 1.0) / realI + (realN - 1.0) / (realI - (realN - 1.0))
        window[start + 1 : start + n - 1] = 1.0 / (np.exp(z) + 1.0)

    # If requested search for num_extrema_end-th peak from end
    if (
        taper_method == "LAL_SIM_INSPIRAL_TAPER_END"
        or taper_method == "LAL_SIM_INSPIRAL_TAPER_STARTEND"
    ):
        i = end - 1
        flag = 0
        while flag < num_extrema_end and i != mid:
            if abs(data[i]) >= abs(data[i + 1]) and abs(data[i]) >= abs(data[i - 1]):
                if abs(data[i]) == abs(data[i - 1]):
                    i -= 1
                # only count local extrema more than 19 samples in
                if end - i > LALSIMULATION_RINGING_EXTENT:
                    flag += 1
                n = end - i
            i -= 1

        # Have we reached the middle without finding `num_extrema_end` peaks?
        if flag < num_extrema_end:
            n = end - mid
            print(
                f"""WARNING: Reached the middle of waveform without finding {num_extrema_end} extrema.
Tapering only till the middle from the end."""
            )

        # Taper to that point
        realN = n
        window[end:] = 0.0
        realI = -np.arange(-n + 2, 0)
        z = (realN - 1.0) / realI + (realN - 1.0) / (realI - (realN - 1.0))
        window[end - n + 2 : end] = 1.0 / (np.exp(z) + 1.0)

    return window


def compute_taper_width(waveform, method="cycles", fixed_duration=0.3, n_cycles=1):
    """
    Compute appropriate taper width for a gravitational waveform.

    Parameters:
    -----------
    waveform : TimeSeries
        The input waveform
    method : str
        'cycles': Based on number of GW cycles at start
        'fixed_time': Fixed time duration in seconds
    fixed_duration : float
        Fixed duration in seconds for 'fixed_time' method (default: 0.3)
    n_cycles : int
        Number of cycles for 'cycles' method (default: 1)

    Returns:
    --------
    float : Taper width in seconds
    """
    if method == "cycles":
        n_samples = min(5000, len(waveform))
        times = waveform.sample_times[:n_samples]
        data = np.abs(waveform.data[:n_samples])

        extrema, _ = find_peaks(data)

        # Check if the first point is an extremum
        if len(data) > 2:
            # Is first point a local maximum (peak)?
            if data[0] > data[1] and data[0] > data[2]:
                extrema = np.insert(extrema, 0, 0)
        else:
            raise ValueError("Waveform data length is less than three. Can't taper!\n")

        n_extrema_needed = 2 * n_cycles + 1

        if len(extrema) >= n_extrema_needed:
            taper_width = times[extrema[n_extrema_needed - 1]] - times[extrema[0]]

        else:
            # Fallback to fixed time
            taper_width = fixed_duration

    elif method == "fixed_time":
        # Simple fixed duration based on total waveform duration
        duration = waveform.sample_times[-1] - waveform.sample_times[0]
        taper_width = min(
            fixed_duration, duration * 0.1
        )  # 10% of waveform or fixed_duration

    else:
        raise ValueError(f"Unknown method: {method}. Use 'cycles' or 'fixed_time'")

    return taper_width


def apply_taper(
    waveform,
    beta=8,
    taper_width=None,
    method="cycles",
    fixed_duration=0.3,
    n_cycles=1,
    verbose=False,
):
    """
    Apply a time-domain taper to the start of the given waveform.

    Parameters:
    -----------
    waveform : TimeSeries
        The input waveform to be tapered
    beta : int
        Kaiser window parameter (default: 8)
    taper_width : float or None
        The width of the taper in seconds. If None, computed automatically
    method : str
        Method for auto-computing taper width ('cycles', 'fixed_time')
    fixed_duration : float
        Fixed duration for 'fixed_time' method (default: 0.3)
    n_cycles : int
        Number of cycles for 'cycles' method (default: 1)

    Returns:
    --------
    TimeSeries : The tapered waveform
    """

    t_start = waveform.sample_times[0]

    # Auto-compute taper width if not provided
    if taper_width is None:
        taper_width = compute_taper_width(
            waveform, method=method, fixed_duration=fixed_duration, n_cycles=n_cycles
        )
        if verbose:
            print(f"Auto-computed taper width: {taper_width:.3f} s (method: {method})")

    t_end_taper = t_start + taper_width

    return td_taper(waveform, t_start, t_end_taper, beta=beta, side="left")


def apply_taper_both_pols(hp, hc, beta=8, method="cycles", n_cycles=1, verbose=False):
    """Apply consistent taper to both polarizations based on hp."""
    taper_width = compute_taper_width(hp, method=method, n_cycles=n_cycles)
    if verbose:
        print(f"Taper width: {taper_width:.3f} s (computed from h+)")
    hp_tapered = apply_taper(
        hp,
        beta=beta,
        taper_width=taper_width,
        method="cycles",
        fixed_duration=0.3,
        n_cycles=1,
        verbose=verbose,
    )
    hc_tapered = apply_taper(
        hc,
        beta=beta,
        taper_width=taper_width,
        method="cycles",
        fixed_duration=0.3,
        n_cycles=1,
        verbose=verbose,
    )

    return (hp_tapered, hc_tapered, taper_width)
