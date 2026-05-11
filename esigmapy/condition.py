# Copyright (C) 2026 Kaushik Paul, Akash Maurya
#
import numpy as np
import warnings
from numba import njit
from scipy.signal import find_peaks
from .utils import extract_waveform_info

# AM: This code is basically the Python version of the Planck tapering C code in LAL:
# https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/_l_a_l_sim_inspiral_waveform_taper_8c_source.html
# The only new thing here is that while the LAL C code is restricted to only 2 extrema wide tapering, the
# following Python code can taper a user-specified number of extrema of the signal

# In the LALSuite implementation, LALSIMULATION_RINGING_EXTENT = 19 was used.
LALSIMULATION_RINGING_EXTENT = 0


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
                # only count local extrema more than
                # LALSIMULATION_RINGING_EXTENT number of samples in
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
                # only count local extrema more than
                # LALSIMULATION_RINGING_EXTENT number of samples in
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


def compute_taper_width(
    waveform, method="cycles", fixed_duration=0.3, n_cycles=1, f_lower=1.0, delta_t=None
):
    """
    Compute appropriate taper width for a gravitational waveform.

    Parameters:
    -----------
    waveform : np.ndarray or TimeSeries
        The input waveform
    method : str
        'cycles': Based on number of GW cycles at start (default)
        'fixed_time': Fixed time duration in seconds
    fixed_duration : float
        Fixed duration in seconds for 'fixed_time' method (default: 0.3). If this is longer
        than 10% of the signal duration, we cap the window width to 10%.
    n_cycles : int
        Number of cycles for 'cycles' method (default: 1)
    f_lower : float
        Lowest frequency supported (Hz). Used to calculate max samples for one
        cycle (default: 1.0)
    delta_t : float or None
        Time step in seconds. Required for numpy array input. Automatically
        extracted for TimeSeries.

    Returns:
    --------
    int : Taper width in samples

    Raises:
    -------
    ValueError : If delta_t is not provided for numpy array or if method is invalid
    """
    try:
        info = extract_waveform_info(waveform, delta_t=delta_t)
    except (ValueError, TypeError) as e:
        raise ValueError(f"Failed to extract waveform info: {str(e)}")

    data = np.abs(info["data"])  # this ensures that minimas are also treated as peaks.
    delta_t = info["delta_t"]

    if len(data) < 3:
        raise ValueError(f"Waveform must have at least 3 samples. Got {len(data)}.")

    if method == "cycles":
        # Calculate max samples for one cycle at the lowest frequency,
        # this scales with f_lower provided.
        max_n_samples = int(1.0 / (f_lower * delta_t))
        n_samples = min(max_n_samples, len(data))

        data_subset = data[:n_samples]
        extrema, _ = find_peaks(data_subset)

        # Check if the first point is an extremum
        if len(data_subset) > 2:
            # Is first point a local maximum (peak)?
            if data_subset[0] > data_subset[1] and data_subset[0] > data_subset[2]:
                extrema = np.insert(extrema, 0, 0)

        n_extrema_needed = 2 * n_cycles + 1

        if len(extrema) >= n_extrema_needed:
            # Calculate taper width in indices, then convert to time
            taper_width = int(extrema[n_extrema_needed - 1] - extrema[0])
        else:
            # Fallback to fixed time
            warnings.warn(
                f"Not enough extrema found ({len(extrema)} < {n_extrema_needed}). "
                f"Falling back to fixed_duration={fixed_duration}s"
            )
            taper_width = int(fixed_duration / delta_t)

    elif method == "fixed_time":
        # Cap at 10% of the waveform length; take the smaller of that and
        # the user-specified fixed_duration (both expressed in samples)
        if int(fixed_duration / delta_t) > int(len(data) * 0.1):
            warnings.warn(
                "Requested tapering window width exceeds 10% of the waveform duration. "
                "Capping it to 10%.",
                UserWarning,
            )
        taper_width = min(int(fixed_duration / delta_t), int(len(data) * 0.1))
    else:
        raise ValueError(f"Unknown method: '{method}'. Use 'cycles' or 'fixed_time'")

    return taper_width


def apply_taper(
    waveform,
    taper_width=None,
    method="cycles",
    fixed_duration=0.3,
    n_cycles=1,
    f_lower=1.0,
    window="planck",
    beta_kaiser=8,
    delta_t=None,
    verbose=False,
):
    """
    Apply a time-domain taper to the start of the given waveform.

    Parameters:
    -----------
    waveform : TimeSeries or np.ndarray
        The input waveform to be tapered
    taper_width : float or None
        The width of the taper in samples. If None, computed automatically
    method : str
        Method for auto-computing taper width ('cycles' or 'fixed_time', default: 'cycles')
    fixed_duration : float
        Fixed duration for 'fixed_time' method (default: 0.3 s)
    n_cycles : int
        Number of cycles for 'cycles' method (default: 1)
    f_lower : float
        Lowest frequency supported (Hz). Used to calculate taper width for 'cycles' method (default: 1.0)
    window : str
        Window function to use: 'kaiser' (Kaiser window) or 'planck' (LAL Planck window)
        (default: 'planck')
    beta_kaiser : int
        Kaiser window parameter for kaiser window (default: 8)
    delta_t : float or None
        Time step. Required if waveform is a numpy array. Automatically extracted if waveform is a TimeSeries.
    verbose : bool
        Verbosity flag (default: False)

    Returns:
    --------
    TimeSeries or np.ndarray : The tapered waveform (same type as input)

    Raises:
    -------
    ValueError : If invalid window type or other parameter issues
    TypeError : If input waveform type is not recognized
    """
    # Extract waveform information
    try:
        info = extract_waveform_info(waveform, delta_t=delta_t)
    except (ValueError, TypeError) as e:
        raise TypeError(f"Invalid waveform input: {str(e)}")

    is_timeseries = info["is_timeseries"]
    data = info["data"].copy()
    delta_t = info["delta_t"]

    # Validate window choice
    if window not in ["kaiser", "planck"]:
        raise ValueError(f"Unknown window: '{window}'. Use 'kaiser' or 'planck'")

    # Auto-compute taper width if not provided
    if taper_width is None:
        try:
            taper_width = compute_taper_width(
                waveform,
                method=method,
                fixed_duration=fixed_duration,
                n_cycles=n_cycles,
                f_lower=f_lower,
                delta_t=delta_t,
            )
            if verbose:
                print(
                    f"Auto-computed taper width: {taper_width} samples"
                    f"({taper_width * delta_t:.6f} s) "
                    f"(method: {method}, window: {window})"
                )
        except Exception as e:
            raise ValueError(f"Failed to compute taper width: {str(e)}")
    if taper_width < 1:
        raise ValueError(
            f"Taper width ({taper_width} samples) must be at least 1 sample."
        )

    if window == "kaiser":
        try:
            from pycbc.waveform.utils import td_taper
            import pycbc.types as pt

            # Convert to TimeSeries, td_taper requires TimeSeries input.
            if not is_timeseries:
                temp_ts = pt.TimeSeries(data, delta_t=delta_t)
            else:
                temp_ts = waveform

            t_start = temp_ts.sample_times[0]
            t_end_taper = t_start + (taper_width * delta_t)
            tapered_ts = td_taper(
                temp_ts, t_start, t_end_taper, beta=beta_kaiser, side="left"
            )
            tapered_data = tapered_ts.data

            if verbose:
                print(
                    f"Applied kaiser window ({taper_width} samples, beta_kaiser={beta_kaiser})"
                )
        except ImportError as e:
            raise ImportError(f"PyCBC td_taper is required for kaiser window: {str(e)}")
        except Exception as e:
            raise RuntimeError(f"Error applying kaiser window: {str(e)}")

    elif window == "planck":
        try:
            planck_window = Planck_window_LAL(
                data, "LAL_SIM_INSPIRAL_TAPER_START", num_extrema_start=2
            )
            tapered_data = data * planck_window

            if verbose:
                print(f"Applied planck window ({taper_width} samples)")
        except Exception as e:
            raise RuntimeError(f"Error applying planck window: {str(e)}")

    # Return in same format as input
    if is_timeseries:
        try:
            import pycbc.types as pt

            return pt.TimeSeries(
                tapered_data, delta_t=delta_t, epoch=waveform.start_time
            )
        except Exception as e:
            raise RuntimeError(f"Error creating output TimeSeries: {str(e)}")
    else:
        return tapered_data


def apply_taper_both_pols(
    hp,
    hc,
    method="cycles",
    n_cycles=1,
    f_lower=1.0,
    window="planck",
    beta_kaiser=8,
    delta_t=None,
    verbose=False,
):
    """
    Apply consistent taper to both polarizations based on hp.

    Parameters:
    -----------
    hp : TimeSeries or np.ndarray
        Plus polarization waveform
    hc : TimeSeries or np.ndarray
        Cross polarization waveform
    method : str
        Taper width computation method: 'cycles' or 'fixed_time' (default: 'cycles')
    n_cycles : int
        Number of cycles for 'cycles' method (default: 1)
    f_lower : float
        Lowest frequency supported (Hz). Used to calculate taper width for 'cycles' method (default: 1.0)
    window : str
        Window function: 'kaiser' or 'planck' (default: 'planck')
    beta_kaiser : int
        Kaiser window parameter (default: 8)
    delta_t : float
        Sampling interval (default: None)
    verbose : bool
        Verbosity flag (default: False)

    Returns:
    --------
    tuple : (hp_tapered, hc_tapered, taper_width)
        Both polarizations tapered with same taper_width and window

    Raises:
    -------
    TypeError : If hp and hc are not the same type or incompatible
    ValueError : If inputs are invalid
    """
    # Validate that hp and hc are compatible types
    hp_is_array = isinstance(hp, np.ndarray)
    hc_is_array = isinstance(hc, np.ndarray)

    if hp_is_array != hc_is_array:
        raise TypeError(
            "hp and hc must be the same type (both numpy array or both TimeSeries)"
        )

    # Extract info from both
    try:
        hp_info = extract_waveform_info(hp, delta_t=delta_t)
        hc_info = extract_waveform_info(hc, delta_t=delta_t)
    except (ValueError, TypeError) as e:
        raise TypeError(f"Failed to extract polarization info: {str(e)}")

    if hp_info["delta_t"] != hc_info["delta_t"]:
        raise ValueError(
            f"hp and hc have different delta_t: {hp_info['delta_t']} vs {hc_info['delta_t']}"
        )

    if len(hp_info["data"]) != len(hc_info["data"]):
        raise ValueError(
            f"hp and hc have different lengths: {len(hp_info['data'])} vs {len(hc_info['data'])}"
        )

    # Compute taper width from hp
    try:
        taper_width = compute_taper_width(
            hp,
            method=method,
            n_cycles=n_cycles,
            f_lower=f_lower,
            delta_t=hp_info["delta_t"],
        )
        if verbose:
            print(
                f"Computed taper width from h+: {taper_width} samples "
                f"({taper_width * hp_info['delta_t']:.6f} s) "
                f"(method: {method}, n_cycles: {n_cycles})"
            )
    except Exception as e:
        raise ValueError(f"Failed to compute taper width: {str(e)}")

    # Apply same taper to both polarizations
    try:
        hp_tapered = apply_taper(
            hp,
            taper_width=taper_width,
            method=method,
            n_cycles=n_cycles,
            f_lower=f_lower,
            window=window,
            beta_kaiser=beta_kaiser,
            delta_t=hp_info["delta_t"],
            verbose=verbose,
        )
        hc_tapered = apply_taper(
            hc,
            taper_width=taper_width,
            method=method,
            n_cycles=n_cycles,
            f_lower=f_lower,
            window=window,
            beta_kaiser=beta_kaiser,
            delta_t=hc_info["delta_t"],
            verbose=verbose,
        )
    except Exception as e:
        raise RuntimeError(f"Error applying taper to polarizations: {str(e)}")

    if verbose:
        print(f"Tapered both polarizations with {window} window")

    return (hp_tapered, hc_tapered, taper_width)
