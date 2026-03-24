def compute_taper_width(waveform, method='cycles', fixed_duration=0.3, n_cycles=1):
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
    import numpy as np
    from scipy.signal import find_peaks
    if method == 'cycles':
        n_samples = min(5000, len(waveform))
        times = waveform.sample_times[:n_samples]
        data = np.abs(waveform.data[:n_samples])
        
        extrema, _ = find_peaks(data)
            
        # Check if the first point is an extremum
        if len(data) > 2:
            # Is first point a local maximum (peak)?
            if (data[0] > data[1] and data[0] > data[2]):
                extrema = np.insert(extrema, 0, 0)
        else:
            raise ValueError("Waveform data length is less than three. Can't taper!\n")

        n_extrema_needed = 2 * n_cycles + 1
        
        if len(extrema) >= n_extrema_needed:
            taper_width = times[extrema[n_extrema_needed-1]] - times[extrema[0]]
            
        else:
            # Fallback to fixed time
            taper_width = fixed_duration
            
    elif method == 'fixed_time':
        # Simple fixed duration based on total waveform duration
        duration = waveform.sample_times[-1] - waveform.sample_times[0]
        taper_width = min(fixed_duration, duration * 0.1)  # 10% of waveform or fixed_duration
        
    else:
        raise ValueError(f"Unknown method: {method}. Use 'cycles' or 'fixed_time'")
    
    return taper_width


def apply_taper(waveform, beta=8, taper_width=None, method='cycles', fixed_duration=0.3, n_cycles=1,verbose=False):
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
    from pycbc.waveform.utils import td_taper
    
    t_start = waveform.sample_times[0]
    
    # Auto-compute taper width if not provided
    if taper_width is None:
        taper_width = compute_taper_width(waveform, method=method, fixed_duration=fixed_duration, n_cycles=n_cycles)
        if verbose:    
            print(f"Auto-computed taper width: {taper_width:.3f} s (method: {method})")
        
    t_end_taper = t_start + taper_width
    
    return td_taper(waveform, t_start, t_end_taper, beta=beta, side='left')

def apply_taper_both_pols(hp, hc, beta=8, method='cycles', n_cycles=1,verbose=False):
    """Apply consistent taper to both polarizations based on hp."""
    taper_width = compute_taper_width(hp, method=method, n_cycles=n_cycles)
    if verbose:    
        print(f"Taper width: {taper_width:.3f} s (computed from h+)")
    hp_tapered = apply_taper(hp, beta=beta, taper_width=taper_width, method='cycles', fixed_duration=0.3, n_cycles=1, verbose=verbose)
    hc_tapered = apply_taper(hc, beta=beta, taper_width=taper_width, method='cycles', fixed_duration=0.3, n_cycles=1, verbose=verbose)
    
    return (hp_tapered, hc_tapered, taper_width)