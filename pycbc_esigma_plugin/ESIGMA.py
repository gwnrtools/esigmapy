from pycbc.types.timeseries import TimeSeries
from copy import deepcopy

try:
    # some versions of pycbc include td_taper in pycbc.waveform
    from pycbc.waveform import td_taper
except:
    # some other versions of pycbc include td_taper in pycbc.waveform.utils
    from pycbc.waveform.utils import td_taper

def taper_signal(signal, beta=5):
    """
    Returns tapered signal between start and start+0.4s.

    Parameters
    ----------

    signal : PyCBC TimeSeries
    beta : The beta parameter to use for the Kaiser window. See scipy.signal.kaiser for details. Default is 5.

    Returns
    -------

    signal_tapered : PyCBC TimeSeries
    """
    signal_length = signal.sample_times[-1] - signal.sample_times[0]
    if signal_length <= 0.4:
        taper_window = signal_length/7
    else:
        taper_window = 0.4
    signal_tapered = td_taper(signal, signal.sample_times[0], signal.sample_times[0]+taper_window, beta=beta)
    return(signal_tapered)

def use_modified_input_params(**input_params):
    # sim_inspiral table format uses alpha, alpha1, alpha2.. for additional parameters
    # hence using alpha as a proxy for eccentricity

    modified_input_params = deepcopy(input_params)
    verbose = modified_input_params.get("verbose", False)
    
    if 'alpha' in input_params:
        eccentricity = float(input_params.get("alpha", 0))
        modified_input_params["eccentricity"] = eccentricity
        if verbose:
                print(f"Using eccentricity from `alpha` column, value = {eccentricity}")
    if 'alpha1' in input_params:
        mean_anomaly = float(input_params.get("alpha1", 0))
        modified_input_params["mean_anomaly"] = mean_anomaly
        if verbose:
                print(f"Using mean_anomaly from `alpha1` column, value = {mean_anomaly}")
    return(modified_input_params)

def IMRESIGMAHM_td(**input_params):
    """
    Returns tapered time domain gravitational polarizations for IMRESIGMAHM waveform containing all (l,|m|) modes available.

    Parameters
    ----------

    Takes the same parameters as pycbc.waveform.get_td_waveform().
    
    Returns
    -------

    hplus : PyCBC TimeSeries
        The plus-polarization of the waveform in time domain tapered from start to 0.4s.
    hcross : PyCBC TimeSeries
        The cross-polarization of the waveform in time domain tapered from start to 0.4s.
    """
    #importing here instead of globally to avoid circular imports
    from gwnr.waveform import esigma_utils

    wf_input_params = use_modified_input_params(**input_params)
    
    hp, hc = esigma_utils.get_imr_esigma_waveform(**wf_input_params)
    hp_ts = TimeSeries(hp, input_params['delta_t'])
    hc_ts = TimeSeries(hc, input_params['delta_t'])
    
    hp_tapered = taper_signal(hp_ts)
    hc_tapered = taper_signal(hc_ts)
    return(hp_tapered, hc_tapered)

def IMRESIGMA_td(**input_params):
    """
    Returns tapered time domain gravitational polarizations for IMRESIGMA waveform containing only the (l,|m|) = (2,2) mode.

    Parameters
    ----------

    Takes the same parameters as pycbc.waveform.get_td_waveform().
    
    Returns
    -------

    hplus : PyCBC TimeSeries
        The plus-polarization of the waveform in time domain tapered from start to 0.4s.
    hcross : PyCBC TimeSeries
        The cross-polarization of the waveform in time domain tapered from start to 0.4s.
    """
    #importing here instead of globally to avoid circular imports
    from gwnr.waveform import esigma_utils

    wf_input_params = use_modified_input_params(**input_params)
    
    hp, hc = esigma_utils.get_imr_esigma_waveform(**wf_input_params, modes_to_use=[(2, 2)])
    hp_ts = TimeSeries(hp, input_params['delta_t'])
    hc_ts = TimeSeries(hc, input_params['delta_t'])
    
    hp_tapered = taper_signal(hp_ts)
    hc_tapered = taper_signal(hc_ts)
    return(hp_tapered, hc_tapered)

def InspiralESIGMAHM_td(**input_params):
    """
    Returns tapered time domain gravitational polarizations for InspiralESIGMAHM waveform containing all (l,|m|) modes available.

    Parameters
    ----------

    Takes the same parameters as pycbc.waveform.get_td_waveform().
    
    Returns
    -------

    hplus : PyCBC TimeSeries
        The plus-polarization of the waveform in time domain tapered from start to 0.4s.
    hcross : PyCBC TimeSeries
        The cross-polarization of the waveform in time domain tapered from start to 0.4s.
    """
    #importing here instead of globally to avoid circular imports
    from gwnr.waveform import esigma_utils

    wf_input_params = use_modified_input_params(**input_params)
    
    _, hp, hc = esigma_utils.get_inspiral_esigma_waveform(**wf_input_params)
    hp_ts = TimeSeries(hp, input_params['delta_t'])
    hc_ts = TimeSeries(hc, input_params['delta_t'])
    
    hp_tapered = taper_signal(hp_ts)
    hc_tapered = taper_signal(hc_ts)
    return(hp_tapered, hc_tapered)

def InspiralESIGMA_td(**input_params):
    """
    Returns tapered time domain gravitational polarizations for InspiralESIGMA waveform containing only the (l,|m|) = (2,2) mode.

    Parameters
    ----------

    Takes the same parameters as pycbc.waveform.get_td_waveform().
    
    Returns
    -------

    hplus : PyCBC TimeSeries
        The plus-polarization of the waveform in time domain tapered from start to 0.4s.
    hcross : PyCBC TimeSeries
        The cross-polarization of the waveform in time domain tapered from start to 0.4s.
    """
    #importing here instead of globally to avoid circular imports
    from gwnr.waveform import esigma_utils

    wf_input_params = use_modified_input_params(**input_params)
    
    _, hp, hc = esigma_utils.get_inspiral_esigma_waveform(**wf_input_params, modes_to_use=[(2, 2)])
    hp_ts = TimeSeries(hp, input_params['delta_t'])
    hc_ts = TimeSeries(hc, input_params['delta_t'])
    
    hp_tapered = taper_signal(hp_ts)
    hc_tapered = taper_signal(hc_ts)
    return(hp_tapered, hc_tapered)
