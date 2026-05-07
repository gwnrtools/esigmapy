# Copyright (C) 2026 Akash Maurya, Prayush Kumar
#
import time
from lal import MSUN_SI, PC_SI
import lalsimulation as ls

# Set of approximants that we support for
# the plunge-merger-ringdown piece of ESIGMA.
LALSIM_APPROXIMANTS = ["NRSur7dq4", "SEOBNRv4PHM"]
PYSEOBNR_APPROXIMANTS = ["SEOBNRv5HM", "SEOBNRv5PHM"]
SUPPORTED_MR_APPROXIMANTS = PYSEOBNR_APPROXIMANTS + LALSIM_APPROXIMANTS


def check_available_mr_approximants(approximant):
    if approximant not in SUPPORTED_MR_APPROXIMANTS:
        raise IOError(
            f"We cannot generate individual modes for {approximant}. "
            f"Try one of: {SUPPORTED_MR_APPROXIMANTS}."
        )


def _generate_lalsim_modes(
    mass1,
    mass2,
    f_lower,
    delta_t,
    spin1z,
    spin2z,
    coa_phase,
    distance,
    f_ref,
    modes_to_use,
    approximant,
):
    try:
        approximant_func = getattr(ls, approximant)
    except AttributeError:
        raise IOError(
            f"{approximant} is not available in your lalsimulation installation."
        )

    hlm_mr = ls.SimInspiralChooseTDModes(
        coa_phase,  # phiRef
        delta_t,  # deltaT
        mass1 * MSUN_SI,
        mass2 * MSUN_SI,
        0,  # spin1x
        0,  # spin1y
        spin1z,
        0,  # spin2x
        0,  # spin2y
        spin2z,
        f_lower,  # f_min
        f_ref,  # f_ref
        distance * PC_SI * 1.0e6,
        None,  # LALpars
        4,  # lmax
        approximant_func,
    )

    # Extracting only the modes we need
    modes_mr_numpy = {}
    while hlm_mr is not None:
        key = (hlm_mr.l, hlm_mr.m)
        if key in modes_to_use:
            modes_mr_numpy[key] = hlm_mr.mode.data.data
        hlm_mr = hlm_mr.next

    return modes_mr_numpy


def _generate_pyseobnr_modes(
    mass1,
    mass2,
    f_lower,
    delta_t,
    spin1z,
    spin2z,
    coa_phase,
    distance,
    f_ref,
    modes_to_use,
    approximant,
):
    try:
        from pyseobnr.generate_waveform import GenerateWaveform
    except ImportError:
        raise IOError(f"{approximant} requires pyseobnr, which is not installed.")

    # Converting to all +ve m modes to be compatible with pyseobnr.
    # pyseobnr outputs all the modes, including the negative m modes.
    mode_array = list({(l, abs(m)) for l, m in modes_to_use})
    wfm_gen = GenerateWaveform(
        dict(
            mass1=mass1,  # in solar masses
            mass2=mass2,
            spin1z=spin1z,
            spin2z=spin2z,
            distance=distance,  # in Mpc
            phi_ref=coa_phase,
            deltaT=delta_t,
            f22_start=f_lower,
            f_ref=f_ref,
            mode_array=mode_array,
            approximant=approximant,
        )
    )
    # Generate mode dictionary
    _, hlm = wfm_gen.generate_td_modes()

    # Extracting only the modes we need
    modes_mr_numpy = {key: hlm[key] for key in modes_to_use}

    return modes_mr_numpy


def get_mr_modes(
    mass1,
    mass2,
    f_lower,
    delta_t,
    spin1z=0.0,
    spin2z=0.0,
    coa_phase=None,
    distance=1.0,
    f_ref=None,
    modes_to_use=[(2, 2), (3, 3), (4, 4)],
    approximant="NRSur7dq4",
    verbose=False,
):
    """
    Returns quasi-circular IMR GW modes to be used as the
    plunge-merger-ringdown (PMR) piece of ESIGMA.
    Available approximants are:
    NRSur7dq4, SEOBNRv4PHM, SEOBNRv5HM, SEOBNRv5PHM

    Parameters:
    -----------
        mass1, mass2    -- Binary's component masses (in solar masses)
        f_lower         -- Starting frequency of the waveform (in Hz)
        f_ref           -- Reference frequency at which to define the
                           waveform parameters.
                           None by default, which means f_ref = f_lower.
        delta_t         -- Waveform's time grid-spacing (in s)
        spin1z, spin2z  -- z-components of component dimensionless
                           spins (lies in (-1,1))
        coa_phase       -- Coalescence phase of the binary (in rad)
        distance        -- Luminosity distance to the binary (in Mpc)
        modes_to_use    -- GW modes to use. List of tuples (l, |m|)
        approximant     -- Choose the plunge-merger-ringdown model.
                           Available choices:
                           NRSur7dq4, SEOBNRv4PHM  (requires `lalsimulation`)
                           SEOBNRv5HM, SEOBNRv5PHM (requires `pyseobnr`)
        verbose         -- Verbosity level
    Returns:
    --------
        modes_mr_numpy  -- Dictionary of PMR GW modes (NumPy arrays)
    """
    if coa_phase is None:
        coa_phase = 0.0
    if f_ref is None:
        f_ref = f_lower

    itime = time.perf_counter()
    if approximant in LALSIM_APPROXIMANTS:
        modes_mr_numpy = _generate_lalsim_modes(
            mass1=mass1,
            mass2=mass2,
            f_lower=f_lower,
            delta_t=delta_t,
            spin1z=spin1z,
            spin2z=spin2z,
            coa_phase=coa_phase,
            distance=distance,
            f_ref=f_ref,
            modes_to_use=modes_to_use,
            approximant=approximant,
        )
        if verbose:
            print(f"Used {approximant} as PMR via lalsimulation.")

    elif approximant in PYSEOBNR_APPROXIMANTS:
        modes_mr_numpy = _generate_pyseobnr_modes(
            mass1=mass1,
            mass2=mass2,
            f_lower=f_lower,
            delta_t=delta_t,
            spin1z=spin1z,
            spin2z=spin2z,
            coa_phase=coa_phase,
            distance=distance,
            f_ref=f_ref,
            modes_to_use=modes_to_use,
            approximant=approximant,
        )
        if verbose:
            print(f"Used {approximant} as PMR via pyseobnr.")

    else:
        raise ValueError(
            f"""Invalid choice of approximant for plunge-merger-ringdown: {approximant}.
Available choices are: {SUPPORTED_MR_APPROXIMANTS}."""
        )

    if verbose:
        print(f"MR modes generation took: {time.perf_counter() - itime} seconds.")
    return modes_mr_numpy
