# Copyright (C) 2026 Akash Maurya, Prayush Kumar

"""Functions for generating ESIGMASur waveforms"""

from __future__ import absolute_import

import numpy as np
import time
import esigmapy
import lal
import lalsimulation as ls
import pycbc.types as pt
from esigmapy.utils import f_ISCO_spin
from esigmapy.generator import _get_transition_frequency_window
from .surrogate import _get_surrogate

ECCENTRICITY_LEVEL_ISCO_WARNING = 0.02
ECCENTRICITY_LEVEL_ISCO_ERROR = 0.1


def get_surrogate_object():
    """
    Returns the surrogate object. Useful for advanced users who want to
    directly use the base surrogate class' functionalities.
    """
    return _get_surrogate()


def get_inspiral_esigmasur_modes(
    mass1,
    mass2,
    reference_eccentricity=0.0,
    reference_mean_anomaly=0.0,
    delta_t=None,
    times=None,
    t_start=None,
    t_end=None,
    distance=1.0,
    include_conjugate_modes=False,
    return_orbital_params=False,
    return_pycbc_timeseries=True,
    verbose=False,
):
    """
    Returns inspiral ESIGMASur GW modes

    Parameters:
    -----------
        mass1, mass2            -- Binary's component masses (in solar masses)
        delta_t                 -- Waveform's time grid-spacing (in s)
        reference_eccentricity  -- Eccentricity at reference time of surrogate
        reference_mean_anomaly  -- Mean anomaly at reference time of surrogate (in rad)
        distance                -- Luminosity distance to the binary (in Mpc)
        t_start, t_end          -- Start and end times of the waveform
                                   to be generated (in seconds).
                                   Note that the surrogate defines t=0 at the end of
                                   inspiral, so t_start and t_end should be negative
                                   and t_start < t_end <= 0.
                                   Defaults to the full duration of the surrogate.
        include_conjugate_modes -- If True, negative "m" modes are included as well
        return_orbital_params   -- If True, returns the orbital evolution of all the orbital elements.
                                   Can also be a list of orbital variable names to
                                   return only those specific variables. Available
                                   orbital variables are:
                                   ['x', 'e', 'l']
        return_pycbc_timeseries -- If True, returns data in the form of PyCBC timeseries.
                                   True by default.
        verbose                 -- Verbosity flag

    Returns:
    --------
        t                 -- Time grid (in seconds).
                             Returned only if return_pycbc_timeseries=False
        orbital_var_dict  -- Dictionary of evolution of orbital elements.
                             Returned only if "return_orbital_params" is specified
        modes             -- Dictionary of GW modes
    """

    if return_orbital_params:
        orbital_var_names = ["x", "e", "l"]
        if return_orbital_params != True:
            for name in return_orbital_params:
                if name not in orbital_var_names:
                    raise Exception(
                        f"{name} is not a valid orbital variable name. Available orbital variable names are: {orbital_var_names}."
                    )

    sur = _get_surrogate()
    total_mass = mass1 + mass2
    q = mass1 / mass2
    if q < 1:
        q = 1 / q

    modes_to_use = [(2, 2)]

    # Calculating the orbital variables
    itime = time.perf_counter()
    modes = {}

    retval = sur(
        M=total_mass,
        params=[
            q,
            reference_eccentricity,
            reference_mean_anomaly,
        ],  # q, e, l at t=t_ref
        delta_t=delta_t,
        t_start=t_start,
        t_end=t_end,
        times=times,
        remove_start_phase=True,
        return_orbital_variables=return_orbital_params,
    )
    # Currently, the surrogate only supports (2, 2) mode, so we directly retrieve that from the returned dictionary.
    el, em = modes_to_use[0]
    if return_orbital_params:
        t, orb_vars, modes[(el, em)] = retval
    else:
        t, modes[(el, em)] = retval
    modes[(el, em)] /= distance

    if include_conjugate_modes:
        for el, em in modes_to_use:
            modes[(el, -em)] = (-1) ** el * np.conjugate(modes[(el, em)])

    if return_pycbc_timeseries:
        if times is None:
            modes = {
                k: pt.TimeSeries(
                    modes[k],
                    delta_t=delta_t,
                    epoch=-delta_t * (len(modes[k]) - 1),
                )
                for k in modes
            }
        else:
            raise ValueError(
                """Cannot return PyCBC TimeSeries when the user provides custom time grid via `times` due to the possibilty of it being a non-uniform time-grid. 
Please set `return_pycbc_timeseries=False` if you want to provide custom time grid."""
            )
    if verbose:
        print(f"Modes generation took: {time.perf_counter() - itime} seconds")

    if return_orbital_params:
        orbital_var_dict = {}
        if return_orbital_params == True:
            return_orbital_params = orbital_var_names

        if return_pycbc_timeseries:
            # No need to check and raise error here if `times` is provided,
            # because the error will be raised while trying to convert modes
            # to PyCBC TimeSeries above.
            for name in return_orbital_params:
                exec(
                    f"orbital_var_dict['{name}'] = pt.TimeSeries(orb_vars['{name}'], delta_t=delta_t, epoch=-delta_t * (len(orb_vars['{name}'])-1))"
                )
            return orbital_var_dict, modes

        for name in return_orbital_params:
            exec(f"orbital_var_dict['{name}'] = orb_vars['{name}']")
        # return (t - t[-1]), orbital_var_dict, modes
        return t, orbital_var_dict, modes

    if return_pycbc_timeseries:
        return modes
    # return (t - t[-1]), modes
    return t, modes


def get_inspiral_esigmasur_waveform(
    mass1,
    mass2,
    reference_eccentricity=0.0,
    reference_mean_anomaly=0.0,
    delta_t=None,
    t_start=None,
    t_end=None,
    times=None,
    inclination=0.0,
    coa_phase=0.0,
    distance=1.0,
    return_orbital_params=False,
    return_pycbc_timeseries=True,
    verbose=False,
    **kwargs,
):
    """
    Returns inspiral ESIGMASur GW polarizations

    Parameters:
    -----------
        mass1, mass2            -- Binary's component masses (in solar masses)
        reference_eccentricity  -- Eccentricity at reference time of surrogate
        reference_mean_anomaly  -- Mean anomaly at reference time of surrogate (in rad)
        delta_t                 -- Waveform's time grid-spacing (in s).
                                   Can be omitted if providing custom time grid via
                                   `times` argument.
        t_start, t_end          -- Start and end times of the waveform
                                   to be generated (in seconds).
                                   Note that the surrogate defines t=0 at the end of
                                   inspiral, so t_start and t_end should be negative
                                   and t_start < t_end <= 0.
                                   Defaults to the full duration of the surrogate.
        times                   -- Custom time grid (can be non-uniform) on which the
                                   waveform should be generated. Should be a numpy
                                   array of time values in seconds.
                                   If provided, `delta_t`, `t_start` and `t_end` are ignored.
                                   Also set `return_pycbc_timeseries=False` to use this option.
        inclination             -- Inclination (in rad), defined as the angle between
                                   the orbital angular momentum L and the line-of-sight
        coa_phase               -- Coalesence phase of the binary (in rad)
        distance                -- Luminosity distance to the binary (in Mpc)
        return_orbital_params   -- If True, returns the orbital evolution of all the
                                   orbital elements (in geometrized units). Can also be
                                   a list of orbital variable names to return only
                                   those specific variables. Available orbital
                                   variables names are:
                                   ['x', 'e', 'l']
        return_pycbc_timeseries -- If True, returns data in the form of PyCBC timeseries.
                                    True by default. Set to False if you want to provide custom time grid via `times` argument, or if you want the output in numpy arrays.
        verbose                 -- Verbosity level.
                                   Available values are: 0, 1, 2

    Returns:
    --------
        t                 -- Time grid (in seconds).
                             Returned only if return_pycbc_timeseries=False
        orbital_var_dict  -- Dictionary of evolution of orbital elements.
                             Returned only if "return_orbital_params" is specified
        hp, hc            -- Plus and cross GW polarizations
    """

    retval = get_inspiral_esigmasur_modes(
        mass1=mass1,
        mass2=mass2,
        reference_eccentricity=reference_eccentricity,
        reference_mean_anomaly=reference_mean_anomaly,
        t_start=t_start,
        t_end=t_end,
        delta_t=delta_t,
        times=times,
        distance=distance,
        include_conjugate_modes=True,  # Always include conjugate modes while generating polarizations
        return_orbital_params=return_orbital_params,
        verbose=verbose,
        return_pycbc_timeseries=False,
    )

    if return_orbital_params:
        t, orbital_var_dict, modes_inspiral = retval
    else:
        t, modes_inspiral = retval

    hp, hc = esigmapy.utils.get_polarizations_from_multipoles(
        modes_inspiral,
        inclination=inclination,
        coa_phase=np.pi / 2 - coa_phase,
        verbose=verbose,
    )

    if return_pycbc_timeseries:
        if times is None:
            hp = pt.TimeSeries(hp, delta_t=delta_t, epoch=-delta_t * (len(hp) - 1))
            hc = pt.TimeSeries(hc, delta_t=delta_t, epoch=-delta_t * (len(hc) - 1))
        else:
            raise ValueError(
                """Cannot return PyCBC TimeSeries when the user provides custom time grid via `times` due to the possibilty of it being a non-uniform time-grid. 
Please set `return_pycbc_timeseries=False` if you want to provide custom time grid."""
            )

    if return_orbital_params:
        if return_pycbc_timeseries:
            for name in orbital_var_dict:
                exec(
                    f"orbital_var_dict['{name}'] = pt.TimeSeries(orbital_var_dict['{name}'], delta_t=delta_t, epoch=-delta_t * (len(orbital_var_dict['{name}'])-1))"
                )
            return orbital_var_dict, hp, hc
        return t, orbital_var_dict, hp, hc

    if return_pycbc_timeseries:
        return hp, hc
    return t, hp, hc


def get_imr_esigmasur_mode(
    mass1,
    mass2,
    delta_t,
    reference_eccentricity=0.0,
    reference_mean_anomaly=0.0,
    t_start=None,
    distance=1.0,
    coa_phase=None,
    include_conjugate_modes=False,
    f_mr_transition=None,
    f_window_mr_transition=None,
    num_hyb_orbits=0.25,
    blend_aligning_merger_to_inspiral=True,
    keep_f_mr_transition_at_center=False,
    merger_ringdown_approximant="NRSur7dq4",
    return_hybridization_info=False,
    return_orbital_params=False,
    failsafe=True,
    verbose=False,
):
    """
    Returns IMR GW modes constructed using ESIGMASur for inspiral and
    NRSur7dq4/SEOBNRv4PHM for merger-ringdown

    Parameters:
    -----------
        mass1, mass2              -- Binary's component masses (in solar masses)
        delta_t                   -- Waveform's time grid-spacing (in s)
        reference_eccentricity    -- Eccentricity at reference time of surrogate
        reference_mean_anomaly    -- Mean anomaly at reference time of surrogate (in rad)
        t_start                   -- Start time of the waveform to be generated (in seconds).
                                     Note that the surrogate defines t=0 at the end of
                                     inspiral, so t_start should be negative
                                     and t_start < 0.
                                     Defaults to the full duration of the surrogate.
        coa_phase                 -- Coalesence phase of the binary (in rad)
        distance                  -- Luminosity distance to the binary (in Mpc)
        include_conjugate_modes   -- If True, (l, -|m|) modes are included as
                                     well
        f_mr_transition           -- Inspiral to merger transition GW frequency
                                     (Hz). Should correspond to the mode
                                     using which alignment is done, i.e. (2,2) mode here.
                                     Defaults to the minimum of the Kerr and
                                     Schwarzschild ISCO frequency equivalent
                                     for the (2,2)-mode.
        f_window_mr_transition    -- Hybridization frequency window (in Hz).
                                     Should correspond to the mode
                                     using which alignment is done, i.e. (2,2) mode here.
                                     Disabled by the default value (None). In
                                     such a case, the hybridization proceeds
                                     over a window of `num_hyb_orbits` orbital
                                     cycles (1 orbital cycle ~ 2 GW cycles)
                                     that ends at the frequency value given by
                                     `f_mr_transition`.
                                     Also see `keep_f_mr_transition_at_center`
                                     to choose the position of `f_mr_transition`
                                     within this window.
        num_hyb_orbits            -- number of orbital cycles to blend over.
                                     Only used if f_window_mr_transition is not
                                     specified.
        blend_aligning_merger_to_inspiral -- (default: False) If True, the
                                     merger-ringdown mode would be phase aligned
                                     to the inspiral
                                     If False, the inspiral is phase aligned
                                     Note: specify the desired
        keep_f_mr_transition_at_center -- If True, `f_mr_transition` is kept at
                                     the center of the hybridization window.
                                     Otherwise, it's kept at the end of the
                                     window (default).
        merger_ringdown_approximant    -- Choose merger-ringdown model. Tested
                                     choices: [NRSur7dq4, SEOBNRv4PHM]
        return_hybridization_info -- If True, returns hybridization related data
        return_orbital_params     -- If True, returns the orbital evolution of
                                     all the orbital elements (in
                                     geometrized units). Can also be a list of
                                     orbital variable names to return
                                     only those specific variables. Available
                                     orbital variables names are:
                                    ['e', 'l', 'x'].
                                     Note that these are available only for the
                                     inspiral portion of the waveform!
        failsafe                  -- If True, we make reasonable choices for the
                                     user, if the inputs to this method lead
                                     into exceptions.
        verbose                   -- Verbosity level. Available values are: 0, 1, 2

    Returns:
    --------
        modes_imr         -- Dictionary of IMR GW modes (PyCBC TimeSeries)
        orbital_var_dict  -- Dictionary of evolution of orbital elements.
                             Returned only if the flag `return_orbital_params`
                             is set
        retval            -- Hybridization related data. Returned only if the
                             flag `return_hybridization_info` is set
    """

    spin1z = 0.0
    spin2z = 0.0
    blend_using_avg_orbital_frequency = True
    modes_to_use = [(2, 2)]
    mode_to_align_by = (2, 2)

    if not hasattr(ls, merger_ringdown_approximant):
        raise IOError(
            """We cannot generate individual modes for {merger_ringdown_approximant}.
Try one of: [NRSur7dq4, SEOBNRv4PHM]"""
        )
    if (reference_mean_anomaly is None) and (coa_phase is None):
        raise IOError(
            f"""Please specify one of the phase angles, either of `reference_mean_anomaly` or `coa_phase`."""
        )
    if blend_aligning_merger_to_inspiral and (reference_mean_anomaly is None):
        raise IOError(
            f"""If you want to attach ESIGMASur inspiral to merger, by phase shifting merger to inspiral, please specify the phase angle `reference_mean_anomaly`"""
        )
    if (not blend_aligning_merger_to_inspiral) and (coa_phase is None):
        raise IOError(
            f"""If you want to attach ESIGMASur inspiral to merger, by phase shifting inspiral to merger, please specify the phase angle `coa_phase`"""
        )
    if reference_mean_anomaly is None:
        reference_mean_anomaly = 0
    if coa_phase is None:
        coa_phase = 0

    available_inspiral_orbital_params = ["e", "l", "x"]
    if return_orbital_params == True:
        return_orbital_params = available_inspiral_orbital_params
        return_orbital_params_user = set(return_orbital_params)
    elif (
        isinstance(return_orbital_params, list)
        or isinstance(return_orbital_params, set)
        or isinstance(return_orbital_params, tuple)
    ):
        return_orbital_params_user = return_orbital_params.copy()
        return_orbital_params_user = set(return_orbital_params_user).intersection(
            set(available_inspiral_orbital_params)
        )
        if return_orbital_params_user != set(return_orbital_params):
            print(
                f"""Warning: You requested the following list of orbital
parameters to be returned: {return_orbital_params}, but we reduce it to
{return_orbital_params_user} as we only have the evolution of the following 
parameters available with us: {available_inspiral_orbital_params}.
                  """
            )
    elif not return_orbital_params:
        return_orbital_params = []
        return_orbital_params_user = False

    return_orbital_params = set(return_orbital_params)
    return_orbital_params = return_orbital_params.union(
        set(["e"])
    )  # "e" needed necessarily for hybridization error printing

    # If the user does not provide the width of hybridization window (in terms
    # of orbital frequency) over which the inspiral should transition to
    # merger-ringdown, we switch schemes and blend over `num_hyb_orbits`
    # orbits instead.
    if f_window_mr_transition is None:
        return_orbital_params = return_orbital_params.union(set(["x"]))

    _, mode_to_align_by_em = mode_to_align_by

    # If the user does not provide the orbital frequency at which the inspiral
    # should transition to merger-ringdown, we use sensible defaults here.
    if f_mr_transition is None:
        # Kerr ISCO frequency
        f_Kerr = f_ISCO_spin(mass1, mass2, spin1z, spin2z)
        # Schwarzschild ISCO frequency
        f_Schwarz = 6.0**-1.5 / (mass1 + mass2) / lal.MTSUN_SI / lal.PI
        f_mr_transition = min(f_Kerr, f_Schwarz) * (mode_to_align_by_em / 2)

    retval = get_inspiral_esigmasur_modes(
        mass1=mass1,
        mass2=mass2,
        reference_eccentricity=reference_eccentricity,
        reference_mean_anomaly=reference_mean_anomaly,
        t_start=t_start,
        t_end=None,  # To generate the full inspiral for attachment to merger
        delta_t=delta_t,
        distance=distance,
        include_conjugate_modes=include_conjugate_modes,
        return_orbital_params=list(return_orbital_params),
        return_pycbc_timeseries=False,
        verbose=verbose,
    )

    # Retrieve modes, orbital phase and frequency from the returned list
    modes_inspiral_numpy = retval[-1]

    if mode_to_align_by not in modes_inspiral_numpy:
        raise RuntimeError(
            f"""The inspiral modes do not contain the primary 
desired {mode_to_align_by} multipole. It currently holds only the following:
{modes_inspiral_numpy.keys()}"""
        )

    orbital_eccentricity = retval[-2]["e"]
    # Throw error if eccentricity at the end of inspiral is definitely unsafe
    if orbital_eccentricity[-1] > ECCENTRICITY_LEVEL_ISCO_ERROR:
        raise IOError(
            f"""ERROR: You entered a very large reference eccentricity
{reference_eccentricity}. The orbital eccentricity at the end of inspiral was
{orbital_eccentricity[-1]}. The merger-ringdown attachment with a
quasicircular will be questionable."""
        )
    # Warn user if eccentricity at the end of inspiral is potentially unsafe
    if orbital_eccentricity[-1] > ECCENTRICITY_LEVEL_ISCO_WARNING and verbose:
        print(
            f"""WARNING: You entered a very large reference eccentricity
{reference_eccentricity}. The orbital eccentricity at the end of inspiral was
{orbital_eccentricity[-1]}. The merger-ringdown attachment with a quasicircular
model might be affected."""
        )

    if (f_window_mr_transition is None) or failsafe or (verbose > 1):
        if blend_using_avg_orbital_frequency:
            orbital_frequency = (
                retval[-2]["x"] ** 1.5 / ((mass1 + mass2) * lal.MTSUN_SI) / (2 * np.pi)
            )
        else:
            NotImplementedError(
                "Can't use any prescription other than the orbit averaged frequency one."
            )

    if return_orbital_params_user:
        orbital_vars_dict = {
            key: pt.TimeSeries(retval[-2][key], delta_t=delta_t, epoch=retval[0][0])
            for key in return_orbital_params_user
        }

    # DEBUG
    if verbose > 5:
        mode_phase = esigmapy.blend.compute_phase(
            modes_inspiral_numpy[mode_to_align_by]
        )
        mode_frequency = esigmapy.blend.compute_frequency(mode_phase, delta_t)
        print(
            f"""DEBUG: Orbital freq at end of inspiral is {orbital_frequency[-1]}Hz,
mode-{mode_to_align_by} freq at the end of inspiral is {mode_frequency[-1]}Hz, max and min
mode-{mode_to_align_by} frequencies are {np.max(mode_frequency)}Hz and
{np.min(mode_frequency)}Hz, and the transition frequency (of {mode_to_align_by}-mode)
requested is {f_mr_transition}Hz, which should be less than the maximum freq of
{mode_to_align_by}-mode: {mode_frequency.max()}Hz."""
        )
        return (
            modes_inspiral_numpy,
            mode_phase,
            mode_frequency,
            orbital_frequency,
            orbital_eccentricity,
            orbital_vars_dict,
        )

    # In case the user-specified transition frequency is too high, and they
    # requested failsafe mode, we reset it to a reasonable value.
    if failsafe:
        mode_phase = esigmapy.blend.compute_phase(
            modes_inspiral_numpy[mode_to_align_by]
        )
        mode_frequency = esigmapy.blend.compute_frequency(mode_phase, delta_t)
        if mode_frequency.max() < f_mr_transition:
            if verbose:
                print(
                    f"""FAILSAFE: Maximum orbital freq during inspiral is
{orbital_frequency.max()}Hz, and max frequency of {mode_to_align_by}-mode is
{mode_frequency.max()}Hz, so we are resetting transition frequency from
{f_mr_transition}Hz to {mode_frequency.max()}Hz."""
                )
            f_mr_transition = mode_frequency.max()

    # If the user does not provide the width of hybridization window (
    # `f_window_mr_transition`) over which the inspiral should transition to
    # merger-ringdown, we switch schemes and blend over `num_hyb_orbits`
    # orbits instead.
    if f_window_mr_transition is None:
        f_window_mr_transition = (
            _get_transition_frequency_window(
                None,
                orbital_frequency,
                delta_t,
                f_mr_transition / mode_to_align_by_em,
                num_hyb_orbits,
                keep_f_mr_transition_at_center,
                blend_using_avg_orbital_frequency,
                failsafe=failsafe,
                verbose=verbose,
            )
            * mode_to_align_by_em
        )

    # This is done to make use of the same hybridization code, that actually
    # assumes f_mr_transition to be at window's midpoint, to keep the
    # hybridization window's end at f_mr_transition
    if not keep_f_mr_transition_at_center:
        f_mr_transition -= f_window_mr_transition / 2.0

    # Generate NR surrogate waveform that will be our merger-ringdown, starting
    # from a frequency = 90% of
    max_retries = 20
    f_lower_mr = (f_mr_transition - f_window_mr_transition / 2) * (
        1.8 / mode_to_align_by_em
    )
    for _ in range(max_retries):
        try:
            if verbose:
                print(f"Generating MR waveform from {f_lower_mr}Hz...")
            hlm_mr = ls.SimInspiralChooseTDModes(
                coa_phase,  # phiRef
                delta_t,  # deltaT
                mass1 * lal.MSUN_SI,
                mass2 * lal.MSUN_SI,
                0,  # spin1x
                0,  # spin1y
                spin1z,
                0,  # spin2x
                0,  # spin2y
                spin2z,
                f_lower_mr,  # f_min
                f_lower_mr,  # f_ref
                distance * lal.PC_SI * 1.0e6,
                None,  # LALpars
                4,  # lmax
                getattr(ls, merger_ringdown_approximant),
            )
            break
        except:
            f_lower_mr *= 0.8
            continue

    # Extracting only the modes we need
    modes_to_use = list(modes_inspiral_numpy.keys())
    modes_mr_numpy = {}
    while hlm_mr is not None:
        key = (hlm_mr.l, hlm_mr.m)
        if key in modes_to_use:
            modes_mr_numpy[key] = hlm_mr.mode.data.data
        hlm_mr = hlm_mr.next

    try:
        retval = esigmapy.blend.blend_modes(
            modes_inspiral_numpy,
            modes_mr_numpy,
            orbital_frequency,
            f_mr_transition,
            frq_width=f_window_mr_transition,
            delta_t=delta_t,
            modes_to_blend=modes_to_use,
            mode_to_align_by=mode_to_align_by,
            blend_using_avg_orbital_frequency=blend_using_avg_orbital_frequency,
            blend_aligning_merger_to_inspiral=blend_aligning_merger_to_inspiral,
            include_conjugate_modes=include_conjugate_modes,
            verbose=verbose,
        )
    except Exception as exc:
        print(
            f"""Inspiral + MergerRingdown attachment failed. Its very likely
that you entered a very large reference eccentricity {reference_eccentricity}. The orbital
eccentricity at the end of inspiral was {orbital_eccentricity[-1]}
              """
        )
        raise exc
    modes_imr_numpy = retval[0]

    # Set t=0 at the end of inspiral surrogate
    if mode_to_align_by not in modes_imr_numpy:
        mode_to_align_by = list(modes_imr_numpy.keys())[0]
    idx_peak = len(modes_inspiral_numpy[mode_to_align_by]) - 1
    t_peak = idx_peak * delta_t

    itime = time.perf_counter()
    modes_imr = {}
    for el, em in modes_imr_numpy:
        modes_imr[(el, em)] = pt.TimeSeries(
            modes_imr_numpy[(el, em)], delta_t=delta_t, epoch=-1 * t_peak
        )
    if verbose > 4:
        print(
            "Time taken to store in pycbc.TimeSeries is {} secs".format(
                time.perf_counter() - itime
            )
        )

    if verbose:
        print("blended.")

    if return_hybridization_info and return_orbital_params_user:
        return retval, orbital_vars_dict, modes_imr
    elif return_orbital_params_user:
        return orbital_vars_dict, modes_imr
    elif return_hybridization_info:
        return retval, modes_imr
    return modes_imr


def get_imr_esigmasur_waveform(
    mass1,
    mass2,
    delta_t,
    reference_eccentricity=0.0,
    reference_mean_anomaly=0.0,
    t_start=None,
    distance=1.0,
    coa_phase=0.0,
    inclination=0.0,
    f_mr_transition=None,
    f_window_mr_transition=None,
    num_hyb_orbits=0.25,
    blend_aligning_merger_to_inspiral=True,
    keep_f_mr_transition_at_center=False,
    merger_ringdown_approximant="NRSur7dq4",
    return_hybridization_info=False,
    return_orbital_params=False,
    failsafe=True,
    verbose=False,
):
    """
    Returns IMR GW polarizations constructed using hybridized IMR ESIGMASur modes

    Parameters:
    -----------
        mass1, mass2              -- Binary's component masses (in solar masses)
        delta_t                   -- Waveform's time grid-spacing (in s)
        reference_eccentricity    -- Eccentricity at reference time of surrogate
        reference_mean_anomaly    -- Mean anomaly at reference time of surrogate (in rad)
        t_start                   -- Start time of the waveform to be generated (in seconds).
                                     Note that the surrogate defines t=0 at the end of
                                     inspiral, so t_start should be negative
                                     and t_start < 0.
                                     Defaults to the full duration of the surrogate.
        distance                  -- Luminosity distance to the binary (in Mpc)
        coa_phase                 -- Coalesence phase of the binary (in rad)
        inclination               -- Inclination (in rad), defined as the angle
                                     between the orbital angular momentum L and
                                     the line-of-sight
        f_mr_transition           -- Inspiral to merger transition GW frequency
                                     (Hz). Should correspond to the mode
                                     using which alignment is done, i.e. (2,2) mode here.
                                     Defaults to the minimum of the Kerr and
                                     Schwarzschild ISCO frequency equivalent
                                     for the (2,2)-mode.
        f_window_mr_transition    -- Hybridization frequency window (in Hz).
                                     Should correspond to the mode
                                     using which alignment is done, i.e. (2,2) mode here.
                                     Disabled by the default value (None). In
                                     such a case, the hybridization proceeds
                                     over a window of `num_hyb_orbits` orbital
                                     cycles (1 orbital cycle ~ 2 GW cycles)
                                     that ends at the frequency value given by
                                     `f_mr_transition`.
                                     Also see `keep_f_mr_transition_at_center`
                                     to choose the position of `f_mr_transition`
                                     within this window.
        num_hyb_orbits            -- number of orbital cycles to blend over.
                                     Only used if f_window_mr_transition is not
                                     specified.
        blend_aligning_merger_to_inspiral -- (default: False) If True, the
                                     merger-ringdown mode would be phase aligned
                                     to the inspiral
                                     If False, the inspiral is phase aligned
                                     Note: specify the desired
        keep_f_mr_transition_at_center -- If True, `f_mr_transition` is kept at
                                     the center of the hybridization window.
                                     Otherwise, it's kept at the end of the
                                     window (default).
        merger_ringdown_approximant    -- Choose merger-ringdown model. Tested
                                     choices: [NRSur7dq4, SEOBNRv4PHM]
        return_hybridization_info -- If True, returns hybridization related data
        return_orbital_params     -- If True, returns the orbital evolution of
                                     all the orbital elements (in
                                     geometrized units). Can also be a list of
                                     orbital variable names to return
                                     only those specific variables. Available
                                     orbital variables names are:
                                    ['e', 'l', 'x'].
                                     Note that these are available only for the
                                     inspiral portion of the waveform!
        failsafe                  -- If True, we make reasonable choices for the
                                     user, if the inputs to this method lead
                                     into exceptions.
        verbose                   -- Verbosity level. Available values are: 0, 1, 2

    Returns:
    --------
        hp, hc       -- Plus and cross IMR GW polarizations PyCBC TimeSeries
        orbital_vars_dict -- Dictionary of evolution of orbital elements.
                        Returned only if return_orbital_params is specified
        retval       -- Hybridization related data.
                        Returned only if return_hybridization_info is True
    """

    retval = get_imr_esigmasur_mode(
        mass1=mass1,
        mass2=mass2,
        reference_eccentricity=reference_eccentricity,
        reference_mean_anomaly=reference_mean_anomaly,
        delta_t=delta_t,
        t_start=t_start,
        distance=distance,
        coa_phase=coa_phase,
        include_conjugate_modes=True,  # Always include conjugate modes while generating polarizations
        f_mr_transition=f_mr_transition,
        f_window_mr_transition=f_window_mr_transition,
        num_hyb_orbits=num_hyb_orbits,
        blend_aligning_merger_to_inspiral=blend_aligning_merger_to_inspiral,
        keep_f_mr_transition_at_center=keep_f_mr_transition_at_center,
        merger_ringdown_approximant=merger_ringdown_approximant,
        return_hybridization_info=return_hybridization_info,
        return_orbital_params=return_orbital_params,
        failsafe=failsafe,
        verbose=verbose,
    )
    if return_hybridization_info and return_orbital_params:
        modes_imr, orbital_vars_dict, retval = retval
    elif return_hybridization_info:
        modes_imr, retval = retval
    elif return_orbital_params:
        modes_imr, orbital_vars_dict = retval
    else:
        modes_imr = retval

    hp, hc = esigmapy.utils.get_polarizations_from_multipoles(
        modes_imr,
        inclination=inclination,
        coa_phase=np.pi / 2 - coa_phase,
        verbose=verbose,
    )

    if return_hybridization_info and return_orbital_params:
        return hp, hc, orbital_vars_dict, retval
    elif return_hybridization_info:
        return hp, hc, retval
    elif return_orbital_params:
        return hp, hc, orbital_vars_dict
    return hp, hc
