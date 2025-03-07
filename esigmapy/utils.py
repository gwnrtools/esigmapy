# Copyright (C) 2025 Kaushik Paul
#
from __future__ import absolute_import, print_function

import lal
import numpy as np


def f_ISCO_spin(mass1, mass2, spin1z, spin2z):
    """
    Kerr ISCO frequency fitting formula for aligned-spins binary black holes

    Parameters:
    ----------
    mass1, mass2   -- Binary's component masses (in solar masses)
    spin1z, spin2z   -- z-components of component dimensionless spins (lies in [0,1))

    Returns: Kerr ISCO frequency (in Hz)
    -------
    """
    Msun = lal.MTSUN_SI
    k00 = -3.821158961
    k01 = -1.2019
    k02 = -1.20764
    k10 = 3.79245
    k11 = 1.18385
    k12 = 4.90494
    Zeta = 0.41616

    m1 = mass1 * Msun
    m2 = mass2 * Msun
    M = m1 + m2
    eta = (m1 * m2) / (M**2)
    z = 0

    atot = (spin1z + spin2z * (m2 / m1) ** 2) / ((1 + (m2 / m1)) ** 2)
    aeff = atot + Zeta * eta * (spin1z + spin2z)

    Z1 = 1 + ((1 - aeff**2) ** (1 / 3)) * (
        ((1 + aeff) ** (1 / 3)) + ((1 - aeff) ** (1 / 3))
    )
    Z2 = np.sqrt(3 * aeff**2 + Z1**2)
    # riscocap = 3 + Z2 - ((aeff) / abs(aeff)) * np.sqrt((3 - Z1) * (3 + Z1 + 2 * Z2))
    riscocap = 3 + Z2 - np.sign(aeff) * np.sqrt((3 - Z1) * (3 + Z1 + 2 * Z2))
    # Equivalent to the above commented expression, and works also for aeff=0.
    # For aeff=0, np.sign(aeff)=0, which is alright because its coefficient
    # np.sqrt((3 - Z1) * (3 + Z1 + 2 * Z2)) is anyway 0 when aeff=0 because Z1 = 3 then.
    Eiscocap = np.sqrt((1 - (2 / (3 * riscocap))))
    Liscocap = (2 / (3 * np.sqrt(3))) * (1 + 2 * np.sqrt(3 * riscocap - 2))

    chichif = (
        atot
        + eta * (Liscocap - 2 * atot * (Eiscocap - 1))
        + k00 * eta**2
        + k01 * eta**2 * aeff
        + k02 * eta**2 * aeff**2
        + k10 * eta**3
        + k11 * eta**3 * aeff
        + k12 * eta**3 * aeff**2
    )
    chichip = chichif

    Z1p = 1 + ((1 - chichip**2) ** (1 / 3)) * (
        ((1 + chichip) ** (1 / 3)) + ((1 - chichip) ** (1 / 3))
    )
    Z2p = np.sqrt(3 * chichip**2 + Z1p**2)
    riscocapp = (
        3 + Z2p - ((chichip) / abs(chichip)) * np.sqrt((3 - Z1p) * (3 + Z1p + 2 * Z2p))
    )

    omegacap = 1 / (riscocapp ** (3 / 2) + chichip)
    Scap = (1 / (1 - 2 * eta)) * (spin1z * m1**2 + spin2z * m2**2) / (M**2)
    SS = (1 + Scap * (-0.00303023 - 2.00661 * eta + 7.70506 * eta**2)) / (
        1 + Scap * (-0.67144 - 1.475698 * eta + 7.30468 * eta**2)
    )

    # Erad = (M* (0.0559745 * eta + 0.580951 * eta**2 - 0.960673 * eta**3 + 3.35241 * eta**4)* SS)
    Erad_by_M = (
        0.0559745 * eta + 0.580951 * eta**2 - 0.960673 * eta**3 + 3.35241 * eta**4
    ) * SS

    # Mfin = M * (1 - Erad / M)
    Mfin = M * (1 - Erad_by_M)
    fre = (1 / (1 + z)) * omegacap / (np.pi * Mfin)

    return 1 / 2 * fre


def get_peak_freqs(freq):
    """
    Inputs
    ------
    freq: Array of similar iterable of frequency values

    Outputs
    -------
    peak_times: Times at which local maxima of frequency are attained
    peak_freqs: Frequency at these times
    """
    peaks, peak_times = [], []
    fvals = freq.data

    for idx, finst in enumerate(fvals):
        if idx == 0 or idx == len(fvals) - 1:
            continue

        if ((fvals[idx - 1]) < (finst)) and ((fvals[idx + 1]) < (finst)):
            peaks.append(finst)
            peak_times.append(freq.sample_times[idx])

    return np.array(peak_times), np.array(peaks)


def get_polarizations_from_multipoles(
    waveform_multipoles, inclination, coa_phase, verbose=False
):
    """
    Returns GW polarizations from complex GW multipoles

    Parameters:
    -----------
        waveform_multipoles     -- Dictionary of complex Waveform Multipoles
                                   indexed by their `(el, em)` values
        inclination             -- Inclination (in rad), defined as the angle between the orbital
                                   angular momentum L and the line-of-sight
        coa_phase               -- Coalesence phase of the binary (in rad)
        delta_t                 -- Waveform's time grid-spacing (in s)
        verbose                 -- Verbosity level. Available values are: 0, 1, 2

    Returns:
    --------
        hp, hc            -- Plus and cross GW polarizations
    """
    available_modes = list(waveform_multipoles.keys())

    # Initialize with zeros, preserving data type and precision
    try:
        hp = waveform_multipoles[available_modes[0]].real * 0
        hc = waveform_multipoles[available_modes[0]].real * 0
    except:
        hp = waveform_multipoles[available_modes[0]].real() * 0
        hc = waveform_multipoles[available_modes[0]].real() * 0

    for el, em in waveform_multipoles:
        ylm = lal.SpinWeightedSphericalHarmonic(inclination, coa_phase, -2, el, em)
        glm = waveform_multipoles[(el, em)] * ylm
        if verbose > 4:
            print(f"Adding mode {el}, {em} with ylm = {ylm}", flush=True)
            print(
                f"... adding {waveform_multipoles[(el, em)]}, {glm}",
                flush=True,
            )
            print(f"... after adding, hp={hp}, hc={hc}", flush=True)
        try:
            hp = hp + glm.real
            hc = hc - glm.imag
        except:
            hp = hp + glm.real()
            hc = hc - glm.imag()

    return hp, hc
