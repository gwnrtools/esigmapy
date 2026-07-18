# Copyright (C) 2026 Akash Maurya

"""JAX implementation of the ESIGMASur inspiral-merger-ringdown (IMR) pipeline.

This module composes the JAX inspiral surrogate (``surrogate_jax``), the JAX
port of the mode hybridization (``esigmapy.blend_jax``), a JAX port of the
transition-frequency-window computation (``esigmapy.generator``), and the JAX
NRSur7dq4 merger-ringdown surrogate from ``gwsurrogate.jax`` into the IMR
waveform construction of ``get_imr_esigmasur_mode`` /
``get_imr_esigmasur_waveform`` (``esigmapy/surrogate/generator.py``).

Design notes
------------
* ``surrogate_jax`` must be imported first in the process (its ``TPI_jax``
  import enables JAX 64-bit mode); this module imports it at the top for that
  reason. ``gwsurrogate.jax`` sets the same flag, so import order between the
  two is safe.
* All data-dependent constructs (transition index, window start/end, blend
  window length, hybrid length) are traced values inside fixed-shape buffers;
  see ``esigmapy.blend_jax`` for the blending kernels. Failsafe clamping
  replaces the NumPy error paths inside traces; host wrappers re-validate on
  concrete values and raise the NumPy backend's errors.
* The merger-ringdown piece is evaluated over the full NRSur7dq4 span (no
  ``f_lower`` truncation or retry loop): the blend discards samples before its
  window and re-aligns time and phase, so the lalsimulation ``f_lower``/
  ``f_ref`` conventions of the NumPy backend need not be reproduced.
"""

import os
import sys

import numpy as np

from . import surrogate_jax  # noqa: F401  (enables jax x64 via TPI_jax)
import jax
import jax.numpy as jnp

import lal

from .. import blend_jax

_TWO_PI = 2.0 * np.pi


def _import_gwsurrogate_jax():
    """Import ``gwsurrogate.jax.surrogate`` (the JAX NRSur7dq4 port).

    The JAX port lives only in the gwsurrogate source tree; if the package is
    not importable (or importable without the ``jax`` subpackage), locate the
    source tree via ``$GWSURROGATE_JAX_PATH`` or the ``gwsurrogate`` directory
    sitting beside the esigmapy repository."""
    try:
        from gwsurrogate.jax import surrogate as gws_jax_surrogate

        return gws_jax_surrogate
    except ModuleNotFoundError:
        pass
    candidates = []
    env = os.environ.get("GWSURROGATE_JAX_PATH")
    if env:
        candidates.append(env)
    repo_root = os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    )
    candidates.append(os.path.join(os.path.dirname(repo_root), "gwsurrogate"))
    for cand in candidates:
        if os.path.isdir(os.path.join(cand, "gwsurrogate", "jax")):
            sys.path.insert(0, cand)
            break
    from gwsurrogate.jax import surrogate as gws_jax_surrogate

    return gws_jax_surrogate


# --- merger-ringdown: NRSur7dq4 ----------------------------------------------

# NRSur7dq4 mode-row ordering for ell_max=2: (2,-2), (2,-1), (2,0), (2,1), (2,2)
_MODE_2M2_ROW = 0
_MODE_22_ROW = 4


class NRSurMergerRingdownJAX:
    """Non-spinning NRSur7dq4 (2, +-2) modes on a uniform physical time grid.

    The merger-ringdown piece of the JAX IMR pipeline. Modes are evaluated
    over the full NRSur7dq4 span (~[-4300, +100] M around the peak) -- no
    ``f_lower`` truncation; the blend discards samples ahead of its window.
    """

    def __init__(self, h5_path=None):
        self._gws = _import_gwsurrogate_jax()
        self._sur = self._gws.NRSur7dq4JAX(h5_path)
        self.data = self._sur.data
        self.t_coorb = np.asarray(self._sur.t_coorb, dtype=np.float64)

    def num_samples(self, M, delta_t):
        """Number of ``delta_t`` samples spanning the full surrogate (host)."""
        span = (self.t_coorb[-1] - self.t_coorb[0]) * M * lal.MTSUN_SI
        return int(np.floor(span / delta_t)) + 1

    def modes_2pm2(self, q, M, delta_t, distance, n_mr):
        """(2,2) and (2,-2) modes; traceable in ``q``, ``M``, ``distance``.

        ``delta_t`` (s) and the static sample count ``n_mr`` (which must keep
        the grid inside the surrogate span, cf. ``num_samples``) define the
        uniform output grid starting at the surrogate's first sample. Returns
        a dict {(2, 2): ..., (2, -2): ...} of complex (n_mr,) arrays scaled to
        ``M`` solar masses at ``distance`` Mpc (lalsimulation convention).
        """
        zeros3 = jnp.zeros(3, dtype=jnp.float64)
        init_quat = jnp.asarray(self._gws._IDENTITY_QUATERNION, dtype=jnp.float64)
        h, _, _ = self._gws._evaluate_dimensionless_modes(
            self.data, q, zeros3, zeros3, init_quat, 0.0, 2
        )
        rows = h[jnp.array([_MODE_22_ROW, _MODE_2M2_ROW])]
        seconds_per_M = M * lal.MTSUN_SI
        times_M = self.t_coorb[0] + jnp.arange(n_mr) * (delta_t / seconds_per_M)
        res = self._gws._resample_modes(self.data, rows, times_M)
        amp_scale = M * lal.MRSUN_SI / (distance * 1.0e6 * lal.PC_SI)
        res = res * amp_scale
        return {(2, 2): res[0], (2, -2): res[1]}


# --- transition frequency window ---------------------------------------------


def transition_frequency_window_jax(
    orbital_freq,
    delta_t,
    f_mr_transition,
    num_hyb_orbits,
    keep_f_mr_transition_at_center=False,
):
    """JAX port of ``esigmapy.generator._get_transition_frequency_window`` for
    the ``blend_using_avg_orbital_frequency=True`` scheme with failsafe
    clamping (the traced pipeline cannot raise; callers with concrete values
    should inspect the returned ``found`` flags to reproduce the non-failsafe
    errors).

    ``keep_f_mr_transition_at_center`` is static; ``orbital_freq`` (Hz,
    orbital), ``f_mr_transition`` (Hz, orbital) and the rest may be traced.

    Returns ``(f_window, window_start_found, window_end_found)``;
    ``window_end_found`` is always True for the end-anchored scheme.
    """
    orbital_freq = jnp.asarray(orbital_freq, dtype=jnp.float64)
    n = orbital_freq.shape[0]
    # last index where orbital_freq drops below the transition frequency
    transition_idx = n - 1 - jnp.argmax(orbital_freq[::-1] < f_mr_transition)

    contrib = 0.5 * (orbital_freq[:-1] + orbital_freq[1:]) * delta_t
    k = jnp.arange(n - 1)

    if keep_f_mr_transition_at_center:
        # backward from transition_idx over orbital_freq[:transition_idx+1]:
        # numba visits k = transition_idx-1 .. 1; masked contributions keep
        # the same accumulation order with zeros outside the slice
        contrib_back = jnp.where(k <= transition_idx - 1, contrib, 0.0)
        vals_back = jnp.cumsum(contrib_back[1:][::-1])
        hit_back = jnp.abs(vals_back) >= jnp.abs(num_hyb_orbits * jnp.pi)
        window_start_idx = n - 2 - jnp.argmax(hit_back)
        start_found = jnp.any(hit_back)
        window_start_idx = jnp.where(start_found, window_start_idx, 0)

        # forward from transition_idx over orbital_freq[transition_idx:]:
        # numba visits slice indices 1.. accumulating contrib[transition_idx + j]
        contrib_fwd = jnp.where(k >= transition_idx, contrib, 0.0)
        vals_fwd = jnp.cumsum(jnp.roll(contrib_fwd, -transition_idx))
        hit_fwd = jnp.abs(vals_fwd) >= jnp.abs(num_hyb_orbits * jnp.pi)
        window_end_idx = transition_idx + jnp.argmax(hit_fwd) + 1
        end_found = jnp.any(hit_fwd)
        window_end_idx = jnp.where(end_found, window_end_idx, n - 1)

        f_transition_at = f_mr_transition
        f_window = 2.0 * jnp.minimum(
            jnp.take(orbital_freq, window_end_idx) - f_transition_at,
            f_transition_at - jnp.take(orbital_freq, window_start_idx),
        )
        return f_window, start_found, end_found

    # end-anchored (default): backward from transition_idx with twice the
    # phase threshold
    contrib_back = jnp.where(k <= transition_idx - 1, contrib, 0.0)
    vals_back = jnp.cumsum(contrib_back[1:][::-1])
    hit_back = jnp.abs(vals_back) >= jnp.abs(2.0 * num_hyb_orbits * jnp.pi)
    window_start_idx = n - 2 - jnp.argmax(hit_back)
    start_found = jnp.any(hit_back)
    window_start_idx = jnp.where(start_found, window_start_idx, 0)

    f_window = jnp.take(orbital_freq, transition_idx) - jnp.take(
        orbital_freq, window_start_idx
    )
    return f_window, start_found, jnp.asarray(True)


# --- IMR pipeline -------------------------------------------------------------


class IMRESIGMASurJAX:
    """JAX inspiral-merger-ringdown ESIGMASur: eccentric inspiral surrogate
    blended with the NRSur7dq4 merger-ringdown.

    The JAX counterpart of ``get_imr_esigmasur_mode``
    (``esigmapy/surrogate/generator.py``): (2,2)-mode, non-spinning, blending
    with the orbit-averaged orbital frequency derived from the surrogate's PN
    parameter ``x``. Differences from the NumPy backend:

    * the merger-ringdown is the (JAX) gwsurrogate NRSur7dq4 evaluated over
      its full span, not lalsimulation's NRSur7dq4 from ``f_lower`` -- the
      blend re-aligns time and phase, and the two NRSur7dq4 implementations
      agree to ~2e-4 in amplitude / ~1.5e-3 rad in phase (see
      tests/test_imr_jax.py Tier C), so blended waveforms agree at that level
      plus discrete window-snapping effects;
    * outputs are plain arrays ``(times, modes_dict)`` (the NumPy backend
      returns pycbc TimeSeries), with the same epoch convention (t=0 at the
      end of inspiral);
    * the circular fallback at ``e_ref <= e_ref_min`` is not supported.

    ``__call__`` is the host-orchestrated path (eager JAX kernels, concrete
    shapes); use ``imr_parameter_space_evaluator`` for a pure traced function
    composable with jit/vmap/grad.
    """

    def __init__(self, ecc_data_dir, circ_data_dir, nrsur_h5_path=None):
        self.ecc_sur = surrogate_jax.EccentricSurrogateJAX(
            ecc_data_dir=ecc_data_dir, circ_data_dir=circ_data_dir
        )
        self.mr_sur = NRSurMergerRingdownJAX(nrsur_h5_path)

    @staticmethod
    def _f_mr_transition_default(M, q):
        """min(Kerr, Schwarzschild) ISCO (2,2)-mode frequency (Hz), as in
        ``get_imr_esigmasur_mode`` with zero spins."""
        from ..utils import f_ISCO_spin

        mass1 = M * q / (1.0 + q)
        mass2 = M / (1.0 + q)
        f_Kerr = f_ISCO_spin(mass1, mass2, 0.0, 0.0)
        f_Schwarz = 6.0**-1.5 / M / lal.MTSUN_SI / lal.PI
        return min(f_Kerr, f_Schwarz)  # * (em / 2) with em = 2

    def __call__(
        self,
        M,
        params,
        delta_t,
        t_start=None,
        distance=1.0,
        coa_phase=None,
        include_conjugate_modes=False,
        f_mr_transition=None,
        f_window_mr_transition=None,
        num_hyb_orbits=0.25,
        blend_aligning_merger_to_inspiral=True,
        keep_f_mr_transition_at_center=False,
        return_hybridization_info=False,
        return_orbital_params=False,
        failsafe=True,
        merger_ringdown_modes=None,
        verbose=False,
    ):
        """IMR (2,2) [and (2,-2)] modes for ``M`` (solar masses) and
        ``params = (q, e_ref, l_ref)`` on a uniform ``delta_t`` grid.

        Mirrors ``get_imr_esigmasur_mode`` semantics (see its docstring for
        the parameters). ``merger_ringdown_modes`` optionally injects
        precomputed MR mode arrays (dict keyed by (l, m), same ``delta_t``) in
        place of the NRSur7dq4-JAX evaluation -- used by parity tests and
        advanced callers.

        Returns ``(times, modes)`` with ``times`` anchored so t=0 is the end
        of inspiral, prepended to the NumPy backend's optional returns:
        ``(times, hyb_info, orbital_var_dict, modes)`` with both flags,
        ``(times, hyb_info, modes)`` / ``(times, orbital_var_dict, modes)``
        with one.
        """
        from ..generator import (
            ECCENTRICITY_LEVEL_ISCO_ERROR,
            ECCENTRICITY_LEVEL_ISCO_WARNING,
        )

        q, e_ref, l_ref = params
        em = 2  # mode_to_align_by = (2, 2)

        if (not blend_aligning_merger_to_inspiral) and (coa_phase is None):
            raise IOError(
                """If you want to attach ESIGMASur inspiral to merger, by phase shifting inspiral to merger, please specify the phase angle `coa_phase`"""
            )
        if coa_phase is None:
            coa_phase = 0.0
        if e_ref <= self.ecc_sur.e_ref_min:
            raise NotImplementedError(
                "The JAX IMR pipeline requires e_ref > e_ref_min (the "
                "circular fallback of the NumPy backend is not supported)."
            )

        available = ["e", "l", "x"]
        if return_orbital_params is True:
            return_orbital_params_user = set(available)
        elif isinstance(return_orbital_params, (list, set, tuple)):
            return_orbital_params_user = set(return_orbital_params) & set(available)
        else:
            return_orbital_params_user = False

        t_insp, orb_vars, h22 = self.ecc_sur(
            M,
            [q, e_ref, l_ref],
            delta_t=delta_t,
            t_start=t_start,
            t_end=None,
            remove_start_phase=True,
            return_orbital_variables=True,
        )
        h22 = np.asarray(h22) / distance
        modes_insp = {(2, 2): h22}
        if include_conjugate_modes:
            modes_insp[(2, -2)] = np.conj(h22)  # (-1)**l conj for l=2

        ecc_end = float(orb_vars["e"][-1])
        if ecc_end > ECCENTRICITY_LEVEL_ISCO_ERROR:
            raise IOError(f"""ERROR: You entered a very large reference eccentricity
{e_ref}. The orbital eccentricity at the end of inspiral was
{ecc_end}. The merger-ringdown attachment with a
quasicircular will be questionable.""")
        if ecc_end > ECCENTRICITY_LEVEL_ISCO_WARNING and verbose:
            print(f"""WARNING: You entered a very large reference eccentricity
{e_ref}. The orbital eccentricity at the end of inspiral was
{ecc_end}. The merger-ringdown attachment with a quasicircular
model might be affected.""")

        orbital_frequency = (
            np.asarray(orb_vars["x"]) ** 1.5 / (M * lal.MTSUN_SI) / (2 * np.pi)
        )

        if f_mr_transition is None:
            f_mr_transition = self._f_mr_transition_default(M, q)

        if failsafe:
            mode_frequency = blend_jax.compute_frequency_jax(
                blend_jax.compute_phase_jax(jnp.asarray(h22)), delta_t
            )
            max_mode_frequency = float(jnp.max(mode_frequency))
            if max_mode_frequency < f_mr_transition:
                if verbose:
                    print(f"""FAILSAFE: resetting transition frequency from
{f_mr_transition}Hz to {max_mode_frequency}Hz.""")
                f_mr_transition = max_mode_frequency

        if f_window_mr_transition is None:
            f_win, start_found, end_found = transition_frequency_window_jax(
                jnp.asarray(orbital_frequency),
                delta_t,
                f_mr_transition / em,
                num_hyb_orbits,
                keep_f_mr_transition_at_center,
            )
            if not failsafe and not bool(start_found):
                raise RuntimeError(
                    f"""Requested number of orbits to blend over not available
in the waveform before the transition frequency. Either decrease the number of
orbits to blend over (currently {num_hyb_orbits}) or increase the
inspiral-to-merger transition frequency. `window_start_idx` is None."""
                )
            if not failsafe and not bool(end_found):
                raise RuntimeError(
                    f"""Requested number of orbits to blend over not available
in the waveform after the transition frequency. Either decrease the number of
orbits to blend over (currently {num_hyb_orbits}) or decrease the
inspiral-to-merger transition frequency. `window_end_idx` is None."""
                )
            f_window_mr_transition = float(f_win) * em

        # Reuse the centered-window blend with the window's right end at
        # f_mr_transition (as in the NumPy driver).
        if not keep_f_mr_transition_at_center:
            f_mr_transition -= f_window_mr_transition / 2.0

        if merger_ringdown_modes is None:
            n_mr = self.mr_sur.num_samples(M, delta_t)
            modes_mr = self.mr_sur.modes_2pm2(q, M, delta_t, distance, n_mr)
            modes_mr = {
                lm: np.asarray(h) * np.exp(-1j * lm[1] * coa_phase)
                for lm, h in modes_mr.items()
                if lm in modes_insp
            }
        else:
            modes_mr = merger_ringdown_modes

        retval = blend_jax.blend_modes_jax(
            modes_insp,
            modes_mr,
            orbital_frequency,
            f_mr_transition,
            frq_width=f_window_mr_transition,
            delta_t=delta_t,
            modes_to_blend=list(modes_insp.keys()),
            mode_to_align_by=(2, 2),
            blend_using_avg_orbital_frequency=True,
            blend_aligning_merger_to_inspiral=blend_aligning_merger_to_inspiral,
            include_conjugate_modes=include_conjugate_modes,
            verbose=verbose,
        )
        modes_imr = retval[0]

        # t=0 at the end of inspiral (epoch convention of the NumPy backend)
        t_peak = (len(h22) - 1) * delta_t
        n_imr = len(modes_imr[(2, 2)])
        times = -t_peak + np.arange(n_imr) * delta_t

        if return_orbital_params_user:
            orbital_vars_dict = {
                key: np.asarray(orb_vars[key]) for key in return_orbital_params_user
            }
            if return_hybridization_info:
                return times, retval, orbital_vars_dict, modes_imr
            return times, orbital_vars_dict, modes_imr
        if return_hybridization_info:
            return times, retval, modes_imr
        return times, modes_imr
