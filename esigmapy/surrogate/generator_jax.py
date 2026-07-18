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

from .. import blend_jax

_TWO_PI = 2.0 * np.pi


# --- transition frequency window ---------------------------------------------


def _window_accumulate_backward_jax(contrib, delta_phi):
    """Backward accumulate-until-threshold of ``_get_window_start_numba``.

    ``contrib[k] = 0.5*(f[k]+f[k+1])*delta_t`` for the series being
    integrated; visits ``k = n-2, n-3, ..., 1`` (never 0, as in the numba
    loop) and stops at the first ``k`` where the accumulated ``|sum|``
    reaches ``|delta_phi|``. Returns ``(idx, found)``.
    """
    vals = jnp.cumsum(contrib[1:][::-1])
    hit = jnp.abs(vals) >= jnp.abs(delta_phi)
    n = contrib.shape[0] + 1
    idx = n - 2 - jnp.argmax(hit)
    return idx, jnp.any(hit)


def _window_accumulate_forward_jax(contrib, delta_phi):
    """Forward variant: visits ``idx = 1..n-1`` accumulating ``contrib[idx-1]``
    and stops when ``|sum|`` reaches ``|delta_phi|``. Returns ``(idx, found)``."""
    vals = jnp.cumsum(contrib)
    hit = jnp.abs(vals) >= jnp.abs(delta_phi)
    return jnp.argmax(hit) + 1, jnp.any(hit)


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
