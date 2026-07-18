# Copyright (C) 2026 Akash Maurya

"""JAX port of ``blend.py``: waveform-mode hybridization as pure, jit-safe
functions.

The port follows ``blend.py`` operation-for-operation so outputs agree at
float64 round-off level, but restructures the data-dependent constructs for
``jax.jit``:

* the bracket searches return traced indices (masked argmax) plus ``found``
  flags instead of raising;
* the blending window lives in fixed-length buffers of static size
  ``n_win_buf + 1`` (``n_win_buf`` >= the traced window length ``n_win``),
  with entries beyond ``n_win`` masked/ignored;
* the hybrid mode is assembled into a fixed-length buffer of static size
  ``len(h_insp) + len(h_mr)`` via clipped gathers, zero-padded beyond the
  data-dependent valid length ``t1_insp + len(h_mr) - t1_mr``.

``blend_modes_jax`` wraps the kernels behind the same dict-in/dict-out
signature (and input validation) as ``blend.blend_modes`` for parity testing;
the traced IMR pipeline calls the kernels directly.
"""

import numpy as np

import jax

jax.config.update("jax_enable_x64", True)  # before any jax array is created
import jax.numpy as jnp

_TWO_PI = 2.0 * np.pi


def compute_phase_jax(waveform):
    """Phase convention of ``blend.compute_phase``: ``-angle`` then unwrap."""
    angle = -jnp.angle(waveform)
    if angle.shape[0] <= 1:
        return angle.astype(jnp.float64)
    diff = jnp.diff(angle)
    diff = jnp.mod(diff + jnp.pi, _TWO_PI) - jnp.pi
    return jnp.concatenate([angle[:1], angle[0] + jnp.cumsum(diff)])


def compute_frequency_jax(phase, delta_t):
    """Central-difference frequency of ``blend.compute_frequency`` (one-sided
    differences at the ends), in cycles per unit time."""
    n = phase.shape[0]
    if n == 0:
        return phase
    if n == 1:
        return jnp.zeros(1, dtype=jnp.float64)
    factor = 1.0 / (_TWO_PI * delta_t)
    first = (phase[1] - phase[0]) * factor
    last = (phase[-1] - phase[-2]) * factor
    interior = (phase[2:] - phase[:-2]) * (0.5 * factor)
    return jnp.concatenate([first[None], interior, last[None]])


def find_first_index_jax(frq_timeseries, frq_desired):
    """First left-to-right bracket of ``frq_desired``, snapped to the nearer
    endpoint (``blend._find_first_value_location_in_series_numba``).

    Returns ``(idx, found)``; ``idx`` is meaningful only where ``found``.
    """
    f1 = frq_timeseries[:-1]
    f2 = frq_timeseries[1:]
    mask = (f1 <= frq_desired) & (f2 >= frq_desired)
    idx = jnp.argmax(mask)
    take_left = jnp.abs(frq_desired - f1[idx]) <= jnp.abs(frq_desired - f2[idx])
    return jnp.where(take_left, idx, idx + 1), jnp.any(mask)


def find_last_index_jax(frq_timeseries, frq_desired):
    """Last (rightmost) bracket of ``frq_desired``, snapped to the nearer
    endpoint (``blend._find_last_value_location_in_series_numba``).

    Returns ``(idx, found)``; ``idx`` is meaningful only where ``found``.
    """
    n = frq_timeseries.shape[0]
    # Bracket at i1 in [1, n-1]: f[i1] >= target and f[i1-1] <= target.
    mask = (frq_timeseries[1:] >= frq_desired) & (frq_timeseries[:-1] <= frq_desired)
    i1 = n - 1 - jnp.argmax(mask[::-1])
    take_right = jnp.abs(frq_desired - frq_timeseries[i1]) <= jnp.abs(
        frq_desired - frq_timeseries[i1 - 1]
    )
    return jnp.where(take_right, i1, i1 - 1), jnp.any(mask)


def locate_blend_indices_jax(
    frq_mr_align, inspiral_orbital_frequency, frq_attach, frq_width, em
):
    """Window indices of ``blend.blend_modes`` (avg-orbital-frequency scheme).

    ``frq_mr_align`` is the alignment mode's frequency over the full
    merger-ringdown series. Returns
    ``(t1_insp, t2_insp, t1_mr, t2_mr, found)`` with traced ints and the
    conjunction of the three bracket-search ``found`` flags.
    """
    t1_mr, ok1 = find_first_index_jax(frq_mr_align, frq_attach - frq_width / 2)
    t2_mr, ok2 = find_first_index_jax(frq_mr_align, frq_attach + frq_width / 2)
    t2_insp, ok3 = find_last_index_jax(
        inspiral_orbital_frequency, (frq_attach + frq_width / 2) / em
    )
    t1_insp = t2_insp - (t2_mr - t1_mr)
    return t1_insp, t2_insp, t1_mr, t2_mr, ok1 & ok2 & ok3


def windowed_mode_frequency_jax(h, delta_t, window_start_idx, n_win_buf):
    """Mode frequency over window positions ``0..n_win`` from a local phase
    reconstruction, replicating blend.py's slice ``[t1-1 : t2+2]`` followed by
    ``compute_frequency(...)[1:-1]``.

    Returns a ``(n_win_buf + 1,)`` buffer; entries beyond the traced window
    length are garbage (clipped gathers) and must be masked by the caller.
    """
    i = jnp.arange(n_win_buf + 3)
    idx = jnp.clip(window_start_idx - 1 + i, 0, h.shape[0] - 1)
    phase = compute_phase_jax(jnp.take(h, idx))
    return compute_frequency_jax(phase, delta_t)[1:-1]


def blend_single_mode_jax(
    h_insp,
    h_mr,
    delta_t,
    t1_insp,
    t1_mr,
    n_win,
    frq_insp_window,
    frq_mr_window,
    n_win_buf,
    align_merger_to_inspiral,
):
    """Blend one mode given window indices and windowed frequencies.

    Follows the per-mode loop of ``blend.blend_modes`` (amplitude and
    frequency sin^2-blended over the window, phase from trapezoidal
    integration of the blended frequency, anchored to the inspiral at the
    window start or to the merger-ringdown at the window end).

    ``n_win_buf`` and ``align_merger_to_inspiral`` are static; all other
    arguments may be traced. Returns ``(hyb, valid_len, amp_insp_w, amp_mr_w,
    amp_h, frq_h, phase_h)`` where ``hyb`` has static length
    ``len(h_insp) + len(h_mr)``, zeroed at and beyond ``valid_len``, and the
    window buffers have static length ``n_win_buf + 1`` (entries beyond
    ``n_win`` are garbage to be masked by the caller).
    """
    n_insp = h_insp.shape[0]
    n_mr = h_mr.shape[0]
    i = jnp.arange(n_win_buf + 1)
    idx_insp = jnp.clip(t1_insp + i, 0, n_insp - 1)
    idx_mr = jnp.clip(t1_mr + i, 0, n_mr - 1)

    amp_insp_w = jnp.abs(jnp.take(h_insp, idx_insp))
    amp_mr_w = jnp.abs(jnp.take(h_mr, idx_mr))

    tau = jnp.sin((jnp.pi / 2) * i / n_win) ** 2
    amp_h = (1.0 - tau) * amp_insp_w + tau * amp_mr_w
    frq_h = (1.0 - tau) * frq_insp_window + tau * frq_mr_window

    # Trapezoidal integration: phase = 2*pi * int f dt
    pair = 0.5 * (frq_h[:-1] + frq_h[1:]) * delta_t
    phase_h = jnp.concatenate([jnp.zeros(1), jnp.cumsum(pair)]) * _TWO_PI

    inspiral_angle_at_window_start = -jnp.angle(jnp.take(h_insp, t1_insp))
    t2_mr = t1_mr + n_win
    mr_angle_at_window_end = -jnp.angle(jnp.take(h_mr, t2_mr))
    phase_h_end = jnp.take(phase_h, n_win)

    if align_merger_to_inspiral:
        phase_h = phase_h + inspiral_angle_at_window_start
        mr_phasor_shift = jnp.exp(
            -1j
            * (phase_h_end + inspiral_angle_at_window_start - mr_angle_at_window_end)
        )
        head_phasor = jnp.ones((), dtype=h_insp.dtype)
        tail_phasor = mr_phasor_shift
    else:
        phase_h = phase_h + (mr_angle_at_window_end - phase_h_end)
        inspiral_phasor_shift = jnp.exp(
            -1j
            * ((mr_angle_at_window_end - phase_h_end) - inspiral_angle_at_window_start)
        )
        head_phasor = inspiral_phasor_shift
        tail_phasor = jnp.ones((), dtype=h_mr.dtype)

    mode_within_window = amp_h * jnp.exp(-1j * phase_h)

    j = jnp.arange(n_insp + n_mr)
    head = jnp.take(h_insp, jnp.clip(j, 0, n_insp - 1)) * head_phasor
    window = jnp.take(mode_within_window, jnp.clip(j - t1_insp, 0, n_win_buf))
    tail = jnp.take(h_mr, jnp.clip(j - t1_insp + t1_mr, 0, n_mr - 1)) * tail_phasor

    hyb = jnp.where(j < t1_insp, head, jnp.where(j <= t1_insp + n_win, window, tail))
    valid_len = t1_insp + n_mr - t1_mr
    hyb = jnp.where(j < valid_len, hyb, 0.0)
    return hyb, valid_len, amp_insp_w, amp_mr_w, amp_h, frq_h, phase_h


def blend_modes_jax(
    inspiral_modes,
    merger_ringdown_modes,
    inspiral_orbital_frequency,
    frq_attach,
    frq_width=10.0,
    delta_t=1.0 / 4096,
    modes_to_blend=[(2, 2), (3, 3), (4, 4)],
    mode_to_align_by=(2, 2),
    blend_using_avg_orbital_frequency=True,
    blend_aligning_merger_to_inspiral=True,
    include_conjugate_modes=True,
    verbose=False,
):
    """Drop-in JAX counterpart of ``blend.blend_modes`` (host wrapper).

    Same signature, validation and return structure (hybrid modes plus window
    indices/diagnostics); computation runs through the jit-safe kernels above
    and results are returned as NumPy arrays cropped to the valid length.
    Only the ``blend_using_avg_orbital_frequency=True`` scheme is supported,
    matching the surrogate IMR path.
    """
    if not blend_using_avg_orbital_frequency:
        raise NotImplementedError(
            "blend_modes_jax supports only blend_using_avg_orbital_frequency=True"
        )
    if frq_width <= 0:
        raise IOError(f"""You are trying to blend over a frequency window of
            negative length (= {frq_width}Hz). Fix this.""")
    if len(inspiral_orbital_frequency) != len(inspiral_modes[mode_to_align_by]):
        raise IOError(
            f"""You asked for hybridization using orbital frequency, but
                the orbital frequency array and inspiral modes array have
                different lengths: {len(inspiral_orbital_frequency)}, {len(inspiral_modes[mode_to_align_by])}"""
        )
    modes_to_blend = list(modes_to_blend)
    if include_conjugate_modes:
        for el, em in modes_to_blend.copy():
            if (el, -em) not in modes_to_blend:
                modes_to_blend.append((el, -em))
    for lm in modes_to_blend:
        if lm not in inspiral_modes or lm not in merger_ringdown_modes:
            raise IOError(
                "We cannot blend {} mode as its missing in the input inspiral modes"
                " ({}) or the merger ringdown modes ({})".format(
                    lm, lm in inspiral_modes, lm in merger_ringdown_modes
                )
            )
    el, em = mode_to_align_by
    if em <= 0:
        raise ValueError(
            f"mode_to_align_by must have positive m, got {mode_to_align_by}"
        )

    h_insp_align = jnp.asarray(inspiral_modes[mode_to_align_by])
    h_mr_align = jnp.asarray(merger_ringdown_modes[mode_to_align_by])
    orb_frq = jnp.asarray(inspiral_orbital_frequency, dtype=jnp.float64)

    # Replicate the out-of-bounds checks of the find_* host wrappers.
    frq_mr_align = compute_frequency_jax(compute_phase_jax(h_mr_align), delta_t)
    for target, series, name in (
        (frq_attach - frq_width / 2, frq_mr_align, "merger-ringdown frequency"),
        (frq_attach + frq_width / 2, frq_mr_align, "merger-ringdown frequency"),
        ((frq_attach + frq_width / 2) / em, orb_frq, "inspiral orbital frequency"),
    ):
        smin, smax = float(jnp.min(series)), float(jnp.max(series))
        if target < smin:
            raise Exception(
                f"Desired frequency {target} out of bounds in {name}, lower than min {smin}"
            )
        if target > smax:
            raise Exception(
                f"Desired frequency {target} out of bounds in {name}, higher than max {smax}"
            )

    t1_insp, t2_insp, t1_mr, t2_mr, found = locate_blend_indices_jax(
        frq_mr_align, orb_frq, frq_attach, frq_width, em
    )
    if not bool(found):
        raise ValueError("Desired frequency bracket not found in timeseries")
    t1_insp, t2_insp = int(t1_insp), int(t2_insp)
    t1_mr, t2_mr = int(t1_mr), int(t2_mr)
    n_win = t2_mr - t1_mr

    if n_win <= 0:
        raise ValueError(
            "Invalid index range for blending. Ensure that t2_index > t1_index for both inspiral and merger-ringdown."
        )
    if t1_insp < 1:
        raise ValueError("""Inspiral too short to accommodate the hybridization window.
Try reducing frq_width or moving frq_attach to a higher value.""")
    if t2_insp + 2 > len(inspiral_modes[mode_to_align_by]):
        raise ValueError(
            """Inspiral too short to accommodate the high attachment frequency.
Try moving frq_attach to a lower value."""
        )
    if t1_mr < 1:
        raise ValueError(
            """Merger-ringdown too short to accommodate the hybridization window.
Try reducing frq_width or moving frq_attach to a higher value."""
        )
    if t2_mr + 2 > len(merger_ringdown_modes[mode_to_align_by]):
        raise ValueError(
            """Merger-ringdown too short to accommodate the high attachment frequency.
Try moving frq_attach to a lower value."""
        )

    frq_insp_window = {}
    frq_mr_window = {}
    amp_insp_window = {}
    amp_mr_window = {}
    amp_hyb_window = {}
    frq_hyb_window = {}
    hybrid_modes = {}

    for el_i, em_i in modes_to_blend:
        h_insp = jnp.asarray(inspiral_modes[(el_i, em_i)])
        h_mr = jnp.asarray(merger_ringdown_modes[(el_i, em_i)])
        if (el_i, em_i) == mode_to_align_by:
            # Alignment mode: MR window frequency from the full-series
            # computation; inspiral window frequency reconstructed locally.
            frq_mr_w = jax.lax.dynamic_slice(frq_mr_align, (t1_mr,), (n_win + 1,))
            frq_insp_w = windowed_mode_frequency_jax(h_insp, delta_t, t1_insp, n_win)
        else:
            frq_insp_w = windowed_mode_frequency_jax(h_insp, delta_t, t1_insp, n_win)
            frq_mr_w = windowed_mode_frequency_jax(h_mr, delta_t, t1_mr, n_win)

        hyb, valid_len, amp_insp_w, amp_mr_w, amp_h, frq_h, phase_h = (
            blend_single_mode_jax(
                h_insp,
                h_mr,
                delta_t,
                t1_insp,
                t1_mr,
                n_win,
                frq_insp_w,
                frq_mr_w,
                n_win,
                blend_aligning_merger_to_inspiral,
            )
        )
        hybrid_modes[(el_i, em_i)] = np.asarray(hyb)[: int(valid_len)]
        frq_insp_window[(el_i, em_i)] = np.asarray(frq_insp_w)
        frq_mr_window[(el_i, em_i)] = np.asarray(frq_mr_w)
        amp_insp_window[(el_i, em_i)] = np.asarray(amp_insp_w)
        amp_mr_window[(el_i, em_i)] = np.asarray(amp_mr_w)
        amp_hyb_window[(el_i, em_i)] = np.asarray(amp_h)
        frq_hyb_window[(el_i, em_i)] = np.asarray(frq_h)

    return (
        hybrid_modes,
        t1_insp,
        t1_mr,
        t2_insp,
        t2_mr,
        frq_insp_window,
        frq_hyb_window,
        frq_mr_window,
        amp_insp_window,
        amp_hyb_window,
        amp_mr_window,
    )
