"""Tier-A parity tests: esigmapy.blend_jax vs esigmapy.blend on identical inputs.

The JAX port mirrors blend.py operation-for-operation; the only expected
deviations come from XLA reassociating the phase cumsums, which shows up at
~1e-12 rad on unwrapped phases (and correspondingly ~1e-10 relative on
frequencies derived from them). Waveform-level outputs must agree to 1e-12 of
the mode amplitude scale; window indices must agree exactly.
"""

import os

import numpy as np
import pytest

pytest.importorskip("jax")

try:
    from esigmapy import blend_jax
except Exception as exc:  # pragma: no cover - environment dependent
    pytest.skip(f"blend_jax unavailable: {exc}", allow_module_level=True)

import jax.numpy as jnp

from esigmapy import blend

_DATA = os.environ.get("ESIGMASUR_DATA_PATH")
if not _DATA or not os.path.isdir(os.path.join(_DATA, "ecc_sur_data")):
    pytest.skip("ESIGMASUR_DATA_PATH not configured", allow_module_level=True)

_DT = 0.000244140625


def _rel(a, b):
    a = np.asarray(a)
    b = np.asarray(b)
    return np.abs(a - b) / np.maximum(np.abs(b), 1e-30)


def _synthetic_pair():
    """Synthetic inspiral/merger-ringdown chirps with overlapping frequency
    bands, exercising the blend independently of surrogate/lalsim data."""
    n_insp, n_mr = 6000, 3000
    t_i = np.arange(n_insp) * _DT
    # inspiral (2,2) frequency 20 -> 80 Hz, mild amplitude growth
    f_i = 20.0 + 60.0 * (t_i / t_i[-1]) ** 2
    phi_i = 2 * np.pi * np.cumsum(f_i) * _DT
    h_i = (1.0 + 0.3 * (t_i / t_i[-1])) * np.exp(-1j * phi_i)
    t_m = np.arange(n_mr) * _DT
    # merger-ringdown frequency 40 -> 300 Hz, decaying amplitude tail
    f_m = 40.0 + 260.0 * (t_m / t_m[-1]) ** 3
    phi_m = 2 * np.pi * np.cumsum(f_m) * _DT + 0.7
    a_m = 1.2 * np.exp(-3.0 * (t_m / t_m[-1]) ** 4)
    h_m = a_m * np.exp(-1j * phi_m)
    insp = {(2, 2): h_i, (2, -2): np.conj(h_i)}
    mr = {(2, 2): h_m, (2, -2): np.conj(h_m)}
    orb_frq = f_i / 2.0
    return insp, mr, orb_frq


@pytest.fixture(scope="module")
def real_inputs():
    """Inspiral surrogate modes + lalsim NRSur7dq4 MR modes with the same
    transition frequency / window the NumPy IMR driver would compute."""
    import lal

    lalsim = pytest.importorskip("lalsimulation")  # noqa: F841
    from esigmapy.generator import _get_transition_frequency_window
    from esigmapy.mr_generator import get_mr_modes
    from esigmapy.surrogate.generator import get_inspiral_esigmasur_modes
    from esigmapy.utils import f_ISCO_spin

    m1, m2 = 41.4, 18.6
    t, orb, insp = get_inspiral_esigmasur_modes(
        mass1=m1,
        mass2=m2,
        delta_t=_DT,
        reference_eccentricity=0.25,
        reference_mean_anomaly=1.3,
        include_conjugate_modes=True,
        return_orbital_params=["x"],
        return_pycbc_timeseries=False,
    )
    orb_frq = orb["x"] ** 1.5 / ((m1 + m2) * lal.MTSUN_SI) / (2 * np.pi)
    f_Kerr = f_ISCO_spin(m1, m2, 0.0, 0.0)
    f_Schwarz = 6.0**-1.5 / (m1 + m2) / lal.MTSUN_SI / lal.PI
    f_tr = min(f_Kerr, f_Schwarz)
    f_win = (
        _get_transition_frequency_window(
            None, orb_frq, _DT, f_tr / 2, 0.25, False, True, True
        )
        * 2
    )
    f_tr -= f_win / 2
    f_lower_mr = (f_tr - f_win / 2) * 0.9
    try:
        mr = get_mr_modes(
            mass1=m1,
            mass2=m2,
            f_lower=f_lower_mr,
            f_ref=f_lower_mr,
            delta_t=_DT,
            modes_to_use=[(2, 2), (2, -2)],
            approximant="NRSur7dq4",
        )
    except Exception as exc:  # pragma: no cover - environment dependent
        pytest.skip(f"lalsim NRSur7dq4 unavailable: {exc}")
    return dict(insp=insp, mr=mr, orb_frq=orb_frq, f_attach=f_tr, f_width=f_win)


# --- phase / frequency primitives -------------------------------------------


def _phase_freq_cases(real_inputs):
    insp, mr = real_inputs["insp"], real_inputs["mr"]
    t = np.arange(4096) * _DT
    chirp = np.exp(-1j * 2 * np.pi * (30.0 * t + 400.0 * t**2))
    return {
        "insp22": insp[(2, 2)],
        "insp2m2": insp[(2, -2)],
        "mr22": mr[(2, 2)],
        "chirp": chirp,
        "constant": np.full(64, 0.3 - 0.4j),
        "short": chirp[:2],
        "single": chirp[:1],
    }


def test_compute_phase_parity(real_inputs):
    for name, w in _phase_freq_cases(real_inputs).items():
        pn = blend.compute_phase(np.asarray(w))
        pj = np.asarray(blend_jax.compute_phase_jax(jnp.asarray(w)))
        assert pj.dtype == np.float64
        assert np.abs(pj - pn).max() < 1e-9, name


def test_compute_frequency_parity(real_inputs):
    for name, w in _phase_freq_cases(real_inputs).items():
        pn = blend.compute_phase(np.asarray(w))
        fn = blend.compute_frequency(pn, _DT)
        fj = np.asarray(
            blend_jax.compute_frequency_jax(
                blend_jax.compute_phase_jax(jnp.asarray(w)), _DT
            )
        )
        assert fj.dtype == np.float64
        # zero-frequency entries (constant case) guard the denominator
        assert _rel(fj, fn).max() < 1e-8 or np.abs(fj - fn).max() < 1e-8, name


# --- bracket searches --------------------------------------------------------


def _sweep_targets(series):
    lo, hi = np.min(series), np.max(series)
    span = hi - lo
    targets = list(np.linspace(lo + 1e-3 * span, hi - 1e-3 * span, 37))
    targets += list(np.asarray(series)[[5, len(series) // 2, -5]])  # exact hits
    return targets


def test_find_first_parity(real_inputs):
    frq_mr = blend.compute_frequency(
        blend.compute_phase(real_inputs["mr"][(2, 2)]), _DT
    )
    fj = jnp.asarray(frq_mr)
    for target in _sweep_targets(frq_mr):
        idx_j, found = blend_jax.find_first_index_jax(fj, target)
        try:
            idx_n = blend.find_first_value_location_in_series(frq_mr, target)
        except ValueError:
            assert not bool(found), f"target={target}"
            continue
        assert bool(found), f"target={target}"
        assert int(idx_j) == idx_n, f"target={target}"


def test_find_last_parity(real_inputs):
    orb_frq = real_inputs["orb_frq"]
    fj = jnp.asarray(orb_frq)
    for target in _sweep_targets(orb_frq):
        idx_j, found = blend_jax.find_last_index_jax(fj, target)
        try:
            idx_n = blend.find_last_value_location_in_series(orb_frq, target)
        except ValueError:
            assert not bool(found), f"target={target}"
            continue
        assert bool(found), f"target={target}"
        assert int(idx_j) == idx_n, f"target={target}"


def test_find_first_tie_break():
    # bracket [1.0, 3.0] around 2.0: equidistant -> left endpoint (numba `<=`)
    series = jnp.asarray([0.0, 1.0, 3.0, 4.0])
    idx, found = blend_jax.find_first_index_jax(series, 2.0)
    assert bool(found) and int(idx) == 1
    ref = blend.find_first_value_location_in_series(np.asarray(series), 2.0)
    assert int(idx) == ref


# --- full blend --------------------------------------------------------------


def _assert_blend_parity(insp, mr, orb_frq, f_attach, f_width, align, conj):
    kw = dict(
        frq_width=f_width,
        delta_t=_DT,
        mode_to_align_by=(2, 2),
        blend_using_avg_orbital_frequency=True,
        blend_aligning_merger_to_inspiral=align,
        include_conjugate_modes=conj,
    )
    # blend.blend_modes mutates modes_to_blend: pass fresh lists
    ret_n = blend.blend_modes(
        {k: np.asarray(v) for k, v in insp.items()},
        {k: np.asarray(v) for k, v in mr.items()},
        np.asarray(orb_frq),
        f_attach,
        modes_to_blend=[(2, 2)],
        **kw,
    )
    ret_j = blend_jax.blend_modes_jax(
        insp, mr, orb_frq, f_attach, modes_to_blend=[(2, 2)], **kw
    )

    # indices (positions 1-4) must match exactly
    assert tuple(ret_j[1:5]) == tuple(int(v) for v in ret_n[1:5])

    hyb_n, hyb_j = ret_n[0], ret_j[0]
    assert set(hyb_j) == set(hyb_n)
    for lm in hyb_n:
        peak = np.abs(hyb_n[lm]).max()
        assert hyb_j[lm].shape == hyb_n[lm].shape, lm
        assert np.abs(hyb_j[lm] - hyb_n[lm]).max() < 1e-12 * peak, lm

    # window diagnostics: frequencies to 1e-8 relative (phase-cumsum
    # reassociation), amplitudes to 1e-12 relative
    for pos, tol in (
        (5, 1e-8),
        (6, 1e-8),
        (7, 1e-8),
        (8, 1e-12),
        (9, 1e-12),
        (10, 1e-12),
    ):
        dn, dj = ret_n[pos], ret_j[pos]
        assert set(dj) == set(dn)
        for lm in dn:
            assert _rel(dj[lm], dn[lm]).max() < tol, (pos, lm)


@pytest.mark.parametrize("align", [True, False])
@pytest.mark.parametrize("conj", [True, False])
def test_blend_modes_parity_real(real_inputs, align, conj):
    _assert_blend_parity(
        real_inputs["insp"],
        real_inputs["mr"],
        real_inputs["orb_frq"],
        real_inputs["f_attach"],
        real_inputs["f_width"],
        align,
        conj,
    )


@pytest.mark.parametrize("align", [True, False])
def test_blend_modes_parity_synthetic(align):
    insp, mr, orb_frq = _synthetic_pair()
    _assert_blend_parity(insp, mr, orb_frq, 60.0, 20.0, align, True)


def test_blend_modes_parity_wide_window(real_inputs):
    # user-specified (wider) frequency window instead of the driver's default;
    # sized so the lower window edge stays inside the MR frequency band
    frq_mr = blend.compute_frequency(
        blend.compute_phase(real_inputs["mr"][(2, 2)]), _DT
    )
    f_attach = real_inputs["f_attach"]
    width_max = 2.0 * (f_attach - frq_mr.min())
    width = min(3.0 * real_inputs["f_width"], 0.9 * width_max)
    assert width > real_inputs["f_width"]
    _assert_blend_parity(
        real_inputs["insp"],
        real_inputs["mr"],
        real_inputs["orb_frq"],
        f_attach,
        width,
        True,
        True,
    )


def test_blend_modes_error_paths(real_inputs):
    insp, mr, orb_frq = (
        real_inputs["insp"],
        real_inputs["mr"],
        real_inputs["orb_frq"],
    )
    with pytest.raises(IOError):
        blend_jax.blend_modes_jax(
            insp,
            mr,
            orb_frq,
            real_inputs["f_attach"],
            frq_width=-1.0,
            modes_to_blend=[(2, 2)],
        )
    with pytest.raises(Exception, match="out of bounds"):
        blend_jax.blend_modes_jax(
            insp, mr, orb_frq, 1e6, frq_width=10.0, modes_to_blend=[(2, 2)]
        )
    with pytest.raises(ValueError, match="positive m"):
        blend_jax.blend_modes_jax(
            insp,
            mr,
            orb_frq,
            real_inputs["f_attach"],
            frq_width=real_inputs["f_width"],
            modes_to_blend=[(2, -2)],
            mode_to_align_by=(2, -2),
            include_conjugate_modes=False,
        )
