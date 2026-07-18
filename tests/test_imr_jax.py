"""Parity tests for the JAX IMR pipeline (esigmapy/surrogate/generator_jax.py)
against the NumPy backend.

Tiered accuracy scheme (see AGENTS.md):
- Tier A (this file + test_blend_jax.py): component ports vs NumPy on
  identical inputs at ~1e-12/exact-index level.
- Tier B: full pipeline with lalsim-generated MR modes injected, sample-level.
- Tier C: gwsurrogate-JAX vs lalsim NRSur7dq4 characterization.
- Tier D: full JAX-vs-NumPy IMR at mismatch/phase level.
"""

import os

import numpy as np
import pytest

pytest.importorskip("jax")

try:
    from esigmapy.surrogate import generator_jax
except Exception as exc:  # pragma: no cover - environment dependent
    pytest.skip(f"JAX IMR backend unavailable: {exc}", allow_module_level=True)

import jax.numpy as jnp
import lal

from esigmapy.generator import _get_transition_frequency_window
from esigmapy.surrogate.generator import get_inspiral_esigmasur_modes

_DATA = os.environ.get("ESIGMASUR_DATA_PATH")
if not _DATA or not os.path.isdir(os.path.join(_DATA, "ecc_sur_data")):
    pytest.skip("ESIGMASUR_DATA_PATH not configured", allow_module_level=True)

_DT = 0.000244140625


def _orbital_frequency(mass1, mass2, e_ref, l_ref, t_start=None):
    _, orb, _ = get_inspiral_esigmasur_modes(
        mass1=mass1,
        mass2=mass2,
        delta_t=_DT,
        reference_eccentricity=e_ref,
        reference_mean_anomaly=l_ref,
        t_start=t_start,
        return_orbital_params=["x"],
        return_pycbc_timeseries=False,
    )
    return orb["x"] ** 1.5 / ((mass1 + mass2) * lal.MTSUN_SI) / (2 * np.pi)


@pytest.fixture(scope="module")
def orbital_freqs():
    return {
        "ecc_high": _orbital_frequency(41.4, 18.6, 0.43, 1.3),
        "ecc_mid": _orbital_frequency(24.0, 16.0, 0.25, 0.0),
        "circ": _orbital_frequency(12.0, 8.0, 1e-5, 0.0),
    }


# --- Tier A: transition frequency window -------------------------------------


@pytest.mark.parametrize("case", ["ecc_high", "ecc_mid", "circ"])
@pytest.mark.parametrize("num_hyb_orbits", [0.25, 1.0, 3.0])
@pytest.mark.parametrize("centered", [False, True])
def test_transition_window_parity(orbital_freqs, case, num_hyb_orbits, centered):
    orb_frq = orbital_freqs[case]
    # orbital transition frequency near (but below) the end-of-inspiral value
    for f_frac in (0.98, 0.7):
        f_tr = f_frac * orb_frq[-1]
        f_win_np = _get_transition_frequency_window(
            None, orb_frq, _DT, f_tr, num_hyb_orbits, centered, True, failsafe=True
        )
        f_win_j, start_found, end_found = generator_jax.transition_frequency_window_jax(
            jnp.asarray(orb_frq), _DT, f_tr, num_hyb_orbits, centered
        )
        assert float(f_win_j) == pytest.approx(f_win_np, rel=1e-10, abs=1e-14), (
            case,
            f_frac,
        )


@pytest.mark.parametrize("case", ["ecc_high", "circ"])
@pytest.mark.parametrize("centered", [False, True])
def test_transition_window_failsafe_clamp(orbital_freqs, case, centered):
    """A window of (absurdly) many orbits exhausts the inspiral: NumPy's
    failsafe clamps the start index to 0 and non-failsafe raises; the JAX port
    clamps identically and reports start_found=False."""
    orb_frq = orbital_freqs[case]
    f_tr = 0.98 * orb_frq[-1]
    num_hyb_orbits = 1.0e6
    f_win_np = _get_transition_frequency_window(
        None, orb_frq, _DT, f_tr, num_hyb_orbits, centered, True, failsafe=True
    )
    with pytest.raises(RuntimeError):
        _get_transition_frequency_window(
            None, orb_frq, _DT, f_tr, num_hyb_orbits, centered, True, failsafe=False
        )
    f_win_j, start_found, end_found = generator_jax.transition_frequency_window_jax(
        jnp.asarray(orb_frq), _DT, f_tr, num_hyb_orbits, centered
    )
    assert not bool(start_found)
    assert float(f_win_j) == pytest.approx(f_win_np, rel=1e-10, abs=1e-14)


def test_transition_window_transition_above_max(orbital_freqs):
    """f_mr_transition above the maximum orbital frequency: the reverse argmax
    lands on the last sample in both backends."""
    orb_frq = orbital_freqs["ecc_mid"]
    f_tr = 2.0 * orb_frq.max()
    f_win_np = _get_transition_frequency_window(
        None, orb_frq, _DT, f_tr, 0.25, False, True, failsafe=True
    )
    f_win_j, _, _ = generator_jax.transition_frequency_window_jax(
        jnp.asarray(orb_frq), _DT, f_tr, 0.25, False
    )
    assert float(f_win_j) == pytest.approx(f_win_np, rel=1e-10)


def test_transition_window_jit_traceable(orbital_freqs):
    """The window computation must jit cleanly with traced frequency inputs."""
    import jax

    orb_frq = jnp.asarray(orbital_freqs["ecc_high"])

    fn = jax.jit(
        lambda f_tr: generator_jax.transition_frequency_window_jax(
            orb_frq, _DT, f_tr, 0.25, False
        )[0]
    )
    f_tr = 0.9 * float(orbital_freqs["ecc_high"][-1])
    expected = _get_transition_frequency_window(
        None, orbital_freqs["ecc_high"], _DT, f_tr, 0.25, False, True, failsafe=True
    )
    assert float(fn(f_tr)) == pytest.approx(expected, rel=1e-10)


# --- Tier C: gwsurrogate-JAX NRSur7dq4 vs lalsim NRSur7dq4 -------------------


def _subsample_peak_time(h, dt):
    a = np.abs(h)
    i = int(np.argmax(a))
    y0, y1, y2 = a[i - 1], a[i], a[i + 1]
    return (i + 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)) * dt


@pytest.fixture(scope="module")
def nrsur_mr():
    try:
        return generator_jax.NRSurMergerRingdownJAX()
    except Exception as exc:  # pragma: no cover - environment dependent
        pytest.skip(f"JAX NRSur7dq4 unavailable: {exc}")


@pytest.mark.parametrize("m1,m2", [(41.4, 18.6), (30.0, 30.0)])
def test_nrsur_jax_vs_lalsim_characterization(nrsur_mr, m1, m2):
    """The two NRSur7dq4 implementations (gwsurrogate-JAX vs lalsimulation)
    read the same model but differ in codebase and grid conditioning; after
    sub-sample peak alignment and constant phase-offset removal they agree to
    ~2e-4 in amplitude (of peak) and ~1.5e-3 rad in phase over the usable
    band (measured 2026-07; bounds hold 2.5x margin). This floor drives the
    Tier D full-IMR tolerances."""
    from esigmapy import blend
    from esigmapy.mr_generator import get_mr_modes

    M, q = m1 + m2, max(m1 / m2, m2 / m1)
    n_mr = nrsur_mr.num_samples(M, _DT)
    h_j = np.asarray(nrsur_mr.modes_2pm2(q, M, _DT, 1.0, n_mr)[(2, 2)])
    try:
        h_n = np.asarray(
            get_mr_modes(
                mass1=m1,
                mass2=m2,
                f_lower=30.0,
                f_ref=30.0,
                delta_t=_DT,
                modes_to_use=[(2, 2), (2, -2)],
                approximant="NRSur7dq4",
            )[(2, 2)]
        )
    except Exception as exc:  # pragma: no cover - environment dependent
        pytest.skip(f"lalsim NRSur7dq4 unavailable: {exc}")

    assert np.abs(h_j).max() == pytest.approx(np.abs(h_n).max(), rel=1e-3)

    t_j = np.arange(len(h_j)) * _DT - _subsample_peak_time(h_j, _DT)
    t_n = np.arange(len(h_n)) * _DT - _subsample_peak_time(h_n, _DT)
    ph_j = blend.compute_phase(h_j)
    ph_n = blend.compute_phase(h_n)

    tau_lo = max(t_j[0], t_n[0]) + 20 * _DT
    tau_hi = min(t_j[-1], t_n[-1], 0.02)
    tau = np.arange(tau_lo, tau_hi, _DT)
    amp_j = np.interp(tau, t_j, np.abs(h_j))
    amp_n = np.interp(tau, t_n, np.abs(h_n))
    phi_diff = np.interp(tau, t_j, ph_j) - np.interp(tau, t_n, ph_n)
    phi_diff -= phi_diff[len(phi_diff) // 2]

    assert np.abs(amp_j - amp_n).max() / amp_n.max() < 5e-4
    assert np.abs(phi_diff).max() < 5e-3


def test_nrsur_jax_mode_symmetry(nrsur_mr):
    """Non-spinning: (2,-2) ~ conj((2,2)), but only to ~3e-4 of peak -- the
    NRSur7dq4 fits carry small mode asymmetries even at zero spin (measured
    3.4e-4 on 2026-07). This is why the pipeline blends (2,-2) as its own
    mode (as the NumPy backend does with lalsim's) instead of conjugating."""
    modes = nrsur_mr.modes_2pm2(2.226, 60.0, _DT, 1.0, 4096)
    h22 = np.asarray(modes[(2, 2)])
    h2m2 = np.asarray(modes[(2, -2)])
    peak = np.abs(h22).max()
    assert np.abs(h2m2 - np.conj(h22)).max() < 1e-3 * peak


# --- Tiers B & D: full IMR pipeline ------------------------------------------

_CIRC_DIR = os.path.join(_DATA, "circ_sur_data")
_ECC_DIR = os.path.join(_DATA, "ecc_sur_data")

# (m1, m2, e_ref, l_ref, t_start): short inspirals for test speed. The Tier B
# configs were verified to sit away from window-index razor edges (the JAX
# inspiral differs from NumPy's at ~1e-9 relative, which at unlucky parameter
# points shifts a frequency-bracket index by one sample -- sometimes with
# compensating shifts that keep the total length equal while displacing the
# window). Such points fail strict sample parity legitimately; e.g.
# (36, 24, 0.1, 0.0) and (36, 24, 0.18, 0.7) are razor-edge cases.
_IMR_CONFIGS = [
    (41.4, 18.6, 0.25, 1.3, -40.0),
    (25.0, 5.0, 0.4, 2.0, -30.0),
    (30.0, 10.0, 0.2, 4.0, -40.0),
]


@pytest.fixture(scope="module")
def imr_jax():
    try:
        return generator_jax.IMRESIGMASurJAX(_ECC_DIR, _CIRC_DIR)
    except Exception as exc:  # pragma: no cover - environment dependent
        pytest.skip(f"IMRESIGMASurJAX unavailable: {exc}")


def _numpy_imr(m1, m2, e0, l0, t_start, **kw):
    from esigmapy.surrogate.generator import get_imr_esigmasur_mode

    return get_imr_esigmasur_mode(
        mass1=m1,
        mass2=m2,
        delta_t=_DT,
        reference_eccentricity=e0,
        reference_mean_anomaly=l0,
        t_start=t_start,
        include_conjugate_modes=True,
        **kw,
    )


def _lalsim_mr_for_driver(m1, m2, e0, l0, t_start):
    """The lalsim NRSur7dq4 modes the NumPy driver generates internally
    (first retry iteration; replicates its f_lower recipe)."""
    from esigmapy.mr_generator import get_mr_modes

    M = m1 + m2
    orb_frq = _orbital_frequency(m1, m2, e0, l0, t_start=t_start)
    f_tr0 = generator_jax.IMRESIGMASurJAX._f_mr_transition_default(M, m1 / m2)
    f_win = (
        _get_transition_frequency_window(
            None, orb_frq, _DT, f_tr0 / 2, 0.25, False, True, True
        )
        * 2
    )
    f_tr = f_tr0 - f_win / 2
    f_lower = (f_tr - f_win / 2) * 0.9
    modes = get_mr_modes(
        mass1=m1,
        mass2=m2,
        f_lower=f_lower,
        f_ref=f_lower,
        delta_t=_DT,
        modes_to_use=[(2, 2), (2, -2)],
        approximant="NRSur7dq4",
    )
    return {k: np.asarray(v) for k, v in modes.items()}


@pytest.mark.parametrize("m1,m2,e0,l0,t_start", _IMR_CONFIGS)
def test_tier_b_injected_mr_parity(imr_jax, m1, m2, e0, l0, t_start):
    """Full JAX driver with the NumPy driver's own lalsim MR modes injected:
    everything except the MR backend must match at the inspiral-backend
    equivalence level (~1e-9 relative; bound 1e-8 of peak).

    The window bracket indices come from the JAX inspiral (1e-9-different from
    NumPy's), so a razor-edge parameter point could shift a window index by
    one sample and legitimately fail the strict comparison; the configs here
    were checked not to sit on such edges."""
    try:
        mr_lal = _lalsim_mr_for_driver(m1, m2, e0, l0, t_start)
    except Exception as exc:  # pragma: no cover - environment dependent
        pytest.skip(f"lalsim NRSur7dq4 unavailable: {exc}")
    modes_np = _numpy_imr(m1, m2, e0, l0, t_start)
    t_j, modes_j = imr_jax(
        m1 + m2,
        (m1 / m2, e0, l0),
        _DT,
        t_start=t_start,
        include_conjugate_modes=True,
        merger_ringdown_modes=mr_lal,
    )
    for lm in ((2, 2), (2, -2)):
        h_np = np.asarray(modes_np[lm].data)
        h_j = np.asarray(modes_j[lm])
        assert len(h_j) == len(h_np), lm
        assert float(modes_np[lm].sample_times[0]) == pytest.approx(t_j[0], abs=1e-9)
        peak = np.abs(h_np).max()
        assert np.abs(h_j - h_np).max() < 1e-8 * peak, lm


def _compare_tier_d(modes_np, t_j, modes_j, remove_const_phase=False):
    """Tier D bounds (measured 2026-07 on full-length configs: inspiral phase
    ~1e-11 rad; to-peak phase <=0.15 rad and amp <=6e-3 of peak, dominated by
    discrete window snapping of the sub-sample-offset MR grids; overall amp
    <=3e-2 of peak; optimized match 1-M <= 2.4e-5). Locked with ~2.5x margin."""
    from esigmapy import blend

    for lm in ((2, 2), (2, -2)):
        h_np = np.asarray(modes_np[lm].data)
        h_j = np.asarray(modes_j[lm])
        assert float(modes_np[lm].sample_times[0]) == pytest.approx(t_j[0], abs=1e-9)
        n = min(len(h_np), len(h_j))
        peak = np.abs(h_np).max()
        a_n, a_j = np.abs(h_np[:n]), np.abs(h_j[:n])
        ph_n = blend.compute_phase(h_np[:n])
        ph_j = blend.compute_phase(h_j[:n])
        dph = ph_j - ph_n
        if remove_const_phase:
            dph -= np.median(dph[: n // 4])
        ipk = int(np.argmax(a_n))
        strong = a_n > 1e-2 * peak
        assert np.abs(dph[: ipk // 2]).max() < 1e-6, lm
        assert np.abs(dph[:ipk]).max() < 0.4, lm
        assert np.abs(dph[strong]).max() < 0.7, lm
        assert np.abs(a_j[:ipk] - a_n[:ipk]).max() < 0.02 * peak, lm
        assert np.abs(a_j - a_n).max() < 0.06 * peak, lm

    # time/phase-optimized match on the real part of (2,2)
    try:
        from pycbc.filter import match
        from pycbc.types import TimeSeries
    except Exception:  # pragma: no cover - environment dependent
        return
    h_np = np.asarray(modes_np[(2, 2)].data)
    h_j = np.asarray(modes_j[(2, 2)])
    ln = 2 ** int(np.ceil(np.log2(max(len(h_np), len(h_j)))))
    tsn = TimeSeries(np.real(np.pad(h_np, (0, ln - len(h_np)))), delta_t=_DT)
    tsj = TimeSeries(np.real(np.pad(h_j, (0, ln - len(h_j)))), delta_t=_DT)
    mm, _ = match(tsn, tsj)
    # short test inspirals weight the merger more than the full-length
    # measurement (1-M <= 2.4e-5); observed up to ~1.6e-4 here
    assert 1 - mm < 5e-4


@pytest.mark.parametrize("m1,m2,e0,l0,t_start", _IMR_CONFIGS)
def test_tier_d_full_imr(imr_jax, m1, m2, e0, l0, t_start):
    modes_np = _numpy_imr(m1, m2, e0, l0, t_start)
    t_j, modes_j = imr_jax(
        m1 + m2,
        (m1 / m2, e0, l0),
        _DT,
        t_start=t_start,
        include_conjugate_modes=True,
    )
    _compare_tier_d(modes_np, t_j, modes_j)


def test_tier_d_align_inspiral_to_merger(imr_jax):
    """blend_aligning_merger_to_inspiral=False: the overall phase is anchored
    to the MR piece, whose raw phase convention differs between the two
    NRSur7dq4 implementations, so compare up to a constant phase offset."""
    m1, m2, e0, l0, t_start = _IMR_CONFIGS[0]
    modes_np = _numpy_imr(
        m1,
        m2,
        e0,
        l0,
        t_start,
        blend_aligning_merger_to_inspiral=False,
        coa_phase=0.3,
    )
    t_j, modes_j = imr_jax(
        m1 + m2,
        (m1 / m2, e0, l0),
        _DT,
        t_start=t_start,
        include_conjugate_modes=True,
        blend_aligning_merger_to_inspiral=False,
        coa_phase=0.3,
    )
    _compare_tier_d(modes_np, t_j, modes_j, remove_const_phase=True)


@pytest.mark.parametrize(
    "kw",
    [
        dict(num_hyb_orbits=1.0),
        # centered windows need a transition below the ISCO default: the
        # window's upper edge must stay within the inspiral's frequency reach
        dict(keep_f_mr_transition_at_center=True, f_mr_transition=62.0),
        dict(f_window_mr_transition=20.0),
        dict(f_mr_transition=1.0e4),  # failsafe clamp engages in both backends
    ],
)
def test_tier_d_option_variants(imr_jax, kw):
    m1, m2, e0, l0, t_start = _IMR_CONFIGS[0]
    modes_np = _numpy_imr(m1, m2, e0, l0, t_start, **kw)
    t_j, modes_j = imr_jax(
        m1 + m2,
        (m1 / m2, e0, l0),
        _DT,
        t_start=t_start,
        include_conjugate_modes=True,
        **kw,
    )
    _compare_tier_d(modes_np, t_j, modes_j)


def test_return_structures(imr_jax):
    """Orbital-variable and hybridization-info returns mirror the NumPy
    backend's content."""
    m1, m2, e0, l0, t_start = _IMR_CONFIGS[0]
    ret_np = _numpy_imr(
        m1,
        m2,
        e0,
        l0,
        t_start,
        return_hybridization_info=True,
        return_orbital_params=True,
    )
    hyb_np, orb_np, modes_np = ret_np
    t_j, hyb_j, orb_j, modes_j = imr_jax(
        m1 + m2,
        (m1 / m2, e0, l0),
        _DT,
        t_start=t_start,
        include_conjugate_modes=True,
        return_hybridization_info=True,
        return_orbital_params=True,
    )
    assert set(orb_j) == {"e", "l", "x"}
    for key in ("e", "x"):
        v_np = np.asarray(orb_np[key].data)
        v_j = np.asarray(orb_j[key])
        assert v_j.shape == v_np.shape
        assert np.abs(v_j - v_np).max() < 1e-6 * max(1.0, np.abs(v_np).max())
    # hybridization info: inspiral-grid window indices within one sample of
    # NumPy's; the MR-grid indices live on different grids (full-span NRSur
    # vs lalsim-from-f_lower), so only the window length is comparable
    for pos in (1, 3):  # t1_index_insp, t2_index_insp
        assert abs(int(hyb_j[pos]) - int(hyb_np[pos])) <= 1
    n_win_j = int(hyb_j[4]) - int(hyb_j[2])
    n_win_np = int(hyb_np[4]) - int(hyb_np[2])
    assert abs(n_win_j - n_win_np) <= 1
    assert set(hyb_j[0]) == set(modes_j)


def test_circular_fallback_unsupported(imr_jax):
    with pytest.raises(NotImplementedError):
        imr_jax(60.0, (2.0, 0.0, 0.0), _DT, t_start=-10.0)


def test_coa_phase_required_for_inspiral_alignment(imr_jax):
    with pytest.raises(IOError):
        imr_jax(
            60.0,
            (2.0, 0.2, 0.0),
            _DT,
            t_start=-10.0,
            blend_aligning_merger_to_inspiral=False,
        )


# --- Tier E: traced evaluator (jit / vmap / grad / GPU) ----------------------


_EVAL_M, _EVAL_TS = 60.0, -40.0
_EVAL_PARAMS = (41.4 / 18.6, 0.25, 1.3)


@pytest.fixture(scope="module")
def imr_evaluator(imr_jax):
    times, fn = imr_jax.imr_parameter_space_evaluator(_EVAL_M, _DT, t_start=_EVAL_TS)
    return times, fn


@pytest.mark.parametrize("params", [_EVAL_PARAMS, (3.0, 0.2, 4.0)])
def test_evaluator_matches_class_path(imr_jax, imr_evaluator, params):
    """The traced evaluator reproduces the host __call__ path to round-off
    (measured 2.1e-13 of peak; same kernels, different orchestration)."""
    times, fn = imr_evaluator
    t_c, modes_c = imr_jax(_EVAL_M, params, _DT, t_start=_EVAL_TS)
    h_c = np.asarray(modes_c[(2, 2)])
    h_e, vlen = fn(*params)
    h_e, vlen = np.asarray(h_e), int(vlen)
    assert h_e.dtype == np.complex128
    assert vlen == len(h_c)
    assert times[0] == pytest.approx(t_c[0], abs=1e-12)
    peak = np.abs(h_c).max()
    assert np.abs(h_e[:vlen] - h_c).max() < 1e-11 * peak
    assert np.abs(h_e[vlen:]).max() == 0.0


def test_evaluator_options_match_class(imr_jax):
    q, e0, l0 = _EVAL_PARAMS
    kw = dict(
        blend_aligning_merger_to_inspiral=False,
        coa_phase=0.7,
        f_window_mr_transition=15.0,
    )
    t_c, modes_c = imr_jax(_EVAL_M, (q, e0, l0), _DT, t_start=_EVAL_TS, **kw)
    h_c = np.asarray(modes_c[(2, 2)])
    _, fn = imr_jax.imr_parameter_space_evaluator(_EVAL_M, _DT, t_start=_EVAL_TS, **kw)
    h_e, vlen = fn(q, e0, l0)
    h_e, vlen = np.asarray(h_e), int(vlen)
    assert vlen == len(h_c)
    assert np.abs(h_e[:vlen] - h_c).max() < 1e-11 * np.abs(h_c).max()


def test_evaluator_amp_phase_output(imr_jax):
    q, e0, l0 = _EVAL_PARAMS
    _, fn_m = imr_jax.imr_parameter_space_evaluator(_EVAL_M, _DT, t_start=_EVAL_TS)
    _, fn_ap = imr_jax.imr_parameter_space_evaluator(
        _EVAL_M, _DT, t_start=_EVAL_TS, output="amp_phase"
    )
    h, vlen = fn_m(q, e0, l0)
    amp, phase, vlen2 = fn_ap(q, e0, l0)
    assert int(vlen) == int(vlen2)
    assert np.asarray(amp).dtype == np.float64
    assert np.asarray(phase).dtype == np.float64
    n = int(vlen)
    h_rec = np.asarray(amp)[:n] * np.exp(-1j * np.asarray(phase)[:n])
    h = np.asarray(h)[:n]
    assert np.abs(h_rec - h).max() < 1e-10 * np.abs(h).max()
    assert np.asarray(amp)[n:].max() == 0.0


def test_evaluator_jit_single_compile(imr_evaluator):
    import jax

    _, fn = imr_evaluator
    jfn = jax.jit(fn)
    q, e0, l0 = _EVAL_PARAMS
    h1, v1 = jfn(q, e0, l0)
    h2, v2 = jfn(q * 1.1, e0 + 0.02, l0 - 0.5)  # new values, same shapes
    jax.block_until_ready(h2)
    assert jfn._cache_size() == 1
    # jit result equals eager result
    h_e, v_e = fn(q, e0, l0)
    assert int(v1) == int(v_e)
    peak = np.abs(np.asarray(h_e)).max()
    assert np.abs(np.asarray(h1) - np.asarray(h_e)).max() < 1e-10 * peak


def test_evaluator_vmap_matches_loop(imr_evaluator):
    import jax
    import jax.numpy as jnp

    _, fn = imr_evaluator
    qs = jnp.array([2.0, 41.4 / 18.6, 3.0])
    es = jnp.array([0.2, 0.25, 0.3])
    ls = jnp.array([0.0, 1.3, 2.0])
    hb, vb = jax.vmap(fn)(qs, es, ls)
    for k in range(3):
        h_k, v_k = fn(qs[k], es[k], ls[k])
        assert int(vb[k]) == int(v_k)
        peak = np.abs(np.asarray(h_k)).max()
        assert np.abs(np.asarray(hb[k]) - np.asarray(h_k)).max() < 1e-10 * peak


def test_evaluator_grad_matches_fd(imr_evaluator):
    """Parameter-space gradients through the full IMR (inspiral surrogate,
    transition window, NRSur7dq4 dynamics, blend) against central finite
    differences (measured rel diffs 1.5e-4 / 1.9e-4 / 3.6e-6 for q/e/l).
    The functional weights the valid region only; window indices are stable
    over the FD step at this parameter point."""
    import jax
    import jax.numpy as jnp

    _, fn = imr_evaluator
    q, e0, l0 = _EVAL_PARAMS

    def g(qq, ee, ll):
        h, _vlen = fn(qq, ee, ll)
        return jnp.sum(jnp.abs(h) ** 2) * 1e36

    grads = jax.grad(g, argnums=(0, 1, 2))(q, e0, l0)
    eps = 1e-6
    fd = (
        (g(q + eps, e0, l0) - g(q - eps, e0, l0)) / (2 * eps),
        (g(q, e0 + eps, l0) - g(q, e0 - eps, l0)) / (2 * eps),
        (g(q, e0, l0 + eps) - g(q, e0, l0 - eps)) / (2 * eps),
    )
    for name, a, b in zip("qel", grads, fd):
        a, b = float(a), float(b)
        assert np.isfinite(a), name
        assert abs(a - b) / max(abs(b), 1e-30) < 5e-3, (name, a, b)


def test_evaluator_gpu_matches_cpu(imr_jax):
    import jax

    try:
        gpu = jax.devices("gpu")[0]
    except Exception:
        pytest.skip("no GPU available")
    q, e0, l0 = _EVAL_PARAMS
    _, fn = imr_jax.imr_parameter_space_evaluator(_EVAL_M, _DT, t_start=_EVAL_TS)
    with jax.default_device(jax.devices("cpu")[0]):
        h_cpu, v_cpu = fn(q, e0, l0)
    with jax.default_device(gpu):
        h_gpu, v_gpu = fn(q, e0, l0)
    assert int(v_cpu) == int(v_gpu)
    peak = np.abs(np.asarray(h_cpu)).max()
    assert np.abs(np.asarray(h_gpu) - np.asarray(h_cpu)).max() < 1e-10 * peak


# --- polarizations -----------------------------------------------------------


def test_ylm_2pm2_vs_lal():
    grid = [
        (0.0, 0.0),
        (0.3, 1.1),
        (np.pi / 2, 2.0),
        (2.5, -0.7),
        (np.pi, 0.4),
    ]
    for incl, az in grid:
        for em in (2, -2):
            y_j = complex(generator_jax.spin_weighted_ylm_2pm2_jax(incl, az, em))
            y_l = lal.SpinWeightedSphericalHarmonic(incl, az, -2, 2, em)
            assert abs(y_j - y_l) < 1e-14, (incl, az, em)


def test_polarizations_combination_vs_numpy(imr_jax):
    """Same modes through both combination steps must agree to round-off."""
    from esigmapy.utils import get_polarizations_from_multipoles

    m1, m2, e0, l0, t_start = _IMR_CONFIGS[0]
    t_j, modes_j = imr_jax(
        m1 + m2,
        (m1 / m2, e0, l0),
        _DT,
        t_start=t_start,
        include_conjugate_modes=True,
    )
    modes_np = {lm: np.asarray(h) for lm, h in modes_j.items()}
    incl, phic = 0.6, 0.9
    hp_n, hc_n = get_polarizations_from_multipoles(modes_np, incl, np.pi / 2 - phic)
    hp_j, hc_j = generator_jax.polarizations_from_modes_jax(
        {lm: jnp.asarray(h) for lm, h in modes_j.items()}, incl, phic
    )
    scale = np.abs(hp_n).max()
    assert np.abs(np.asarray(hp_j) - hp_n).max() < 1e-13 * scale
    assert np.abs(np.asarray(hc_j) - hc_n).max() < 1e-13 * scale


def test_polarizations_method_vs_numpy_waveform(imr_jax):
    """Class polarizations vs get_imr_esigmasur_waveform at Tier D bounds."""
    from esigmapy.surrogate.generator import get_imr_esigmasur_waveform

    m1, m2, e0, l0, t_start = _IMR_CONFIGS[0]
    incl, phic = 0.6, 0.0
    hp_n, hc_n = get_imr_esigmasur_waveform(
        mass1=m1,
        mass2=m2,
        delta_t=_DT,
        reference_eccentricity=e0,
        reference_mean_anomaly=l0,
        t_start=t_start,
        inclination=incl,
        coa_phase=phic,
    )
    times, hp_j, hc_j = imr_jax.polarizations(
        m1 + m2,
        (m1 / m2, e0, l0),
        _DT,
        t_start=t_start,
        inclination=incl,
        coa_phase=phic,
    )
    assert float(hp_n.sample_times[0]) == pytest.approx(times[0], abs=1e-9)
    n = min(len(hp_j), len(hp_n))
    scale = np.abs(np.asarray(hp_n.data)).max()
    # inspiral region agrees at backend-equivalence level; overall at the
    # Tier D window-snapping level
    assert np.abs(hp_j[: n // 3] - np.asarray(hp_n.data)[: n // 3]).max() < 1e-6 * scale
    assert np.abs(hp_j[:n] - np.asarray(hp_n.data)[:n]).max() < 0.25 * scale
    assert np.abs(hc_j[:n] - np.asarray(hc_n.data)[:n]).max() < 0.25 * scale


def test_polarizations_inclination_gradient(imr_jax, imr_evaluator):
    """d(hp)/d(inclination) through evaluator + polarization combine."""
    import jax
    import jax.numpy as jnp

    _, fn = imr_evaluator
    q, e0, l0 = _EVAL_PARAMS

    def g(incl):
        h, _ = fn(q, e0, l0)
        modes = {(2, 2): h, (2, -2): jnp.conj(h)}
        hp, _hc = generator_jax.polarizations_from_modes_jax(modes, incl, 0.0)
        return jnp.sum(hp**2) * 1e36

    incl = 0.6
    ad = float(jax.grad(g)(incl))
    eps = 1e-6
    fd = float((g(incl + eps) - g(incl - eps)) / (2 * eps))
    assert np.isfinite(ad)
    assert abs(ad - fd) / max(abs(fd), 1e-30) < 1e-5
