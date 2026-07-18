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
