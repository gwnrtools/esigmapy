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
