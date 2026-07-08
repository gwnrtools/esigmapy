"""Tests for the numpy surrogate backend (esigmapy.surrogate.surrogate).

Focus areas:

* the vector-valued TPI parameter-space evaluation (_eval_shared_ref /
  _shared_ref_interp) agrees component-by-component with the per-fit scalar
  TPI interpolants it was stacked from -- this pins the FromComponentSplines
  component ordering;
* the fused uniform-grid fast path agrees with the generic np.interp path
  (exercised via user-supplied `times`) to the fastmath-carve-out level;
* basic waveform / time-grid / phase invariants across the parameter space,
  including the near-circular fallback and the orbital-variables path.

Skipped when ESIGMASUR_DATA_PATH is not configured.
"""

import os

import numpy as np
import pytest

from esigmapy.surrogate.surrogate import CircularSurrogate, EccentricSurrogate

_DATA = os.environ.get("ESIGMASUR_DATA_PATH")
if not _DATA or not os.path.isdir(os.path.join(_DATA, "circ_sur_data")):
    pytest.skip("ESIGMASUR_DATA_PATH not configured", allow_module_level=True)

_CIRC_DIR = os.path.join(_DATA, "circ_sur_data")
_ECC_DIR = os.path.join(_DATA, "ecc_sur_data")

_DT = 0.000244140625
# Fused numba kernels use fastmath and best-fit-uniform index arithmetic, so
# fast-path vs generic-path agreement is held to the AGENTS.md fastmath
# carve-out level rather than machine epsilon.
_FASTMATH_RTOL = 1e-8


def _rel(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return np.abs(a - b) / np.maximum(np.abs(b), 1e-30)


@pytest.fixture(scope="module")
def ecc():
    return EccentricSurrogate(ecc_data_dir=_ECC_DIR, circ_data_dir=_CIRC_DIR)


@pytest.fixture(scope="module")
def circ(ecc):
    return ecc.circ_sur


# --- parameter-space evaluation: vector interpolant vs scalar fits ---------


@pytest.mark.parametrize(
    "point",
    [
        (0.2455, 0.25, 2.1),  # generic interior point
        (0.16528925619834708, 0.431, 6.28),  # near q_max / e_ref_max / l_ref_max
        (0.25, 5.0e-7, 0.0),  # q=1, smallest supported e_ref, l_ref=0
    ],
)
def test_shared_ref_vector_matches_scalar_fits(ecc, point):
    """The stacked vector interpolant must reproduce every per-fit scalar TPI
    interpolant at the same reference point, in the recorded component order."""
    eta, e_ref, l_ref = point
    vals = ecc._eval_shared_ref(eta, e_ref, l_ref)
    x = np.array([eta, e_ref, l_ref])
    scalar_fits = {
        "e": ecc.fit["e"],
        "res_circ_phase": ecc.fit["res_circ_phase"],
        "shifted_mean_anomaly": ecc.fit["shifted_mean_anomaly"],
        "mean_anomaly_offset": ecc.fit["mean_anomaly_offset_fit"],
        "x": ecc.fit["x"],
    }
    assert set(vals) == set(scalar_fits)
    for name, fits in scalar_fits.items():
        expected = np.array([fit.TPInterpolationND(x) for fit in fits])
        assert vals[name].shape == expected.shape
        np.testing.assert_allclose(
            vals[name], expected, rtol=1e-14, atol=1e-15, err_msg=name
        )


def test_shared_ref_slices_cover_all_components(ecc):
    """The per-piece slices must tile the full component vector exactly."""
    n_total = int(np.prod(ecc._shared_ref_interp.GetSplineCoefficientsND().shape[3:]))
    covered = sorted(
        i for sl in ecc._shared_ref_slices.values() for i in range(sl.start, sl.stop)
    )
    assert covered == list(range(n_total))


def test_res_amp_res_phase_fit_counts_match_eim(ecc):
    """One parametric fit per EIM node for the residual pieces."""
    assert len(ecc.fit["res_amp"]) == len(ecc.ei_indices["res_amp"])
    assert len(ecc.fit["res_phase"]) == len(ecc.ei_indices["res_phase"])


# --- fast path vs generic path ----------------------------------------------


@pytest.mark.parametrize("t_start", [-136.0, -1.0])
def test_ecc_fast_path_matches_generic_path(ecc, t_start):
    """The fused uniform-grid kernel must agree with the generic np.interp
    path when handed the identical (uniform) time grid via `times`."""
    kw = dict(M=10.0, params=(2.3, 0.25, 2.1), delta_t=_DT, t_start=t_start)
    t_fast, h_fast = ecc(**kw)
    t_gen, h_gen = ecc(M=10.0, params=(2.3, 0.25, 2.1), times=t_fast)
    assert np.array_equal(t_fast, t_gen)
    assert _rel(h_fast, h_gen).max() < _FASTMATH_RTOL


@pytest.mark.parametrize("t_start", [-136.0, -1.0])
def test_circ_fast_path_matches_generic_path(circ, t_start):
    kw = dict(M=10.0, q=2.3, delta_t=_DT, t_start=t_start, remove_initial_phase=True)
    t_fast, h_fast = circ(**kw)
    t_gen, h_gen = circ(M=10.0, q=2.3, times=t_fast, remove_initial_phase=True)
    assert np.array_equal(t_fast, t_gen)
    assert _rel(h_fast, h_gen).max() < _FASTMATH_RTOL


# --- waveform / grid invariants ----------------------------------------------


@pytest.mark.parametrize(
    "M,params,t_start",
    [
        (10.0, (2.3, 0.25, 2.1), -136.0),
        (10.0, (2.3, 0.25, 2.1), -1.0),
        (10.0, (1.1, 0.43, 6.2), -20.0),
        (60.0, (4.0, 0.1, 3.14), -300.0),
        (10.0, (5.9, 0.05, 0.0), -5.0),
    ],
)
def test_ecc_waveform_invariants(ecc, M, params, t_start):
    t, h = ecc(M=M, params=params, delta_t=_DT, t_start=t_start)
    # Deterministic uniform grid starting exactly at t_start.
    assert t[0] == t_start
    assert t.shape == h.shape
    dt = np.diff(t)
    assert np.allclose(dt, _DT, rtol=0, atol=1e-12)
    assert np.all(np.isfinite(h.real)) and np.all(np.isfinite(h.imag))
    assert np.all(np.abs(h) > 0)
    # remove_start_phase=True (default): first sample has zero phase.
    assert abs(h[0].imag) < 1e-12 * abs(h[0].real) + 1e-30
    assert h[0].real > 0


def test_ecc_low_e_falls_back_to_circular(ecc, circ):
    """Below e_ref_min the eccentric surrogate must delegate to the circular
    surrogate evaluated on the same output grid. (The two surrogates have
    different native t_end, so the circular reference is evaluated via `times`
    on the eccentric call's grid.)"""
    t_e, h_e = ecc(M=10.0, params=(2.3, 0.0, 1.0), delta_t=_DT, t_start=-2.0)
    t_c, h_c = circ(M=10.0, q=2.3, times=t_e, remove_initial_phase=True)
    assert np.array_equal(t_e, t_c)
    np.testing.assert_allclose(h_e, h_c, rtol=1e-13, atol=0)


def test_ecc_small_e_close_to_circular(ecc, circ):
    """Just above e_ref_min the eccentric waveform should stay close to the
    circular one (the residual pieces scale with eccentricity)."""
    t_e, h_e = ecc(M=10.0, params=(2.3, 1.0e-4, 1.0), delta_t=_DT, t_start=-2.0)
    _, h_c = circ(M=10.0, q=2.3, times=t_e, remove_initial_phase=True)
    assert _rel(np.abs(h_e), np.abs(h_c)).max() < 1e-2


def test_ecc_param_range_validation(ecc):
    with pytest.raises(ValueError):
        ecc(M=10.0, params=(0.5, 0.1, 0.0), delta_t=_DT, t_start=-1.0)
    with pytest.raises(ValueError):
        ecc(M=10.0, params=(2.0, 0.6, 0.0), delta_t=_DT, t_start=-1.0)
    with pytest.raises(ValueError):
        ecc(M=10.0, params=(2.0, 0.1, 7.0), delta_t=_DT, t_start=-1.0)


def test_ecc_time_range_validation(ecc):
    with pytest.raises(ValueError):
        ecc(M=10.0, params=(2.3, 0.1, 0.0), delta_t=_DT, t_start=-1.0e4)


def test_ecc_requires_delta_t_or_times(ecc):
    with pytest.raises(ValueError):
        ecc(M=10.0, params=(2.3, 0.1, 0.0))


# --- orbital variables --------------------------------------------------------


def test_ecc_orbital_variables(ecc):
    t, orb, h = ecc(
        M=10.0,
        params=(2.3, 0.25, 2.1),
        delta_t=_DT,
        t_start=-10.0,
        return_orbital_variables=True,
    )
    assert set(orb) == {"e", "l", "x"}
    for key in ("e", "l", "x"):
        assert orb[key].shape == t.shape
        assert np.all(np.isfinite(orb[key]))
    # eccentricity decays during inspiral and stays physical
    assert np.all(orb["e"] >= 0.0) and np.all(orb["e"] < 0.5)
    assert orb["e"][0] > orb["e"][-1]
    # mean anomaly starts in [0, 2pi) and is increasing
    assert 0.0 <= orb["l"][0] < 2 * np.pi
    assert np.all(np.diff(orb["l"]) > 0)
    # x = (M*omega)^(2/3) proxy: positive and growing toward merger
    assert np.all(orb["x"] > 0)
    assert orb["x"][-1] > orb["x"][0]


# --- circular surrogate -------------------------------------------------------


@pytest.mark.parametrize("t_start", [-136.0, -1.0])
def test_circ_amp_phase_consistency(circ, t_start):
    """amp/phase-only output must reassemble into the returned mode."""
    kw = dict(M=10.0, q=2.3, delta_t=_DT, t_start=t_start)
    amp, phase = circ(return_amp_phase_only=True, **kw)
    t, h = circ(**kw)
    np.testing.assert_allclose(amp * np.exp(-1j * phase), h, rtol=1e-12, atol=1e-30)
