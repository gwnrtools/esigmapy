"""Compare the JAX surrogate backend against the numpy backend.

Skipped unless jax and the TPI JAX backend (TPI_jax) are importable and the
surrogate data is available. The JAX backend reproduces the numpy backend to
the fastmath-carve-out level: it uses the same best-fit-uniform index-arithmetic
interpolation, and XLA reassociates the gemvs/transcendentals, so differences
sit at ~1e-10 relative on the waveform rather than at machine epsilon.
"""

import os

import numpy as np
import pytest

pytest.importorskip("jax")

try:
    from esigmapy.surrogate.surrogate_jax import (
        CircularSurrogateJAX,
        EccentricSurrogateJAX,
    )
except Exception as exc:  # pragma: no cover - environment dependent
    pytest.skip(f"JAX surrogate backend unavailable: {exc}", allow_module_level=True)

from esigmapy.surrogate.surrogate import CircularSurrogate, EccentricSurrogate

_DATA = os.environ.get("ESIGMASUR_DATA_PATH")
if not _DATA or not os.path.isdir(os.path.join(_DATA, "circ_sur_data")):
    pytest.skip("ESIGMASUR_DATA_PATH not configured", allow_module_level=True)

_CIRC_DIR = os.path.join(_DATA, "circ_sur_data")
_ECC_DIR = os.path.join(_DATA, "ecc_sur_data")

# Waveform (physical output) tolerance and the large-magnitude phase tolerance.
_WF_RTOL = 1e-9
_PHASE_ATOL = 1e-8


def _rel(a, b):
    a = np.asarray(a)
    b = np.asarray(b)
    return np.abs(a - b) / np.maximum(np.abs(b), 1e-30)


@pytest.fixture(scope="module")
def circ_np():
    return CircularSurrogate(_CIRC_DIR)


@pytest.fixture(scope="module")
def circ_jax():
    return CircularSurrogateJAX(_CIRC_DIR)


@pytest.mark.parametrize("t_start", [-136.0, -1.0])
@pytest.mark.parametrize("q", [1.5, 2.3, 5.0])
def test_circular_full_waveform(circ_np, circ_jax, t_start, q):
    kw = dict(M=10.0, q=q, delta_t=0.000244140625, t_start=t_start, t_end=None)
    tn, hn = circ_np(**kw)
    tj, hj = circ_jax(**kw)
    assert np.array_equal(tn, tj), "time grids must match exactly"
    assert _rel(hj, hn).max() < _WF_RTOL


@pytest.mark.parametrize("t_start", [-136.0, -1.0])
def test_circular_amp_phase(circ_np, circ_jax, t_start):
    kw = dict(
        M=20.0,
        q=2.3,
        delta_t=0.000244140625,
        t_start=t_start,
        t_end=None,
        return_amp_phase_only=True,
        remove_initial_phase=True,
    )
    an, pn = circ_np(**kw)
    aj, pj = circ_jax(**kw)
    assert _rel(aj, an).max() < _WF_RTOL
    # initial phase removed -> first sample is exactly zero in both
    assert abs(pj[0]) < 1e-12
    assert np.abs(np.asarray(pj) - np.asarray(pn)).max() < _PHASE_ATOL


def test_circular_user_times(circ_np, circ_jax):
    times = np.linspace(-50.0, -0.5, 8192)
    tn, hn = circ_np(M=15.0, q=3.0, times=times)
    tj, hj = circ_jax(M=15.0, q=3.0, times=times)
    assert np.array_equal(tn, tj)
    assert _rel(hj, hn).max() < _WF_RTOL


@pytest.fixture(scope="module")
def ecc_np():
    return EccentricSurrogate(ecc_data_dir=_ECC_DIR, circ_data_dir=_CIRC_DIR)


@pytest.fixture(scope="module")
def ecc_jax():
    return EccentricSurrogateJAX(ecc_data_dir=_ECC_DIR, circ_data_dir=_CIRC_DIR)


@pytest.mark.parametrize("t_start", [-136.0, -1.0])
@pytest.mark.parametrize(
    "q,e_ref,l_ref",
    [(2.3, 0.43, 1.3), (1.5, 0.1, 0.0), (5.0, 0.25, 2.0 * np.pi)],
)
def test_eccentric_full_waveform(ecc_np, ecc_jax, t_start, q, e_ref, l_ref):
    kw = dict(
        M=10.0,
        params=(q, e_ref, l_ref),
        delta_t=0.000244140625,
        t_start=t_start,
        t_end=None,
    )
    tn, hn = ecc_np(**kw)
    tj, hj = ecc_jax(**kw)
    assert np.array_equal(tn, tj)
    assert _rel(hj, hn).max() < _WF_RTOL


def test_eccentric_e_ref_fallback(ecc_np, ecc_jax):
    kw = dict(
        M=10.0,
        params=(2.3, 1e-8, 1.3),
        delta_t=0.000244140625,
        t_start=-1.0,
        t_end=None,
    )
    tn, hn = ecc_np(**kw)
    tj, hj = ecc_jax(**kw)
    assert np.array_equal(tn, tj)
    assert _rel(hj, hn).max() < _WF_RTOL


def test_eccentric_user_times(ecc_np, ecc_jax):
    times = np.linspace(-30.0, -0.5, 6000)
    tn, hn = ecc_np(M=12.0, params=(3.0, 0.2, 1.0), times=times)
    tj, hj = ecc_jax(M=12.0, params=(3.0, 0.2, 1.0), times=times)
    assert np.array_equal(tn, tj)
    assert _rel(hj, hn).max() < _WF_RTOL
