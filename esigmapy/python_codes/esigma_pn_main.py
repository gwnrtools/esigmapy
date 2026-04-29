import numpy as np
from scipy.integrate import ode
from scipy.interpolate import CubicSpline
import math
from .esigma_pn_inspiral import *
from .esigma_go_terms_draft import *
import lal

# Constants (LAL equivalents)
LAL_PI = math.pi
LAL_MTSUN_SI = 4.925491025543576e-6  # Solar mass in seconds
RadiationPNOrderDefault = 8  # Default radiation reaction PN order (3PN)

import os
# from dataclasses import dataclass, field

# ------------------------------------------------------------------ #
# Constants / defaults (set these to match your C macros)
# ------------------------------------------------------------------ #
L_MIN = 2
L_MAX = 4
ONLY_LeqM_MODES = False
ModePNOrderDefault = 8
LAL_MRSUN_SI = 1.476625061404649e3   # solar mass in metres

from dataclasses import dataclass

@dataclass
class Params:
    eta: float
    radiation_pn_order: int
    m1: float
    m2: float
    S1z: float
    S2z: float
# # ------------------------------------------------------------------ #
# # Helper: populate powers of orbital variables
# # (mirrors PopulateKeplerParams + kepler_vars struct)
# # ------------------------------------------------------------------ #
# @dataclass
# class KeplerVars:
#     eta:   float = 0.0
#     eta2:  float = 0.0
#     eta3:  float = 0.0
#     eta4:  float = 0.0
#     eta5:  float = 0.0
#     eta6:  float = 0.0

#     Mtot:  float = 0.0
#     Mtot2: float = 0.0
#     Mtot3: float = 0.0
#     Mtot4: float = 0.0
#     Mtot5: float = 0.0
#     Mtot6: float = 0.0

#     S1z:   float = 0.0
#     S1z2:  float = 0.0
#     S1z3:  float = 0.0
#     S1z4:  float = 0.0
#     S1z5:  float = 0.0
#     S1z6:  float = 0.0
#     S1z7:  float = 0.0
#     S1z8:  float = 0.0
#     S1z9:  float = 0.0
#     S1z10: float = 0.0
#     S1z11: float = 0.0
#     S1z12: float = 0.0
#     S1z13: float = 0.0
#     S1z14: float = 0.0
#     S1z15: float = 0.0
#     S1z16: float = 0.0

#     S2z:  float = 0.0
#     S2z2: float = 0.0
#     S2z3: float = 0.0
#     S2z4: float = 0.0
#     S2z5: float = 0.0
#     S2z6: float = 0.0


# def _build_kepler_vars(eta: float, total_mass: float,
#                        S1z: float, S2z: float) -> KeplerVars:
#     """Pre-compute and cache all integer powers used by the mode functions."""
#     kv = KeplerVars()
#     kv.eta  = eta
#     kv.eta2 = eta**2;  kv.eta3 = eta**3;  kv.eta4 = eta**4
#     kv.eta5 = eta**5;  kv.eta6 = eta**6

#     kv.Mtot  = total_mass
#     kv.Mtot2 = total_mass**2;  kv.Mtot3 = total_mass**3
#     kv.Mtot4 = total_mass**4;  kv.Mtot5 = total_mass**5
#     kv.Mtot6 = total_mass**6

#     kv.S1z = S1z
#     for p in range(2, 17):
#         setattr(kv, f"S1z{p}", S1z**p)

#     kv.S2z = S2z
#     for p in range(2, 7):
#         setattr(kv, f"S2z{p}", S2z**p)

#     return kv


# ------------------------------------------------------------------ #
# 1.  compute_mode_from_dynamics
# ------------------------------------------------------------------ #
def compute_mode_from_dynamics(
    l: int,
    m: int,
    x_vec:       np.ndarray,   # PN expansion parameter (length N)
    phi_vec:     np.ndarray,   # orbital phase
    phi_dot_vec: np.ndarray,   # d(phi)/dt
    r_vec:       np.ndarray,   # orbital separation
    r_dot_vec:   np.ndarray,   # d(r)/dt
    mass1:  float,
    mass2:  float,
    S1z:    float,
    S2z:    float,
    R:      float,             # source distance (m)
    vpnorder: int,
) -> np.ndarray:               # complex128, length N
    """
    Compute the (l, m) spin-weighted spherical harmonic mode h_lm
    at every time step from the orbital dynamics arrays.

    Mirrors the static C function compute_mode_from_dynamics().
    """
    total_mass = mass1 + mass2
    eta        = (mass1 * mass2) / total_mass**2

    # kv = _build_kepler_vars(eta, total_mass, S1z, S2z)

    length  = len(x_vec)
    h_lm    = np.zeros(length, dtype=complex)

    for i in range(length):
        # populate_kepler_params(
        #     kv,
        #     e=0.0,
        #     x=x_vec[i],
        #     r=r_vec[i] * total_mass,
        #     r_dot=r_dot_vec[i],
        #     phi_dot=phi_dot_vec[i] / total_mass,
        # )
        h_lm[i] = (
            hlmGOresult(
                l, m, total_mass, eta,
                r_vec[i] * total_mass,
                r_dot_vec[i],
                phi_vec[i],
                phi_dot_vec[i] / total_mass,
                R, vpnorder, S1z, S2z,
                x_vec[i],
            )
            * LAL_MRSUN_SI
        )

    return h_lm


# ------------------------------------------------------------------ #
# 2.  compute_strain_from_dynamics
# ------------------------------------------------------------------ #
def compute_strain_from_dynamics(
    x_vec:       np.ndarray,
    phi_vec:     np.ndarray,
    phi_dot_vec: np.ndarray,
    r_vec:       np.ndarray,
    r_dot_vec:   np.ndarray,
    mass1:  float,
    mass2:  float,
    S1z:    float,
    S2z:    float,
    euler_iota: float,
    euler_beta: float,
    R:          float,
    vpnorder:   int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Sum all (l, m) modes weighted by spin-weighted spherical harmonics
    to produce the + and × strain polarizations.

    Mirrors the static C function compute_strain_from_dynamics().
    """
    length  = len(x_vec)
    h_plus  = np.zeros(length, dtype=float)
    h_cross = np.zeros(length, dtype=float)

    for ell in range(L_MIN, L_MAX + 1):
        for em in range(-ell, ell + 1):

            if ONLY_LeqM_MODES and ell != abs(em):
                continue

            hlm = compute_mode_from_dynamics(
                ell, em,
                x_vec, phi_vec, phi_dot_vec, r_vec, r_dot_vec,
                mass1, mass2, S1z, S2z, R, vpnorder,
            )

            ylm = lal.SpinWeightedSphericalHarmonic(
                euler_iota, euler_beta, -2, ell, em
            )

            hlm_times_ylm = hlm * ylm
            h_plus  += hlm_times_ylm.real
            h_cross -= hlm_times_ylm.imag

    return h_plus, h_cross


# ------------------------------------------------------------------ #
# 3.  XLALSimInspiralesigmaModeFromDynamics  (public)
# ------------------------------------------------------------------ #
def inspiral_esigma_mode_from_dynamics(
    l: int,
    m: int,
    t_vector:       np.ndarray,
    x_vector:       np.ndarray,
    phi_vector:     np.ndarray,
    phi_dot_vector: np.ndarray,
    r_vector:       np.ndarray,
    r_dot_vector:   np.ndarray,
    mass1: float,
    mass2: float,
    S1z:   float,
    S2z:   float,
    R:     float,
) -> np.ndarray:               # complex128
    """
    Public wrapper: compute a single (l, m) waveform mode from dynamics.

    Note: in the C version this modifies t_vec, r_vec, phi_dot_vec in place
    by scaling by total_mass.  Here we keep the arrays immutable and do the
    scaling inside compute_mode_from_dynamics (matches the C end result).
    """
    mode_pn_order = int(os.environ.get("ModePNOrder", ModePNOrderDefault))
    print(mode_pn_order)

    return compute_mode_from_dynamics(
        l, m,
        x_vector, phi_vector, phi_dot_vector,
        r_vector, r_dot_vector,
        mass1, mass2, S1z, S2z, R, mode_pn_order,
    )


# ------------------------------------------------------------------ #
# 4.  XLALSimInspiralesigmaStrainFromDynamics  (public)
# ------------------------------------------------------------------ #
def esigma_strain_from_dynamics(
    t_vector:       np.ndarray,
    x_vector:       np.ndarray,
    phi_vector:     np.ndarray,
    phi_dot_vector: np.ndarray,
    r_vector:       np.ndarray,
    r_dot_vector:   np.ndarray,
    mass1:      float,
    mass2:      float,
    S1z:        float,
    S2z:        float,
    euler_iota: float,
    euler_beta: float,
    R:          float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Public wrapper: compute h+ and hx polarizations from dynamics arrays.
    Returns (h_plus, h_cross) as 1-D float64 arrays.
    """
    mode_pn_order = int(os.environ.get("ModePNOrder", ModePNOrderDefault))

    return compute_strain_from_dynamics(
        x_vector, phi_vector, phi_dot_vector,
        r_vector, r_dot_vector,
        mass1, mass2, S1z, S2z,
        euler_iota, euler_beta, R, mode_pn_order,
    )


# ------------------------------------------------------------------ #
# 5.  x_model_eccbbh_inspiral_waveform  (internal, called by esigma)
# ------------------------------------------------------------------ #
def x_model_eccbbh_inspiral_waveform(
    mass1:          float,   # solar masses
    mass2:          float,
    S1z:            float,
    S2z:            float,
    e_init:         float,
    f_gw_init:      float,   # Hz
    distance:       float,   # metres
    mean_anom_init: float,
    ode_eps:        float,
    euler_iota:     float,
    euler_beta:     float,
    sampling_rate:  float,   # Hz
) -> tuple[np.ndarray, np.ndarray]:
    """
    Drive the full esigma inspiral:
      1. integrate orbital dynamics,
      2. project onto polarizations.

    Returns (h_plus, h_cross) as float64 arrays.
    """
    # --- Step 1: orbital dynamics ----------------------------------- #
    dyn = inspiral_esigma_dynamics(
        mass1, mass2, S1z, S2z,
        e_init, f_gw_init, mean_anom_init,
        ode_eps, sampling_rate,
    )

    # --- Step 2: strain from dynamics ------------------------------- #
    h_plus, h_cross = esigma_strain_from_dynamics(
        dyn["time_evol"],
        dyn["x_evol"],
        dyn["phi_evol"],
        dyn["phi_dot_evol"],
        dyn["r_evol"],
        dyn["r_dot_evol"],
        mass1, mass2, S1z, S2z,
        euler_iota, euler_beta, distance,
    )

    return h_plus, h_cross

def inspiral_esigma_dynamics(
    mass1,          # mass1 in solar mass
    mass2,          # mass2 in solar mass
    S1z,            # z-component of spin of companion 1
    S2z,            # z-component of spin of companion 2
    e_init,         # initial eccentricity
    f_gw_init,      # initial GW frequency
    mean_anom_init, # initial mean anomaly
    ode_eps,        # tolerance (relative)
    sampling_rate,  # sample rate in Hz
):
    """
    Compute ESIGMA orbital dynamics via ODE integration, then interpolate
    to a uniform time grid.

    Returns a dict with keys:
        time_evol, x_evol, eccentricity_evol, mean_ano_evol,
        phi_evol, phi_dot_evol, r_evol, r_dot_evol
    Each value is a 1D numpy array sampled at 1/sampling_rate intervals.
    """

    # ------------------------------------------------------------------ #
    # Input validation
    # ------------------------------------------------------------------ #
    if not (0.0 <= e_init < 1.0):
        raise ValueError("Invalid eccentricity, must be in range [0, 1)")

    # ------------------------------------------------------------------ #
    # Mass / PN bookkeeping
    # ------------------------------------------------------------------ #
    total_mass    = mass1 + mass2
    reduced_mass  = mass1 * mass2 / total_mass
    eta           = reduced_mass / total_mass   # symmetric mass ratio

    omega_init = LAL_PI * f_gw_init * LAL_MTSUN_SI
    x_init     = (total_mass * omega_init) ** (2.0 / 3.0)

    # PN / radiation order (mirror the C env-var logic; default hardcoded)
    import os
    rad_pn_order = int(os.environ.get("RadiationPNOrder", RadiationPNOrderDefault))

    params = Params(
                    eta=eta,
                    radiation_pn_order=rad_pn_order,
                    m1=mass1,
                    m2=mass2,
                    S1z=S1z,
                    S2z=S2z,
                    )

    # ------------------------------------------------------------------ #
    # Termination condition: ISCO
    # ------------------------------------------------------------------ #
    TRANS = float(os.environ.get("InspiralEndRadius", 4.0))
    f_gw_isco = 1.0 / (TRANS * math.sqrt(TRANS) * LAL_PI * total_mass)
    x_final   = (LAL_PI * total_mass * f_gw_isco) ** (2.0 / 3.0)

    # ------------------------------------------------------------------ #
    # Time step in geometric units
    # ------------------------------------------------------------------ #
    dt_sec = 1.0 / sampling_rate                        # seconds
    dt     = dt_sec / (total_mass * LAL_MTSUN_SI)       # geometric (M)

    # ------------------------------------------------------------------ #
    # Initial conditions  y = [x, e, l (mean anomaly), phi]
    # ------------------------------------------------------------------ #
    u0      = pn_kepler_equation(eta, x_init, e_init, mean_anom_init)
    r0      = separation(u0, eta, x_init, e_init, mass1, mass2, S1z, S2z)
    y0_dot  = eccentric_x_model_odes(0.0, [x_init, e_init, mean_anom_init, 0.0], params)
    phi_dot0 = y0_dot[3]

    y0 = np.array([x_init, e_init, mean_anom_init, 0.0])

    # ------------------------------------------------------------------ #
    # ODE integration (RK45 via scipy)
    # ------------------------------------------------------------------ #
    def rhs(t, y):
        return eccentric_x_model_odes(t, y, params)

    solver = ode(rhs).set_integrator(
        "dopri5",           # RK4(5) Dormand-Prince ≈ gsl rkf45
        rtol=ode_eps,
        atol=0.0,
        nsteps=10_000,
        # max_hnil=0,
        first_step=2 * LAL_PI / phi_dot0 / 100,
    )
    solver.set_initial_value(y0, 0.0)

    MAX_SAMPLES = 2048 * 16384

    t_list       = [0.0]
    x_list       = [x_init]
    e_list       = [e_init]
    l_list       = [mean_anom_init]
    phi_list     = [0.0]
    phi_dot_list = [phi_dot0]
    r_list       = [r0]
    u_list       = [u0]

    bad_number = False

    while solver.successful():
        solver.integrate(solver.t + dt)

        y = solver.y

        # NaN / Inf guard
        if not np.all(np.isfinite(y)):
            bad_number = True
            break

        t_list.append(solver.t)
        x_list.append(y[0])
        e_list.append(y[1])
        l_list.append(y[2])
        phi_list.append(y[3])

        yd = eccentric_x_model_odes(solver.t, y, params)
        phi_dot_list.append(yd[3])

        ui = pn_kepler_equation(eta, y[0], y[1], y[2])
        ri = separation(ui, eta, y[0], y[1], mass1, mass2, S1z, S2z)
        u_list.append(ui)
        r_list.append(ri)

        if len(t_list) >= MAX_SAMPLES:
            break

        if y[0] >= x_final:
            break

    # Convert to arrays
    t_arr       = np.array(t_list)
    x_arr       = np.array(x_list)
    e_arr       = np.array(e_list)
    l_arr       = np.array(l_list)
    phi_arr     = np.array(phi_list)
    phi_dot_arr = np.array(phi_dot_list)
    r_arr       = np.array(r_list)

    final_i = len(t_arr) if not bad_number else len(t_arr)  # same either way here
    if final_i < 4:
        raise RuntimeError("Integration produced fewer than 4 points; cannot interpolate.")

    # ------------------------------------------------------------------ #
    # Uniform-grid interpolation
    # ------------------------------------------------------------------ #
    t_final = t_arr[-1]
    Length  = int(math.ceil(t_final / dt))

    if Length < 2:
        raise RuntimeError("Output length < 2; waveform too short.")

    uniform_t = np.arange(Length) * dt   # [0, dt, 2·dt, …]

    def interp_uniform(t_raw, y_raw):
        cs = CubicSpline(t_raw, y_raw)
        return cs(uniform_t)

    def interp_deriv_uniform(t_raw, y_raw):
        cs = CubicSpline(t_raw, y_raw)
        return cs(uniform_t, 1)   # first derivative

    uniform_x       = interp_uniform(t_arr, x_arr)
    uniform_phi     = interp_uniform(t_arr, phi_arr)
    uniform_phi_dot = interp_uniform(t_arr, phi_dot_arr)
    uniform_r       = interp_uniform(t_arr, r_arr)
    uniform_r_dot   = interp_deriv_uniform(t_arr, r_arr)
    uniform_e       = interp_uniform(t_arr, e_arr)
    uniform_l       = interp_uniform(t_arr, l_arr)

    # Convert time back to seconds for the caller
    uniform_t_sec = uniform_t * (total_mass * LAL_MTSUN_SI)

    return {
        "time_evol":          uniform_t_sec,
        "x_evol":             uniform_x,
        "eccentricity_evol":  uniform_e,
        "mean_ano_evol":      uniform_l,
        "phi_evol":           uniform_phi,
        "phi_dot_evol":       uniform_phi_dot,
        "r_evol":             uniform_r,
        "r_dot_evol":         uniform_r_dot,
    }