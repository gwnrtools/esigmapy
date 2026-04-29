import numpy as np
from dataclasses import dataclass

# Storing constants 
M_PI = np.pi
M_PI2 = M_PI ** 2
LOG2 = np.log(2)

#----------------------------------------------------
def rect(r: float, phi: float) -> complex:
    """Convert polar coordinates to rectangular form."""
    return r * np.exp(1j * phi)
#----------------------------------------------------

@dataclass
class CommonVars:
    xp5: float
    logx: float
    b0: float
    r0: float
    logb0: float
    logr0: float
    delta: float
#----------------------------------------------------

H_GO_LM_FUNCS = {}
H_QC_LM_FUNCS = {}
def register_hgo_lm(l, m):
    """
    Decorator to register a GO waveform mode function for a given (l, m).

    This decorator adds the decorated function to the global HLM_FUNCS
    registry, using the tuple (l, m) as the key. It enables clean and
    scalable dispatch of mode-specific functions without requiring a
    large manual lookup table or switch-case logic.

    Parameters
    ----------
    l : int
        Spherical harmonic index (typically 2 ≤ l ≤ 8).
    m : int
        Azimuthal index (-l ≤ m ≤ l, excluding m = 0 if unused).

    Returns
    -------
    decorator : callable
        A decorator that takes a function and registers it in HLM_FUNCS.

    Notes
    -----
    - The decorated function should have a signature compatible with:
          func(params, pno, ...)
      where `pno` is the PN order.
    - If a function is already registered for the same (l, m), it may be
      overwritten unless explicitly guarded against.

    Examples
    --------
    >>> @register_hlm(2, 2)
    ... def hl_2_m_2(params, pno):
    ...     return 0j

    >>> H_GO_LM_FUNCS[(2, 2)] is hGO_2_m_2
    True
    """
    def decorator(func):
        H_GO_LM_FUNCS[(l, m)] = func
        return func
    return decorator

def register_hqc_lm(l, m):
    """Similarly for QC functions"""
    def decorator(func):
        H_QC_LM_FUNCS[(l, m)] = func
        return func
    return decorator

# H22

@register_hgo_lm(2, 2)
def hGO_2_m_2(total_mass: float, eta: float, r: float, rDOT: float,
              PhiDOT: float, vpnorder: int, S1z: float, S2z: float,
              x: float, params: CommonVars) -> complex:
    
    # For black holes kappa and lambda is 1
    kappa1 = 1.0 
    kappa2 = 1.0
    lambda1 = 1.0
    lambda2 = 1.0
    
    # delta = np.sqrt(1 - 4 * eta)
    delta = params.delta
    
    combination_a = (PhiDOT * r + complex(0, 1) * rDOT)
    combination_a3 = combination_a * combination_a * combination_a
    combination_a4 = combination_a3 * combination_a
    combination_a5 = combination_a4 * combination_a

    combination_b = (PhiDOT * r - complex(0, 1) * rDOT)
    combination_b2 = combination_b * combination_b
    combination_b3 = combination_b2 * combination_b

    rDOT2 = rDOT * rDOT
    rDOT3 = rDOT2 * rDOT
    rDOT4 = rDOT3 * rDOT
    rDOT5 = rDOT4 * rDOT
    rDOT6 = rDOT5 * rDOT

    total_mass2 = total_mass * total_mass
    total_mass3 = total_mass2 * total_mass
    total_mass4 = total_mass3 * total_mass

    PhiDOT2 = PhiDOT * PhiDOT
    PhiDOT3 = PhiDOT2 * PhiDOT
    PhiDOT4 = PhiDOT3 * PhiDOT
    PhiDOT5 = PhiDOT4 * PhiDOT
    PhiDOT6 = PhiDOT5 * PhiDOT

    r2 = r * r
    r3 = r2 * r
    r4 = r3 * r
    r5 = r4 * r
    r6 = r5 * r
    
    eta2 = eta * eta
    eta3 = eta2 * eta

    S1z2 = S1z * S1z
    S1z3 = S1z2 * S1z
 
    S2z2 = S2z * S2z
    S2z3 = S2z2 * S2z

    if vpnorder == 0:
        return (total_mass / r + PhiDOT2 * r2 +
                complex(0, 2) * PhiDOT * r * rDOT - rDOT2)

    elif vpnorder == 2:
        return ((21 * total_mass2 * (-10 + eta) -
                 27 * (-1 + 3 * eta) * r2 * combination_b * combination_a3 +
                 total_mass * r *
                     ((11 + 156 * eta) * PhiDOT2 * r2 +
                      complex(0, 10) * (5 + 27 * eta) * PhiDOT * r * rDOT -
                      3 * (15 + 32 * eta) * rDOT2)) /
                (42. * r2))

    elif vpnorder == 3:
        return ((total_mass2 * (complex(0, -1) * rDOT *
                                 ((3 + 3 * delta - 8 * eta) * S1z +
                                  (3 - 3 * delta - 8 * eta) * S2z) +
                             PhiDOT * r *
                                 ((-3 - 3 * delta + 5 * eta) * S1z +
                                  (-3 + 3 * delta + 5 * eta) * S2z))) /
                (3. * r2)) 
        # (<--This is the general orbit term)
        # (This is the quasi-circular limit of the general orbit term-->)
        # - ((-4 * ((1 + delta - eta) * S1z + S2z - (delta + eta) * S2z) * params.x2p5) / 3.) +
        # ((-4 * (S1z + delta * S1z + S2z - delta * S2z - eta * (S1z + S2z)) * params.x2p5) / 3.) 
        # (<--This is Quentins quasi-circular term)

    elif vpnorder == 4:
        return ((6 * total_mass3 * (3028 + 1267 * eta + 158 * eta2) +
                 9 * (83 - 589 * eta + 1111 * eta2) * r3 *
                     combination_b2 * combination_a4 +
                 total_mass2 * r *
                     ((-11891 - 36575 * eta + 13133 * eta2) * PhiDOT2 *
                          r2 +
                      complex(0, 8) * (-773 - 3767 * eta + 2852 * eta2) *
                          PhiDOT * r * rDOT -
                      6 * (-619 + 2789 * eta + 934 * eta2) * rDOT2) -
                 3 * total_mass * r2 *
                     (2 * (-835 - 19 * eta + 2995 * eta2) * PhiDOT4 *
                          r4 +
                      complex(0, 6) * (-433 - 721 * eta + 1703 * eta2) *
                          PhiDOT3 * r3 * rDOT +
                      6 * (-33 + 1014 * eta + 232 * eta2) * PhiDOT2 *
                          r2 * rDOT2 +
                      complex(0, 4) * (-863 + 1462 * eta + 2954 * eta2) *
                          PhiDOT * r * rDOT3 -
                      3 * (-557 + 664 * eta + 1712 * eta2) * rDOT4)) /
                    (1512. * r3) +
                (3 * total_mass3 *
                 (S1z * (4 * eta * S2z + (1 + delta - 2 * eta) * S1z * kappa1) -
                  (-1 + delta + 2 * eta) * S2z2 * kappa2)) /
                    (4. * r3))
                # This is where circular limit terms added and subtracted
                # - ((kappa1 * (1 + delta - 2 * eta) * S1z2 + S2z * (4 * eta * S1z - kappa2 * (-1 + delta + 2 * eta) * S2z)) * params.x3) +
                # ((kappa1 * (1 + delta - 2 * eta) * S1z2 + S2z * (4 * eta * S1z - kappa2 * (-1 + delta + 2 * eta) * S2z)) * params.x3)

    elif vpnorder == 5:
        return ((total_mass2 * eta *
                 (2 * total_mass * (complex(0, -702) * PhiDOT * r + rDOT) +
                  3 * r *
                      (complex(0, -316) * PhiDOT3 * r3 -
                       847 * PhiDOT2 * r2 * rDOT +
                       complex(0, 184) * PhiDOT * r * rDOT2 -
                       122 * rDOT3))) /
                (105. * r3) 
                # Henry et al. QC spin terms
                # + ((2*(56*delta*eta*(-S1z + S2z) + 101*eta*(S1z + S2z) + 132*eta2*(S1z + S2z) - 80*(S1z + delta*S1z + S2z - delta*S2z))*params.x3p5)/63.)
                +
                # Henry et al. ecc spin terms
                ((total_mass2 *
                 (total_mass *
                      ((238 + delta * (238 - 141 * eta) + eta * (-181 + 474 * eta)) *
                           PhiDOT * r * S1z +
                       complex(0, 8) *
                           (55 + delta * (55 - 19 * eta) +
                            2 * eta * (-50 + 43 * eta)) *
                           rDOT * S1z +
                       (238 + delta * (-238 + 141 * eta) + eta * (-181 + 474 * eta)) *
                           PhiDOT * r * S2z +
                       complex(0, 8) *
                           (55 + delta * (-55 + 19 * eta) +
                            2 * eta * (-50 + 43 * eta)) *
                           rDOT * S2z) -
                  r * (PhiDOT * r * rDOT2 *
                           (-((18 * (1 + delta) + 5 * (-63 + 55 * delta) * eta +
                               188 * eta2) *
                              S1z) +
                            (18 * (-1 + delta) + 5 * (63 + 55 * delta) * eta -
                             188 * eta2) *
                                S2z) -
                       complex(0, 2) * rDOT3 *
                           ((-27 * (1 + delta) + 6 * (5 + 7 * delta) * eta -
                             4 * eta2) *
                                S1z +
                            (-27 + 27 * delta + 30 * eta - 42 * delta * eta -
                             4 * eta2) *
                                S2z) +
                       complex(0, 2) * PhiDOT2 * r2 * rDOT *
                           ((51 + 88 * eta * (-3 + 5 * eta) +
                             delta * (51 + 62 * eta)) *
                                S1z +
                            (51 + 88 * eta * (-3 + 5 * eta) -
                             delta * (51 + 62 * eta)) *
                                S2z) +
                       PhiDOT3 * r3 *
                           ((120 * (1 + delta) + (-483 + 83 * delta) * eta +
                             234 * eta2) *
                                S1z +
                            (120 + 3 * eta * (-161 + 78 * eta) -
                             delta * (120 + 83 * eta)) *
                                S2z)))) /
                (84. * r3)))

    elif vpnorder == 6:
        return ((4 * total_mass4 *
                 (-8203424 + 2180250 * eta2 + 592600 * eta3 +
                  15 * eta * (-5503804 + 142065 * M_PI2)) -
             2700 * (-507 + 6101 * eta - 25050 * eta2 + 34525 * eta3) *
                 r4 * combination_b3 * combination_a5 +
             total_mass3 * r *
                 (PhiDOT2 *
                      (337510808 - 198882000 * eta2 +
                       56294600 * eta3 +
                       eta * (183074880 - 6392925 * M_PI2)) *
                      r2 +
                  complex(0, 110) * PhiDOT *
                      (-5498800 - 785120 * eta2 + 909200 * eta3 +
                       3 * eta * (-1849216 + 38745 * M_PI2)) *
                      r * rDOT +
                  2 *
                      (51172744 - 94929000 * eta2 - 5092400 * eta3 +
                       45 * eta * (2794864 + 142065 * M_PI2)) *
                      rDOT2) -
             20 * total_mass2 * r2 *
                 ((-986439 + 1873255 * eta - 9961400 * eta2 +
                   6704345 * eta3) *
                      PhiDOT4 * r4 +
                  complex(0, 4) *
                      (-273687 - 978610 * eta - 4599055 * eta2 +
                       2783005 * eta3) *
                      PhiDOT3 * r3 * rDOT +
                  (-181719 + 19395325 * eta + 8237980 * eta2 +
                   2612735 * eta3) *
                      PhiDOT2 * r2 * rDOT2 +
                  complex(0, 8) *
                      (-234312 + 1541140 * eta + 1230325 * eta2 +
                       1828625 * eta3) *
                      PhiDOT * r * rDOT3 -
                  3 *
                      (-370268 + 1085140 * eta + 2004715 * eta2 +
                       1810425 * eta3) *
                      rDOT4) +
             300 * total_mass * r3 *
                 (4 *
                      (12203 - 36427 * eta - 27334 * eta2 +
                       149187 * eta3) *
                      PhiDOT6 * r6 +
                  complex(0, 2) *
                      (44093 - 68279 * eta - 295346 * eta2 +
                       541693 * eta3) *
                      PhiDOT5 * r5 * rDOT +
                  2 *
                      (27432 - 202474 * eta + 247505 * eta2 +
                       394771 * eta3) *
                      PhiDOT4 * r4 * rDOT2 +
                  complex(0, 2) *
                      (97069 - 383990 * eta - 8741 * eta2 +
                       1264800 * eta3) *
                      PhiDOT3 * r3 * rDOT3 +
                  (-42811 + 53992 * eta + 309136 * eta2 -
                   470840 * eta3) *
                      PhiDOT2 * r2 * rDOT4 +
                  complex(0, 2) *
                      (51699 - 252256 * eta + 131150 * eta2 +
                       681160 * eta3) *
                      PhiDOT * r * rDOT5 -
                  3 *
                      (16743 - 75104 * eta + 26920 * eta2 +
                       207200 * eta3) *
                      rDOT6)) /
                (3.3264e6 * r4)
            # Henry et al. QC spin terms
            # + (((4*(1 + delta)*(-7 + 9*kappa1) - 7*(9 + 17*delta)*eta - 9*(15 + 7*delta)*kappa1*eta + 12*(7 - 17*kappa1)*eta2)*S1z2 + 2*S1z*(complex(0,-42)*(1 + delta - 2*eta) - 84*(1 + delta - eta)*M_PI + eta*(-271 + 288*eta)*S2z) + S2z*(12*(7 - 17*kappa2)*eta2*S2z + 4*(-1 + delta)*(complex(0,21) + 42*M_PI + 7*S2z - 9*kappa2*S2z) + eta*(168*(complex(0,1) + M_PI) + 7*delta*(17 + 9*kappa2)*S2z - 9*(7 + 15*kappa2)*S2z)))*params.x4)/63.
            +
            # Henry et al. ecc spin terms
            (-0.005952380952380952 *
                  (total_mass3 *
                   (2 * total_mass *
                        (S1z * (complex(0, 14) * (1 + delta - 2 * eta) +
                                42 * (1 + delta - 2 * eta) * eta * S1z +
                                kappa1 *
                                    (438 + delta * (438 + 103 * eta) +
                                     eta * (-773 + 108 * eta)) *
                                    S1z) +
                         2 *
                             (complex(0, -7) * (-1 + delta + 2 * eta) +
                              (995 - 192 * eta) * eta * S1z) *
                             S2z -
                         (42 * eta * (-1 + delta + 2 * eta) +
                          kappa2 * (-438 + (773 - 108 * eta) * eta +
                                    delta * (438 + 103 * eta))) *
                             S2z2) +
                    r * (rDOT2 *
                             (complex(0, 56) * (1 + delta - 2 * eta) * S1z +
                              (-56 * (1 + delta - 2 * eta) * eta +
                               kappa1 * (291 * (1 + delta) -
                                         2 * (445 + 154 * delta) * eta +
                                         24 * eta2)) *
                                  S1z2 +
                              4 * eta * (-3 + 44 * eta) * S1z * S2z +
                              S2z * (complex(0, -56) * (-1 + delta + 2 * eta) +
                                     56 * eta * (-1 + delta + 2 * eta) * S2z +
                                     kappa2 *
                                         (291 - 890 * eta + 24 * eta2 +
                                          delta * (-291 + 308 * eta)) *
                                         S2z)) +
                         PhiDOT2 * r2 *
                             (complex(0, 196) * (1 + delta - 2 * eta) * S1z +
                              (56 * eta * (7 + 7 * delta + eta) +
                               kappa1 * (-153 * (1 + delta) -
                                         2 * (62 + 215 * delta) * eta +
                                         804 * eta2)) *
                                  S1z2 +
                              8 * (60 - 187 * eta) * eta * S1z * S2z +
                              S2z * (complex(0, -196) * (-1 + delta + 2 * eta) +
                                     56 * eta * (7 - 7 * delta + eta) * S2z +
                                     kappa2 *
                                         (-153 + 4 * eta * (-31 + 201 * eta) +
                                          delta * (153 + 430 * eta)) *
                                         S2z)) +
                         2 * PhiDOT * r * rDOT *
                             (complex(0, -1) *
                                  (117 * (1 + delta) * kappa1 -
                                   434 * (1 + delta) * eta +
                                   (-23 + 211 * delta) * kappa1 * eta +
                                   2 * (14 - 195 * kappa1) * eta2) *
                                  S1z2 +
                              4 * S1z *
                                  (35 * (1 + delta - 2 * eta) +
                                   complex(0, 1) * (9 - 209 * eta) * eta * S2z) +
                              S2z * (-140 * (-1 + delta + 2 * eta) +
                                     complex(0, 1) *
                                         (-14 * eta * (-31 + 31 * delta + 2 * eta) +
                                          kappa2 * (117 * (-1 + delta) +
                                                    (23 + 211 * delta) * eta +
                                                    390 * eta2)) *
                                         S2z))))) /
                  r4))
            # + Henry et al. QC spinning hereditary terms
            # (((-8 * M_PI * ((1 + delta - eta) * S1z + S2z - (delta + eta) * S2z) * params.x4) / 3.))

    elif vpnorder == 7:
        return (
            # Henry et al QC spin terms
            # ((3318*eta3*(S1z + S2z) + eta*(-504*((7 + delta)*kappa1 - 3*(3 + delta)*lambda1)*S1z3 - 1008*S1z2*(3*kappa1*M_PI - 3*(1 + delta)*S2z + 2*(1 + delta)*kappa1*S2z) + S1z*(17387 + 20761*delta + 1008*S2z*(6*M_PI + (-1 + delta)*(-3 + 2*kappa2)*S2z)) + S2z*(17387 - 20761*delta + 504*S2z*(-6*kappa2*M_PI + (-7 + delta)*kappa2*S2z - 3*(-3 + delta)*lambda2*S2z))) + 2*(2809*(1 + delta)*S1z + 756*(1 + delta)*kappa1*M_PI*S1z2 + 756*(1 + delta)*(kappa1 - lambda1)*S1z3 - (-1 + delta)*S2z*(2809 + 756*S2z*(-(lambda2*S2z) + kappa2*(M_PI + S2z)))) - 2*eta2*(708*delta*(-S1z + S2z) + (S1z + S2z)*(4427 + 1008*(kappa1*S1z2 + S2z*(-2*S1z + kappa2*S2z)))))*params.x4p5)/756.
            # +
            # Henry et al. ecc+spin terms
            ((total_mass2 *
                 (-3 * total_mass * r *
                      (complex(0, -16) * rDOT3 *
                           (complex(0, 12) * eta * (-16703 + 4427 * eta) +
                            35 * eta *
                                (4578 + eta * (-4288 + 765 * eta) +
                                 delta * (3748 + 802 * eta)) *
                                S1z +
                            35 *
                                (delta * (942 - 2 * eta * (1874 + 401 * eta)) +
                                 eta * (4578 + eta * (-4288 + 765 * eta))) *
                                S2z -
                            32970 * (S1z + delta * S1z + S2z)) +
                       4 * PhiDOT * r * rDOT2 *
                           (-338520 * eta3 * (S1z + S2z) +
                            48930 * (S1z + delta * S1z + S2z - delta * S2z) +
                            eta * (complex(0, 3420696) -
                                  35 * (14154 + 21167 * delta) * S1z -
                                  495390 * S2z + 740845 * delta * S2z) +
                            eta2 * (complex(0, -612336) +
                                           245 * (3566 - 1565 * delta) * S1z +
                                           245 * (3566 + 1565 * delta) * S2z)) +
                       PhiDOT3 * r3 *
                           (2515380 * eta3 * (S1z + S2z) -
                            5 * eta *
                                (complex(0, 1859936) +
                                 7 * (82329 + 37061 * delta) * S1z +
                                 7 * (82329 - 37061 * delta) * S2z) -
                            128100 * (S1z + delta * S1z + S2z - delta * S2z) +
                            4 * eta2 *
                                (complex(0, 381348) +
                                 35 * (-18505 + 1777 * delta) * S1z -
                                 35 * (18505 + 1777 * delta) * S2z)) +
                       complex(0, 8) * PhiDOT2 * r2 * rDOT *
                           (779100 * eta3 * (S1z + S2z) +
                            5 * eta *
                                (complex(0, 828806) +
                                 7 * (4839 + 5971 * delta) * S1z +
                                 7 * (4839 - 5971 * delta) * S2z) -
                            62475 * (S1z + delta * S1z + S2z - delta * S2z) +
                            eta2 * (complex(0, -976002) +
                                           35 * (-29599 + 3109 * delta) * S1z -
                                           35 * (29599 + 3109 * delta) * S2z))) +
                  3 * r2 *
                      (complex(0, 4) * PhiDOT2 * r2 * rDOT3 *
                           (complex(0, 2) * (65451 - 350563 * eta) * eta +
                            105 * eta *
                                (6020 + delta * (3114 + 411 * eta) +
                                 eta * (-4513 + 9136 * eta)) *
                                S1z +
                            105 *
                                (-3 * delta * (-408 + eta * (1038 + 137 * eta)) +
                                 eta * (6020 + eta * (-4513 + 9136 * eta))) *
                                S2z -
                            128520 * (S1z + delta * S1z + S2z)) +
                       PhiDOT * r * rDOT4 *
                           (complex(0, -128) * eta * (2487 + 18334 * eta) -
                            105 * eta *
                                (8689 + 8 * eta * (-687 + 1402 * eta) +
                                 delta * (2143 + 8212 * eta)) *
                                S1z +
                            105 *
                                (eta * (-8689 + 8 * (687 - 1402 * eta) * eta) +
                                 delta * (-1470 + eta * (2143 + 8212 * eta))) *
                                S2z +
                            154350 * (S1z + delta * S1z + S2z)) -
                       complex(0, 6) * rDOT5 *
                           (-11760 * eta3 * (S1z + S2z) +
                            4 * eta2 *
                                (complex(0, 57338) +
                                 35 * (569 + 301 * delta) * S1z +
                                 35 * (569 - 301 * delta) * S2z) -
                            16 * eta *
                                (complex(0, -2517) + 35 * (72 + 71 * delta) * S1z +
                                 35 * (72 - 71 * delta) * S2z) +
                            8715 * (S1z + delta * S1z + S2z - delta * S2z)) +
                       PhiDOT5 * r5 *
                           (2263380 * eta3 * (S1z + S2z) +
                            3 * eta *
                                (complex(0, 653432) +
                                 35 * (13283 + 7839 * delta) * S1z +
                                 35 * (13283 - 7839 * delta) * S2z) -
                            219240 * (S1z + delta * S1z + S2z - delta * S2z) +
                            4 * eta2 *
                                (complex(0, -291268) +
                                 105 * (-7669 + 165 * delta) * S1z -
                                 105 * (7669 + 165 * delta) * S2z)) +
                       complex(0, 14) * PhiDOT4 * r4 * rDOT *
                           (385170 * eta3 * (S1z + S2z) +
                            15 * eta *
                                (complex(0, 16914) + (8633 + 2267 * delta) * S1z +
                                 8633 * S2z - 2267 * delta * S2z) -
                            18630 * (S1z + delta * S1z + S2z - delta * S2z) +
                            2 * eta2 *
                                (complex(0, -67904) +
                                 15 * (-13932 + 679 * delta) * S1z -
                                 15 * (13932 + 679 * delta) * S2z)) +
                       6 * PhiDOT3 * r3 * rDOT2 *
                           (-6720 * eta3 * (S1z + S2z) +
                            5 * eta *
                                (complex(0, -241704) +
                                 7 * (4377 + 3083 * delta) * S1z +
                                 7 * (4377 - 3083 * delta) * S2z) -
                            36960 * (S1z + delta * S1z + S2z - delta * S2z) +
                            eta2 * (complex(0, -364384) +
                                           35 * (5059 - 4753 * delta) * S1z +
                                           35 * (5059 + 4753 * delta) * S2z))) +
                  14 * total_mass2 *
                      (complex(0, 30) * rDOT *
                           (15280 * eta3 * (S1z + S2z) -
                            4 * eta2 *
                                (complex(0, -111) + 3562 * S1z + 682 * delta * S1z +
                                 945 * kappa1 * S1z3 +
                                 (3562 - 682 * delta +
                                  945 * (-2 + kappa1) * S1z2) *
                                     S2z +
                                 945 * (-2 + kappa2) * S1z * S2z2 +
                                 945 * kappa2 * S2z3) +
                            eta * (complex(0, -27520) +
                                  378 *
                                      (5 * (1 + delta) * kappa1 +
                                       3 * (3 + delta) * lambda1) *
                                      S1z3 +
                                  (29749 - 13605 * delta) * S2z -
                                  1512 * (1 + delta) * kappa1 * S1z2 * S2z -
                                  378 *
                                      (5 * (-1 + delta) * kappa2 +
                                       3 * (-3 + delta) * lambda2) *
                                      S2z3 +
                                  S1z * (29749 + 13605 * delta +
                                         1512 * (-1 + delta) * kappa2 *
                                             S2z2)) +
                            2 * (-8009 * (1 + delta) * S1z -
                                 567 * (1 + delta) * lambda1 * S1z3 +
                                 (-1 + delta) * S2z *
                                     (8009 + 567 * lambda2 * S2z2))) +
                       PhiDOT * r *
                           (285840 * eta3 * (S1z + S2z) -
                            30 * eta2 *
                                (complex(0, 11504) + 1890 * kappa1 * S1z3 +
                                 23090 * S2z - 1823 * delta * S2z +
                                 1890 * (-2 + kappa1) * S1z2 * S2z +
                                 1890 * kappa2 * S2z3 +
                                 S1z * (23090 + 1823 * delta +
                                        1890 * (-2 + kappa2) * S2z2)) +
                            30 * (689 * (1 + delta) * S1z -
                                  1134 * (1 + delta) * lambda1 * S1z3 +
                                  (-1 + delta) * S2z *
                                      (-689 + 1134 * lambda2 * S2z2)) +
                            2 * eta *
                                (complex(0, 415432) - 66840 * S2z +
                                 15 * (8 * (-557 + 1342 * delta) * S1z +
                                       189 *
                                           (5 * (1 + delta) * kappa1 +
                                            6 * (3 + delta) * lambda1) *
                                           S1z3 -
                                       10736 * delta * S2z -
                                       2457 * (1 + delta) * kappa1 * S1z2 *
                                           S2z +
                                       2457 * (-1 + delta) * kappa2 * S1z *
                                           S2z2 -
                                       189 *
                                           (5 * (-1 + delta) * kappa2 +
                                            6 * (-3 + delta) * lambda2) *
                                           S2z3)))))) /
                (317520. * r4))
            # + Henry et al. QC spinning hereditary terms
            # (2 * M_PI * (kappa1 * (1 + delta - 2 * eta) * S1z2 + S2z * (4 * eta * S1z - kappa2 * (-1 + delta + 2 * eta) * S2z)) * params.x4p5)
        )

    else:
        return complex(0, 0)

@register_hqc_lm(2, 2)    
def hQC_2_m_2(total_mass: float, eta: float, vpnorder: int, x: float, S1z: float, S2z: float,
              params: CommonVars) -> complex:
    
    EulerGamma = 0.5772156649015329
    x0 = 0.4375079683656479
    # b0 = 2 * total_mass / np.exp(0.5)
    # r0 = b0
    kappa1 = 1.0 # for black holes kappa and lambda is 1
    kappa2 = 1.0
    # delta = np.sqrt(1 - 4 * eta)
    xp5 = params.xp5
    logx = params.logx
    logb0 = params.logb0
    logr0 = params.logr0
    delta = params.delta
    # b0 = params.b0
    # r0 = params.r0
    x2p5 = x * x * xp5
    x3p5 = x * x2p5
    x4p5 = x * x3p5
    x4 = x*x*x*x
    x5 = x*x4

    eta2 = eta * eta

    S1z2 = S1z * S1z
 
    S2z2 = S2z * S2z

    # keeping only the hereditary terms
    # if vpnorder == 0:
    #     return (2 * x)

    # elif vpnorder == 2:
    #     return ((-5.095238095238095 + (55 * eta) / 21.) * x2)

    if vpnorder == 3:
        # return (4 * M_PI * x2p5)
        return (complex(0, 0.6666666666666666) * x2p5 * (-11 + 12 * EulerGamma - complex(0, 6) * M_PI + 24 * LOG2 + 12 * logb0 + 18 * logx))

    # elif vpnorder == 4:
    #     return ((-2.874338624338624 - (1069 * eta) / 108. + (2047 * eta2) / 756.) * x3)

    elif vpnorder == 5:
        # return ((complex(0,-48)*eta - (214*M_PI)/21. + (68*eta*M_PI)/21.)*x3p5)
        # return ((2 * (-107 + 34 * eta) * M_PI * x3p5) / 21.) # This is old implementation
        return (complex(0, 0.015873015873015872) * (-107 + 34 * eta) * x3p5 * (-11 + 12 * EulerGamma - complex(0, 6) * M_PI + 24 * LOG2 + 12 * logb0 + 18 * logx))

    elif vpnorder == 6:
        # return ((x4 * (-27392 * EulerGamma + M_PI * (complex(0, 13696) + 35 * (64 + 41 * eta) * M_PI) - 13696 * np.log(16 * x))) / 1680.) # This is the old implementation.

        return (x4 * (-45.216145124716554 + (456 * EulerGamma) / 35. - 16 * EulerGamma * EulerGamma - 
                complex(0, 6.514285714285714) * M_PI + complex(0, 16) * EulerGamma * M_PI + (4 * M_PI2) / 3. + 
                complex(0, 4.888888888888889) * S1z - complex(0, 5.333333333333333) * EulerGamma * S1z - 
                complex(0, 4.888888888888889) * eta * S1z + complex(0, 5.333333333333333) * EulerGamma * eta * S1z - 
                (8 * M_PI * S1z) / 3. + (8 * eta * M_PI * S1z) / 3. + complex(0, 4.888888888888889) * S2z - 
                complex(0, 5.333333333333333) * EulerGamma * S2z - complex(0, 4.888888888888889) * eta * S2z + 
                complex(0, 5.333333333333333) * EulerGamma * eta * S2z - (8 * M_PI * S2z) / 3. + (8 * eta * M_PI * S2z) / 3. + 
                (912 * LOG2) / 35. - 64 * EulerGamma * LOG2 + complex(0, 32) * M_PI * LOG2 - 
                complex(0, 10.666666666666666) * S1z * LOG2 + complex(0, 10.666666666666666) * eta * S1z * LOG2 - 
                complex(0, 10.666666666666666) * S2z * LOG2 + complex(0, 10.666666666666666) * eta * S2z * LOG2 - 
                64 * LOG2 * LOG2 + (88 * logb0) / 3. - 32 * EulerGamma * logb0 + 
                complex(0, 16) * M_PI * logb0 - complex(0, 5.333333333333333) * S1z * logb0 + 
                complex(0, 5.333333333333333) * eta * S1z * logb0 - complex(0, 5.333333333333333) * S2z * logb0 + 
                complex(0, 5.333333333333333) * eta * S2z * logb0 - 64 * LOG2 * logb0 - 
                16 * logb0 * logb0 - (1712 * logr0) / 105. + (684 * logx) / 35. - 
                48 * EulerGamma * logx + complex(0, 24) * M_PI * logx - complex(0, 8) * S1z * logx + 
                complex(0, 8) * eta * S1z * logx - complex(0, 8) * S2z * logx + complex(0, 8) * eta * S2z * logx - 
                96 * LOG2 * logx - 48 * logb0 * logx - 36 * logx * logx - 
                complex(0, 0.4444444444444444) * delta * (S1z - S2z) * (-11 + 12 * EulerGamma - complex(0, 6) * M_PI + 
                24 * LOG2 + 12 * logb0 + 18 * logx)))

    elif vpnorder == 7:
        # return (((-2173 - 4990 * eta + 1120 * eta2) * M_PI * x4p5) / 378.) # This is the old implementation.
        return (complex(0, 0.0004409171075837742) * x4p5 * (23903 - 26076 * EulerGamma + 57452 * eta - 
                59880 * EulerGamma * eta - 18728 * eta2 + 13440 * EulerGamma * eta2 + 
                complex(0, 13038) * M_PI + complex(0, 29940) * eta * M_PI - complex(0, 6720) * eta2 * M_PI - 
                8316 * kappa1 * S1z2 - 8316 * delta * kappa1 * S1z2 + 9072 * EulerGamma * kappa1 * S1z2 + 
                9072 * delta * EulerGamma * kappa1 * S1z2 + 16632 * kappa1 * eta * S1z2 - 
                18144 * EulerGamma * kappa1 * eta * S1z2 - complex(0, 4536) * kappa1 * M_PI * S1z2 - 
                complex(0, 4536) * delta * kappa1 * M_PI * S1z2 + complex(0, 9072) * kappa1 * eta * M_PI * S1z2 - 
                33264 * eta * S1z * S2z + 36288 * EulerGamma * eta * S1z * S2z - complex(0, 18144) * eta * M_PI * S1z * S2z - 
                8316 * kappa2 * S2z2 + 8316 * delta * kappa2 * S2z2 + 9072 * EulerGamma * kappa2 * S2z2 - 
                9072 * delta * EulerGamma * kappa2 * S2z2 + 16632 * kappa2 * eta * S2z2 - 
                18144 * EulerGamma * kappa2 * eta * S2z2 - complex(0, 4536) * kappa2 * M_PI * S2z2 + 
                complex(0, 4536) * delta * kappa2 * M_PI * S2z2 + complex(0, 9072) * kappa2 * eta * M_PI * S2z2 - 
                52152 * LOG2 - 119760 * eta * LOG2 + 26880 * eta2 * LOG2 + 
                18144 * kappa1 * S1z2 * LOG2 + 18144 * delta * kappa1 * S1z2 * LOG2 - 
                36288 * kappa1 * eta * S1z2 * LOG2 + 72576 * eta * S1z * S2z * LOG2 + 
                18144 * kappa2 * S2z2 * LOG2 - 18144 * delta * kappa2 * S2z2 * LOG2 - 
                36288 * kappa2 * eta * S2z2 * LOG2 + 12 * (-2173 - 4990 * eta + 1120 * eta2 + 
                756 * kappa1 * (1 + delta - 2 * eta) * S1z2 + 3024 * eta * S1z * S2z - 
                756 * kappa2 * (-1 + delta + 2 * eta) * S2z2) * logb0 + 18 * (-2173 - 4990 * eta + 
                1120 * eta2 + 756 * kappa1 * (1 + delta - 2 * eta) * S1z2 + 3024 * eta * S1z * S2z - 
                756 * kappa2 * (-1 + delta + 2 * eta) * S2z2) * logx))

    # 4PN non-spinning quasi-circular (2,2) mode has been obtained from Blanchet
    # et al. arXiv:2304.11185. This is written in terms of phi. Please refer to file shared by Quentin.

    elif vpnorder == 8:
        return (-1.3109423540262542e-11 * (x5 * (29059430400 * EulerGamma * EulerGamma * (-107 + 13 * eta) + 
                12 * (628830397253 + 1854914893791 * eta + 421984442880 * LOG2) + 415134720 * EulerGamma * (6099 - 40277 * eta + complex(0, 70) * (107 - 13 * eta) * M_PI + 280 * (-107 + 13 * eta) * LOG2) + 
                35 * (-28 * eta2 * (5385456111 + 5 * eta * (-163158374 + 26251249 * eta)) + 
                135135 * (54784 + 5 * eta * (1951 + 6560 * eta)) * M_PI2 - 955450349568 * eta * LOG2 + 
                3321077760 * (-107 + 13 * eta) * LOG2 * LOG2 - complex(0, 5930496) * M_PI * (6099 - 74773 * eta + 280 * (-107 + 13 * eta) * LOG2)) - 5700491596800 * np.log(total_mass) + 
                6966444925440 * logx + 69189120 * (420 * (-107 + 13 * eta) * logb0 * logb0 + 
                420 * (-107 + 13 * eta) * np.log(total_mass) * np.log(total_mass) - 14 * np.log(total_mass) * (-11803 * eta + 
                60 * EulerGamma * (-107 + 13 * eta) + complex(0, 30) * (107 - 13 * eta) * M_PI + 
                120 * (-107 + 13 * eta) * LOG2 + 90 * (-107 + 13 * eta) * logx) + 14 * logb0 * (5885 - 6420 * EulerGamma - 11803 * eta + 780 * EulerGamma * eta + complex(0, 3210) * M_PI - 
                complex(0, 390) * eta * M_PI + 120 * (-107 + 13 * eta) * LOG2 + (6420 - 780 * eta) * np.log(total_mass) + 
                90 * (-107 + 13 * eta) * logx) + 5 * logx * (-74149 * eta + 252 * EulerGamma * (-107 + 13 * eta) - 
                complex(0, 126) * (-107 + 13 * eta) * M_PI + 504 * (-107 + 13 * eta) * LOG2 + 
                189 * (-107 + 13 * eta) * logx) + 84672 * eta * np.log(x0)))))
        
        # return ((x5 * (276756480 * EulerGamma * (11449 + 19105 * eta) - 12 * (846557506853 + 1008017482431 * eta) + 35 * (28 * eta2 * (5385456111 + 5 * eta * (-163158374 + 26251249 * eta)) - complex(0, 3953664) * (11449 + 109657 * eta) * M_PI - 135135 * (54784 + 5 * eta * (1951 + 6560 * eta)) * M_PI2) + 138378240 * (11449 + 19105 * eta) * np.log(16 * x))) / 7.62810048e10)

    else:
        return complex(0, 0)

# H21

@register_hgo_lm(2, 1)
def hGO_2_m_1(
                total_mass: float,
                eta: float,
                r: float,
                rDOT: float,
                PhiDOT: float,
                vpnorder: int,
                S1z: float,
                S2z: float,
                x: float,
                params: CommonVars
            ) -> complex:
    
    delta = params.delta
    kappa1 = 1.0
    kappa2 = 1.0
    # r0 = 2 * total_mass / np.exp(0.5)
    r0 = params.r0

    rDOT2 = rDOT * rDOT
    rDOT3 = rDOT2 * rDOT
    rDOT4 = rDOT3 * rDOT
    rDOT5 = rDOT4 * rDOT
    rDOT6 = rDOT5 * rDOT

    total_mass2 = total_mass * total_mass
    total_mass3 = total_mass2 * total_mass

    PhiDOT2 = PhiDOT * PhiDOT
    PhiDOT3 = PhiDOT2 * PhiDOT
    PhiDOT4 = PhiDOT3 * PhiDOT
    PhiDOT5 = PhiDOT4 * PhiDOT
    PhiDOT6 = PhiDOT5 * PhiDOT

    r2 = r * r
    r3 = r2 * r
    r4 = r3 * r
    r5 = r4 * r
    r6 = r5 * r
    
    eta2 = eta * eta
    eta3 = eta2 * eta

    S1z2 = S1z * S1z
    S1z3 = S1z2 * S1z
 
    S2z2 = S2z * S2z

    if vpnorder == 1:
        return 0.6666666666666666j * delta * total_mass * PhiDOT

    elif vpnorder == 2:
        return (
            (-0.5j * total_mass2 *
             ((1 + delta) * S1z + (-1 + delta) * S2z))
            / r2
        )

    elif vpnorder == 3:
        return (
            (0.023809523809523808j * delta * total_mass * PhiDOT *
             (4 * total_mass * (-9 + 11 * eta) +
              r * ((19 - 24 * eta) * PhiDOT2 * r2 +
                   2j * (83 + 2 * eta) * PhiDOT * r * rDOT +
                   2 * (-33 + 10 * eta) * rDOT2)))
            / r
        )

    elif vpnorder == 4:
        return (
            (0.011904761904761904j * total_mass2 *
             (2 * total_mass *
              ((77 + 59 * eta + 11 * delta * (7 + eta)) * S1z +
               (-77 - 59 * eta + 11 * delta * (7 + eta)) * S2z) +
              r * (-2j * PhiDOT * r * rDOT *
                   (147 * (1 + delta) * S1z + (-83 + 13 * delta) * eta * S1z +
                    147 * (-1 + delta) * S2z + (83 + 13 * delta) * eta * S2z) +
                   rDOT2 *
                   ((105 * (1 + delta) - 4 * (13 + 15 * delta) * eta) * S1z +
                    (-105 + 15 * delta * (7 - 4 * eta) + 52 * eta) * S2z) +
                   4 * PhiDOT2 * r2 *
                   ((-21 - 21 * delta + 66 * eta + 4 * delta * eta) * S1z +
                    (21 - 21 * delta - 66 * eta + 4 * delta * eta) * S2z))))
            / r3
        )

    elif vpnorder == 5:
        term1 = (
            (0.0013227513227513227j * delta * total_mass * PhiDOT *
             (10 * total_mass2 * (31 - 205 * eta + 111 * eta2) -
              2 * total_mass * r *
              ((-197 + 5 * eta + 660 * eta2) * PhiDOT2 * r2 +
               1j * (-3167 - 5278 * eta + 201 * eta2) * PhiDOT * r * rDOT +
               8 * (202 + 587 * eta - 177 * eta2) * rDOT2) +
              3 * r2 *
              ((152 - 692 * eta + 333 * eta2) * PhiDOT4 * r4 +
               2j * (308 - 1607 * eta + 111 * eta2) * PhiDOT3 * r3 * rDOT -
               3 * (75 - 560 * eta + 68 * eta2) * PhiDOT2 * r2 * rDOT2 -
               2j * (-265 + 526 * eta + 18 * eta2) * PhiDOT * r * rDOT3 +
               (-241 + 550 * eta - 264 * eta2) * rDOT4)))
            / r2
        )

        # Henry et al. ecc spin terms
        term2 = (
            (total_mass3 *
             (-4 * rDOT *
              (kappa1 * (1 + delta - 2 * eta) * S1z2 +
               kappa2 * (-1 + delta + 2 * eta) * S2z2) -
              1j * PhiDOT * r *
              ((-((1 + delta) * (9 + kappa1)) +
                2 * (9 + (4 + 3 * delta) * kappa1) * eta) * S1z2 -
               12 * delta * eta * S1z * S2z +
               (9 + kappa2 - 2 * (9 + 4 * kappa2) * eta +
                delta * (-9 - kappa2 + 6 * kappa2 * eta)) * S2z2)))
            / (6.0 * r3)
        )

        return term1 + term2

    elif vpnorder == 6:
        term1 = (
            (delta * total_mass2 * eta * PhiDOT *
             (total_mass * (195 * PhiDOT * r - 946j * rDOT) +
              9 * r *
              (270 * PhiDOT3 * r3 -
               483j * PhiDOT2 * r2 * rDOT -
               580 * PhiDOT * r * rDOT2 +
               42j * rDOT3)))
            / (315.0 * r2)
        )

        # Henry et al. ecc spin terms
        term2 = (
            0.00033068783068783067j * total_mass2 *
            (3 * r2 *
             (4j * PhiDOT * r * rDOT3 *
              ((-315 * (1 + delta) + 2 * (251 + 463 * delta) * eta +
                4 * (15 + delta) * eta2) * S1z +
               (315 - 315 * delta - 502 * eta + 926 * delta * eta +
                4 * (-15 + delta) * eta2) * S2z) +
              12 * PhiDOT2 * r2 * rDOT2 *
              ((189 * (1 + delta) - 2 * (521 + 293 * delta) * eta +
                7 * (55 + 23 * delta) * eta2) * S1z +
               (189 * (-1 + delta) + 2 * (521 - 293 * delta) * eta +
                7 * (-55 + 23 * delta) * eta2) * S2z) +
              rDOT4 *
              ((567 * (1 + delta) - 16 * (77 + 64 * delta) * eta +
                8 * (177 + 173 * delta) * eta2) * S1z +
               (567 * (-1 + delta) + 16 * (77 - 64 * delta) * eta +
                8 * (-177 + 173 * delta) * eta2) * S2z) -
              4j * PhiDOT3 * r3 * rDOT *
              ((936 * (1 + delta) - 5 * (979 + 215 * delta) * eta +
                2 * (1353 + 293 * delta) * eta2) * S1z +
               (936 * (-1 + delta) + 5 * (979 - 215 * delta) * eta +
                2 * (-1353 + 293 * delta) * eta2) * S2z) +
              4 * PhiDOT4 * r4 *
              ((-252 * (1 + delta) + (1315 + 857 * delta) * eta +
                4 * (-285 + 43 * delta) * eta2) * S1z +
               (252 + 5 * eta * (-263 + 228 * eta) +
                delta * (-252 + eta * (857 + 172 * eta))) * S2z)) -
             2 * total_mass * r *
             (-1j * PhiDOT * r * rDOT *
              ((2043 * (1 + delta) + (37 + 2597 * delta) * eta +
                (10635 + 139 * delta) * eta2) * S1z +
               (2043 * (-1 + delta) + (-37 + 2597 * delta) * eta +
                (-10635 + 139 * delta) * eta2) * S2z) +
              PhiDOT2 * r2 *
              ((-765 - eta * (667 + 7773 * eta) +
                delta * (-765 + 7 * eta * (-533 + 245 * eta))) * S1z +
               (765 + eta * (667 + 7773 * eta) +
                delta * (-765 + 7 * eta * (-533 + 245 * eta))) * S2z) +
              4 * rDOT2 *
              ((-234 * (1 + delta) - 4 * (560 + 901 * delta) * eta +
                (483 + 1111 * delta) * eta2) * S1z +
               (234 + 7 * (320 - 69 * eta) * eta +
                delta * (-234 + eta * (-3604 + 1111 * eta))) * S2z)) +
             2 * total_mass2 *
             (1134 * kappa1 * (-1 - delta + (3 + delta) * eta) * S1z3 +
              1134 * (1 + delta) * (-2 + kappa1) * eta * S1z2 * S2z +
              S1z *
              (-5661 - 5661 * delta - 17156 * eta - 9172 * delta * eta +
               231 * eta2 + 775 * delta * eta2 +
               1134 * (-1 + delta) * (-2 + kappa2) * eta * S2z2) +
              S2z * (5661 - 5661 * delta + 17156 * eta - 9172 * delta * eta -
                     231 * eta2 + 775 * delta * eta2 +
                     1134 * kappa2 * (1 - delta + (-3 + delta) * eta) *
                     S2z2)))
            ) / r4
        

        return term1 + term2

    elif vpnorder == 7:

        spin_terms = (
            -23760 * total_mass3 *
            (PhiDOT3 * r4 *
             (-12j * (-35 + 107 * eta) * S1z +
              5 * (126 - 462 * eta + 80 * eta2 +
                   kappa1 * (-2 - 79 * eta + 153 * eta2)) * S1z2 +
              S2z * (-420j + 1284j * eta +
                     10 * (-63 + kappa2) * S2z +
                     5 * (462 + 79 * kappa2) * eta * S2z -
                     5 * (80 + 153 * kappa2) * eta2 * S2z)) -
             1j * PhiDOT2 * r3 * rDOT *
             (-24j * (-35 + 202 * eta) * S1z +
              5 * (-14 * eta * (3 + 58 * eta) +
                   kappa1 * (-409 + 1096 * eta + 6 * eta2)) * S1z2 +
              S2z * (-840j + 2045 * kappa2 * S2z +
                     10 * (406 - 3 * kappa2) * eta2 * S2z +
                     eta * (4848j + 210 * S2z - 5480 * kappa2 * S2z))) -
             4j * r * rDOT3 *
             (3j * (26 + 57 * eta) * S1z +
              5 * (4 * eta2 +
                   kappa1 * (-7 + 49 * eta + 15 * eta2)) * S1z2 -
              S2z * (78j - 35 * kappa2 * S2z +
                     5 * (4 + 15 * kappa2) * eta2 * S2z +
                     eta * (171j + 245 * kappa2 * S2z))) +
             2j * total_mass * rDOT *
             (-4j * (-78 + 769 * eta) * S1z +
              5 * ((14 - 59 * eta) * eta +
                   kappa1 * (-70 - 77 * eta + 18 * eta2)) * S1z2 +
              S2z * (-312j + 350 * kappa2 * S2z +
                     5 * (59 - 18 * kappa2) * eta2 * S2z +
                     eta * (3076j - 70 * S2z + 385 * kappa2 * S2z))) +
             PhiDOT * r2 * rDOT2 *
             (6j * (-62 + 219 * eta) * S1z -
              5 * (189 - 756 * eta + 388 * eta2 +
                   kappa1 * (98 - 266 * eta + 276 * eta2)) * S1z2 +
              S2z * (372j + 945 * S2z + 490 * kappa2 * S2z +
                     20 * (97 + 69 * kappa2) * eta2 * S2z -
                     2 * eta * (657j + 35 * (54 + 19 * kappa2) * S2z))) +
             total_mass * PhiDOT * r *
             (8j * (-61 + 480 * eta) * S1z +
              5 * (-392 + 448 * eta + 474 * eta2 +
                   kappa1 * (-11 + 150 * eta + 58 * eta2)) * S1z2 -
              S2z * (-488j - 5 * (392 + 11 * kappa2) * S2z +
                     10 * (237 + 29 * kappa2) * eta2 * S2z +
                     10 * eta * (384j + (224 + 75 * kappa2) * S2z))))
        )

        orbital_terms = (
            delta * total_mass *
            (-240 * total_mass * PhiDOT * r3 *
             ((197936 - 139360 * eta - 367105 * eta2 +
               253245 * eta3) * PhiDOT4 * r4 +
              1j * (279236 - 483940 * eta - 2817805 * eta2 +
                    459180 * eta3) * PhiDOT3 * r3 * rDOT -
              6 * (38627 + 89295 * eta - 492740 * eta2 +
                   75975 * eta3) * PhiDOT2 * r2 * rDOT2 -
              1j * (-731008 + 2287930 * eta + 981060 * eta2 +
                    10275 * eta3) * PhiDOT * r * rDOT3 +
              (-327667 + 436705 * eta + 659790 * eta2 -
               438255 * eta3) * rDOT4) +
             900 * PhiDOT * r4 *
             (2 * (-2594 + 27609 * eta - 74032 * eta2 +
                   25974 * eta3) * PhiDOT6 * r6 +
              4j * (-5730 + 58833 * eta - 137842 * eta2 +
                    17123 * eta3) * PhiDOT5 * r5 * rDOT +
              2 * (-114 - 41622 * eta + 147569 * eta2 +
                   4196 * eta3) * PhiDOT4 * r4 * rDOT2 +
              4j * (-9554 + 70788 * eta - 156227 * eta2 +
                    5810 * eta3) * PhiDOT3 * r3 * rDOT3 +
              (17619 - 138450 * eta + 322600 * eta2 -
               80816 * eta3) * PhiDOT2 * r2 * rDOT4 -
              2j * (8793 - 52230 * eta + 69340 * eta2 +
                    2536 * eta3) * PhiDOT * r * rDOT5 +
              2 * (3957 - 24534 * eta + 42584 * eta2 -
                   20800 * eta3) * rDOT6) -
             2 * total_mass3 *
             (-23760j * rDOT *
              (5 * (eta * (-14 + 31 * eta) + 7 * kappa1 * (10 + 31 * eta)) * S1z2 +
               2 * S1z * (-156j + 155 * eta2 * S2z +
                          2 * eta * (613j + 390 * S2z)) +
               S2z * (-312j + 350 * kappa2 * S2z +
                      155 * eta2 * S2z +
                      eta * (2452j + 35 * (-2 + 31 * kappa2) * S2z))) +
              PhiDOT * r *
              (8946400 * eta3 -
               8 * (6991786 + 724680j * S1z +
                    7425 * (392 + 11 * kappa1) * S1z2 +
                    724680j * S2z +
                    7425 * (392 + 11 * kappa2) * S2z2) -
               3600 * eta2 *
               (-628 + 33 * (-19 + 92 * kappa1) * S1z2 -
                7326 * S1z * S2z +
                33 * (-19 + 92 * kappa2) * S2z2) +
               15 * eta *
               (994455 * M_PI2 +
                8 * (-2249485 +
                     7920 * (-21 + 8 * kappa1) * S1z2 +
                     283536j * S2z +
                     7920 * (-21 + 8 * kappa2) * S2z2 -
                     1584 * S1z * (-179j + 170 * S2z))))) +
             3 * total_mass2 * r *
             (31680j * rDOT3 *
              (5 * (4 * eta2 + 7 * kappa1 * (-1 + 5 * eta)) * S1z2 +
               S1z * (78j + 40 * eta2 * S2z +
                      eta * (327j + 420 * S2z)) +
               S2z * (78j - 35 * kappa2 * S2z +
                      20 * eta2 * S2z +
                      eta * (327j + 175 * kappa2 * S2z))) -
              22 * PhiDOT * r * rDOT2 *
              (2553200 * eta3 -
               24 * (268267 + 5580j * S1z +
                     525 * (27 + 14 * kappa1) * S1z2 +
                     5580j * S2z +
                     525 * (27 + 14 * kappa2) * S2z2) -
               200 * eta2 *
               (39445 + 72 * (-4 + 21 * kappa1) * S1z2 -
                3600 * S1z * S2z +
                72 * (-4 + 21 * kappa2) * S2z2) +
               25 * eta *
               (23247 * M_PI2 +
                8 * (-69259 + 1026j * S1z +
                     126 * (27 + 5 * kappa1) * S1z2 +
                     1026j * S2z +
                     126 * (27 + 5 * kappa2) * S2z2))) +
              PhiDOT3 * r3 *
              (10071200 * eta3 +
               96 * (-421183 - 34650j * S1z +
                     825 * (-63 + kappa1) * S1z2 -
                     34650j * S2z +
                     825 * (-63 + kappa2) * S2z2) -
               400 * eta2 *
               (64177 + 792 * (-5 + 6 * kappa1) * S1z2 -
                17424 * S1z * S2z +
                792 * (-5 + 6 * kappa2) * S2z2) +
               15 * eta *
               (426195 * M_PI2 +
                8 * (-509635 +
                     330 * (210 + 83 * kappa1) * S1z2 +
                     29304j * S2z +
                     330 * (210 + 83 * kappa2) * S2z2 -
                     792 * S1z * (-37j + 70 * S2z)))) -
              2j * PhiDOT2 * r2 * rDOT *
              (-8330400 * eta3 +
               8 * (-2810116 - 415800j * S1z +
                    1012275 * kappa1 * S1z2 -
                    415800j * S2z +
                    1012275 * kappa2 * S2z2) +
               4800 * eta2 *
               (13411 + 33 * (19 + 12 * kappa1) * S1z2 +
                462 * S1z * S2z +
                33 * (19 + 12 * kappa2) * S2z2) +
               5 * eta *
               (1278585 * M_PI2 -
                8 * (5139685 +
                     990 * (-21 + 139 * kappa1) * S1z2 -
                     313632j * S2z +
                     990 * (-21 + 139 * kappa2) * S2z2 -
                     3564 * S1z * (88j + 185 * S2z))))))
        )

        log_term = (
            -13559040 * delta * total_mass3 * PhiDOT * r *
            (2 * total_mass - 3 * PhiDOT2 * r3 +
             6j * PhiDOT * r2 * rDOT +
             6 * r * rDOT2) *
            np.log(r / r0)
        )

        return (
            -1.0020843354176688e-7j *
            (spin_terms + orbital_terms + log_term)
        ) / r4

    else:
        return complex(0.0, 0.0)

@register_hqc_lm(2, 1)
def hQC_2_m_1(
                total_mass: float, 
                eta: float, 
                vpnorder: int, 
                x: float, 
                S1z: float, 
                S2z: float, 
                params:CommonVars
            ) -> complex:
    
    # delta: float = np.sqrt(1.0 - 4.0 * eta)
    EulerGamma: float = 0.5772156649015329
    # b0: float = 2.0 * total_mass / np.exp(0.5)
    # r0: float = b0

    delta = params.delta
    xp5 = params.xp5
    logx = params.logx
    logb0 = params.logb0
    logr0 = params.logr0

    x3p5 = x**3 * xp5
    x4p5 = x**4 * xp5
    x4 = x**4
    x3 = x**3

    # Common log terms to simplify the messy 4PN-7PN expressions
    log_term_base = (-7.0 + 6.0 * EulerGamma - 3j * M_PI + np.log(64.0) + 6.0 * logb0 + 9.0 * logx)

    if vpnorder == 4:
        return (-2.0 * delta * x3 * log_term_base) / 9.0

    elif vpnorder == 5:
        spin_factor = ((1.0 + delta) * S1z + (-1.0 + delta) * S2z)
        return (spin_factor * x3p5 * log_term_base) / 6.0

    elif vpnorder == 6:
        return -0.007936507936507936 * (delta * (-17.0 + 6.0 * eta) * x4 * log_term_base)

    elif vpnorder == 7:
        # Heavily nested 3.5PN (order 7) term
        term1 = (x4p5 * (
            -98 * S1z + 84 * EulerGamma * S1z + 3017 * eta * S1z - 2586 * EulerGamma * eta * S1z - 
            42j * M_PI * S1z + 1293j * eta * M_PI * S1z + 98 * S2z - 84 * EulerGamma * S2z - 
            3017 * eta * S2z + 2586 * EulerGamma * eta * S2z + 42j * M_PI * S2z - 1293j * eta * M_PI * S2z + 
            84 * S1z * LOG2 - 2586 * eta * S1z * LOG2 - 84 * S2z * LOG2 + 2586 * eta * S2z * LOG2 + 
            84 * S1z * logb0 - 2586 * eta * S1z * logb0 - 84 * S2z * logb0 + 2586 * eta * S2z * logb0 + 
            126 * S1z * logx - 3879 * eta * S1z * logx - 126 * S2z * logx + 3879 * eta * S2z * logx + 
            252 * delta * (
                1.5875434618291762j + 1.7523809523809524j * EulerGamma - 1.3333333333333333j * EulerGamma**2 + 
                (92 * M_PI) / 105.0 - (4 * EulerGamma * M_PI) / 3.0 + 0.1111111111111111j * M_PI**2 - 
                (7 * S1z) / 18.0 + (EulerGamma * S1z) / 3.0 + (29 * eta * S1z) / 12.0 - (29 * EulerGamma * eta * S1z) / 14.0 - 
                0.16666666666666666j * M_PI * S1z + 1.0357142857142858j * eta * M_PI * S1z - 
                (7 * S2z) / 18.0 + (EulerGamma * S2z) / 3.0 + (29 * eta * S2z) / 12.0 - (29 * EulerGamma * eta * S2z) / 14.0 - 
                0.16666666666666666j * M_PI * S2z + 1.0357142857142858j * eta * M_PI * S2z + 
                1.7523809523809524j * LOG2 - 2.6666666666666665j * EulerGamma * LOG2 - 
                (4 * M_PI * LOG2) / 3.0 + (S1z * LOG2) / 3.0 - (29 * eta * S1z * LOG2) / 14.0 + 
                (S2z * LOG2) / 3.0 - (29 * eta * S2z * LOG2) / 14.0 - 1.3333333333333333j * LOG2**2 - 
                1.3333333333333333j * logb0**2 - 1.3587301587301588j * logr0 + 
                (logb0 * (392j - 336j * EulerGamma - 168 * M_PI + 42 * S1z - 261 * eta * S1z + 42 * S2z - 261 * eta * S2z - 336j * LOG2 - 504j * logx)) / 126.0 + 
                2.6285714285714286j * logx - 4j * EulerGamma * logx - 2 * M_PI * logx + 
                (S1z * logx) / 2.0 - (87 * eta * S1z * logx) / 28.0 + (S2z * logx) / 2.0 - 
                (87 * eta * S2z * logx) / 28.0 - 4j * LOG2 * logx - 3j * logx**2
            )
        )) / 252.0
        return term1

    else:
        return complex(0.0, 0.0)


# H33

# def hGO_3_m_3(
#                 total_mass: float,
#                 eta: float,
#                 r: float,
#                 rDOT: float,
#                 PhiDOT: float,
#                 vpnorder: int,
#                 S1z: float,
#                 S2z: float,
#                 x: float,
#             ) -> complex:
#     delta = np.sqrt(1 - 4 * eta)
#     kappa1 = 1.0
#     kappa2 = 1.0
#     r0 = 2 * total_mass / np.exp(0.5)

#     combination_a  = 1j * PhiDOT * r - rDOT
#     combination_a3 = combination_a ** 3

#     combination_b  = PhiDOT * r + 1j * rDOT
#     combination_b2 = combination_b ** 2
#     combination_b4 = combination_b2 ** 2
#     combination_b6 = combination_b ** 6

#     combination_c  = PhiDOT * r - 1j * rDOT
#     combination_c2 = combination_c ** 2

#     combination_d  = -1j * PhiDOT * r + rDOT
#     combination_d5 = combination_d ** 5

#     combination_e  = 1j * PhiDOT * r + rDOT
#     combination_e3 = combination_e ** 3

#     rDOT2 = rDOT * rDOT
#     rDOT3 = rDOT2 * rDOT
#     rDOT4 = rDOT3 * rDOT
#     rDOT5 = rDOT4 * rDOT
#     rDOT6 = rDOT5 * rDOT
#     rDOT7 = rDOT6 * rDOT

#     total_mass2 = total_mass * total_mass
#     total_mass3 = total_mass2 * total_mass
#     total_mass4 = total_mass3 * total_mass

#     PhiDOT2 = PhiDOT * PhiDOT
#     PhiDOT3 = PhiDOT2 * PhiDOT
#     PhiDOT4 = PhiDOT3 * PhiDOT
#     PhiDOT5 = PhiDOT4 * PhiDOT
#     PhiDOT6 = PhiDOT5 * PhiDOT
#     PhiDOT7 = PhiDOT6 * PhiDOT

#     r2 = r * r
#     r3 = r2 * r
#     r4 = r3 * r
#     r5 = r4 * r
#     r6 = r5 * r
#     r7 = r6 * r
    
#     eta2 = eta * eta
#     eta3 = eta2 * eta

#     S1z2 = S1z * S1z
 
#     S2z2 = S2z * S2z


#     if vpnorder == 1:
#         return (
#             (np.sqrt(0.11904761904761904) * delta *
#              (2 * r * combination_a3 +
#               total_mass * (-7j * PhiDOT * r + 4 * rDOT)))
#             / (2.0 * r)
#         )

#     elif vpnorder == 3:
#         return (
#             (np.sqrt(0.11904761904761904) * delta *
#              (6 * (-5 + 19 * eta) * r2 * combination_b4 *
#               (1j * PhiDOT * r + rDOT) +
#               2 * total_mass2 *
#               (-3j * (-101 + 43 * eta) * PhiDOT * r +
#                (-109 + 86 * eta) * rDOT) +
#               3 * total_mass * r *
#               (-12j * (1 + 4 * eta) * PhiDOT3 * r3 +
#                6 * (14 + 31 * eta) * PhiDOT2 * r2 * rDOT +
#                3j * (33 + 62 * eta) * PhiDOT * r * rDOT2 -
#                4 * (8 + 17 * eta) * rDOT3)))
#             / (36.0 * r2)
#         )

#     elif vpnorder == 4:
#         return (
#             (-0.125j * np.sqrt(0.11904761904761904) * total_mass2 *
#              (4 * total_mass * (-1 + 5 * eta) * ((1 + delta) * S1z + (-1 + delta) * S2z) +
#               r * (2 * rDOT2 *
#                    (6 * (1 + delta) * S1z - 5 * (5 + 3 * delta) * eta * S1z +
#                     (-6 + delta * (6 - 15 * eta) + 25 * eta) * S2z) +
#                    PhiDOT2 * r2 *
#                    (-24 * (1 + delta) * S1z + (119 + 33 * delta) * eta * S1z +
#                     (24 - 119 * eta + 3 * delta * (-8 + 11 * eta)) * S2z) +
#                    2j * PhiDOT * r * rDOT *
#                    (-18 * (1 + delta) * S1z + (77 + 39 * delta) * eta * S1z +
#                     (18 - 77 * eta + 3 * delta * (-6 + 13 * eta)) * S2z))))
#             / r3
#         )

#     elif vpnorder == 5:
#         term1 = (
#             delta *
#             (30 * (183 - 1579 * eta + 3387 * eta2) * r3 *
#              combination_c2 * combination_d5 +
#              10 * total_mass3 *
#              (-1j * (26473 - 27451 * eta + 9921 * eta2) * PhiDOT * r +
#               4 * (623 - 732 * eta + 1913 * eta2) * rDOT) +
#              2 * total_mass2 * r *
#              (-11j * (-5353 - 13493 * eta + 4671 * eta2) * PhiDOT3 * r3 +
#               (-75243 - 142713 * eta + 192821 * eta2) * PhiDOT2 * r2 * rDOT +
#               220j * (-256 + 781 * eta + 840 * eta2) * PhiDOT * r * rDOT2 -
#               10 * (-756 + 8238 * eta + 7357 * eta2) * rDOT3) +
#              3 * total_mass * r2 *
#              (2j * (-7633 + 9137 * eta + 28911 * eta2) * PhiDOT5 * r5 -
#               4 * (-8149 + 1576 * eta + 43533 * eta2) * PhiDOT4 * r4 * rDOT -
#               2j * (-9297 - 19517 * eta + 64839 * eta2) * PhiDOT3 * r3 * rDOT2 -
#               32 * (-1288 + 3667 * eta + 4056 * eta2) * PhiDOT2 * r2 * rDOT3 -
#               5j * (-9851 + 17954 * eta + 40968 * eta2) * PhiDOT * r * rDOT4 +
#               20 * (-771 + 1126 * eta + 3616 * eta2) * rDOT5))
#             / (1584.0 * np.sqrt(210) * r3)
#         )

#         # Henry et al. ecc spin terms
#         term2 = (
#             0.125j * np.sqrt(1.0714285714285714) * total_mass3 *
#             (7 * PhiDOT * r + 2j * rDOT) *
#             (kappa1 * (-1 - delta + 2 * (2 + delta) * eta) * S1z2 +
#              S2z * (-4 * delta * eta * S1z +
#                     kappa2 * (1 - delta + 2 * (-2 + delta) * eta) * S2z))
#             / r3
#         )

#         return term1 + term2

#     elif vpnorder == 6:
#         term1 = (
#             -(delta * total_mass2 * eta *
#               (668 * total_mass2 +
#                2 * total_mass * r *
#                (4081 * PhiDOT2 * r2 +
#                 297j * PhiDOT * r * rDOT - 452 * rDOT2) +
#                5 * r2 *
#                (1329 * PhiDOT4 * r4 -
#                 2926j * PhiDOT3 * r3 * rDOT -
#                 384 * PhiDOT2 * r2 * rDOT2 -
#                 408j * PhiDOT * r * rDOT3 +
#                 200 * rDOT4)))
#             / (36.0 * np.sqrt(210) * r4)
#         )

#         # Henry et al. ecc spin terms
#         term2 = (
#             -0.006944444444444444j * total_mass2 *
#             (10 * total_mass2 *
#              ((252 * (1 + delta) - (1277 + 1279 * delta) * eta +
#                8 * (12 + 47 * delta) * eta2) * S1z +
#               (252 * (-1 + delta) + (1277 - 1279 * delta) * eta +
#                8 * (-12 + 47 * delta) * eta2) * S2z) +
#              2 * total_mass * r *
#              (2 * PhiDOT2 * r2 *
#               ((1320 * (1 + delta) - 2 * (4469 + 211 * delta) * eta +
#                 (8709 + 2777 * delta) * eta2) * S1z +
#                (1320 * (-1 + delta) + 8938 * eta - 422 * delta * eta +
#                 (-8709 + 2777 * delta) * eta2) * S2z) +
#               3j * PhiDOT * r * rDOT *
#               ((2000 * (1 + delta) - (9147 + 3173 * delta) * eta +
#                 (8911 + 5273 * delta) * eta2) * S1z +
#                (2000 * (-1 + delta) + (9147 - 3173 * delta) * eta +
#                 (-8911 + 5273 * delta) * eta2) * S2z) +
#               10 * rDOT2 *
#               ((-105 * (1 + delta) + (541 + 77 * delta) * eta -
#                 2 * (462 + 247 * delta) * eta2) * S1z +
#                (105 + eta * (-541 + 924 * eta) +
#                 delta * (-105 + (77 - 494 * eta) * eta)) * S2z)) -
#              3 * r2 *
#              (-3 * PhiDOT2 * r2 * rDOT2 *
#               ((480 * (1 + delta) - (1711 + 1889 * delta) * eta +
#                 2 * (-1161 + 757 * delta) * eta2) * S1z +
#                (480 * (-1 + delta) + (1711 - 1889 * delta) * eta +
#                 2 * (1161 + 757 * delta) * eta2) * S2z) +
#               2 * PhiDOT4 * r4 *
#               ((350 * (1 + delta) - 4 * (404 + 461 * delta) * eta +
#                 (883 + 769 * delta) * eta2) * S1z +
#                (350 * (-1 + delta) + 4 * (404 - 461 * delta) * eta +
#                 (-883 + 769 * delta) * eta2) * S2z) +
#               2j * PhiDOT3 * r3 * rDOT *
#               ((660 * (1 + delta) - (4061 + 2899 * delta) * eta +
#                 (2643 + 4789 * delta) * eta2) * S1z +
#                (660 * (-1 + delta) + (4061 - 2899 * delta) * eta +
#                 (-2643 + 4789 * delta) * eta2) * S2z) +
#               10 * rDOT4 *
#               ((-30 * (1 + delta) + (187 + 101 * delta) * eta -
#                 2 * (159 + 61 * delta) * eta2) * S1z +
#                (30 + eta * (-187 + 318 * eta) +
#                 delta * (-30 + (101 - 122 * eta) * eta)) * S2z) +
#               2j * PhiDOT * r * rDOT3 *
#               ((90 + eta * (-1321 + 5118 * eta) +
#                 delta * (90 + eta * (-319 + 714 * eta))) * S1z +
#                (-90 + (1321 - 5118 * eta) * eta +
#                 delta * (90 + eta * (-319 + 714 * eta))) * S2z)))
#             / (np.sqrt(210) * r4)
#         )

#         return term1 + term2

#     elif vpnorder == 7:
#         M_PI2 = M_PI ** 2

#         spin_terms = (
#             504504 * total_mass3 *
#             (2 * total_mass *
#              (rDOT *
#               (S1z * (108 - 498 * eta +
#                       5j * (-24 - 5 * kappa1 + 3 * (76 + kappa1) * eta +
#                             4 * (-111 + 13 * kappa1) * eta2) * S1z) +
#                6 * (-18 + 83 * eta) * S2z -
#                5j * (-24 - 5 * kappa2 + 3 * (76 + kappa2) * eta +
#                      4 * (-111 + 13 * kappa2) * eta2) * S2z2) +
#               PhiDOT * r *
#               (S1z * (-3j * (-99 + 581 * eta) +
#                       5 * (-24 + 399 * kappa1 + 48 * (7 - 19 * eta) * eta +
#                            kappa1 * eta * (-1629 + 188 * eta)) * S1z) +
#                3j * (-99 + 581 * eta) * S2z -
#                5 * (-24 + 399 * kappa2 + 48 * (7 - 19 * eta) * eta +
#                     kappa2 * eta * (-1629 + 188 * eta)) * S2z2)) +
#              r * (3 * PhiDOT * r * rDOT2 *
#                   (S1z * (216j + 545 * kappa1 * S1z +
#                           40 * (45 + 8 * kappa1) * eta2 * S1z -
#                           30 * eta * (50j + 20 * S1z + 73 * kappa1 * S1z)) +
#                    12j * (-18 + 125 * eta) * S2z -
#                    5 * (109 * kappa2 - 6 * (20 + 73 * kappa2) * eta +
#                         8 * (45 + 8 * kappa2) * eta2) * S2z2) +
#                   2 * rDOT3 *
#                   (S1z * (-54 + 145j * kappa1 * S1z -
#                           30j * eta * (11j + 2 * eta * S1z +
#                                       kappa1 * (17 + 8 * eta) * S1z)) +
#                    6 * (9 - 55 * eta) * S2z +
#                    5j * (12 * eta2 +
#                          kappa2 * (-29 + 6 * eta * (17 + 8 * eta))) * S2z2) +
#                   6 * PhiDOT2 * r2 * rDOT *
#                   (S1z * (297 - 465 * eta +
#                           5j * (6 * (20 - 87 * eta) * eta +
#                                 kappa1 * (-50 + 3 * eta * (47 + 76 * eta))) * S1z) +
#                    3 * (-99 + 155 * eta) * S2z -
#                    5j * (6 * (20 - 87 * eta) * eta +
#                          kappa2 * (-50 + 3 * eta * (47 + 76 * eta))) * S2z2) +
#                   PhiDOT3 * r3 *
#                   (3j * (531 - 1295 * eta) * S1z +
#                    10 * (-33 * kappa1 + 6 * (30 + 13 * kappa1) * eta +
#                          4 * (-96 + 67 * kappa1) * eta2) * S1z2 +
#                    S2z * (-1593j + 330 * kappa2 * S2z +
#                           5 * eta * (777j -
#                                     4 * (90 + 39 * kappa2 - 192 * eta +
#                                          134 * kappa2 * eta) * S2z)))))
#         )

#         orbital_terms = (
#             delta *
#             (-17640 * (-4083 + eta * (58311 + eta * (-269240 + 405617 * eta))) *
#              r4 * combination_b6 * combination_e3 +
#              168 * total_mass2 * r2 *
#              (1j * (-7508635 + 7 * eta * (-1318438 + eta * (-10231834 + 9667755 * eta))) *
#               PhiDOT5 * r5 +
#               7 * (1235591 + eta * (884445 + (23935218 - 26913443 * eta) * eta)) *
#               PhiDOT4 * r4 * rDOT -
#               1j * (8961149 + 7 * eta * (-31755709 + eta * (-11134798 + 22187331 * eta))) *
#               PhiDOT3 * r3 * rDOT2 -
#               (-36806435 + 7 * eta * (33178545 + eta * (24565078 + 22873537 * eta))) *
#               PhiDOT2 * r2 * rDOT3 -
#               5j * (-7761899 + 7 * eta * (2892563 + 5998602 * eta + 7493619 * eta2)) *
#               PhiDOT * r * rDOT4 +
#               5 * (-2422057 + 7 * eta * (501045 + eta * (2033141 + 2771816 * eta))) *
#               rDOT5) +
#              1764 * total_mass * r3 *
#              (-2j * (239087 + eta * (-1206515 + eta * (422631 + 3979375 * eta))) *
#               PhiDOT7 * r7 +
#               2 * (621284 + eta * (-2279907 + 2 * eta * (-1180187 + 5876531 * eta))) *
#               PhiDOT6 * r6 * rDOT +
#               2j * (39270 + eta * (1235486 - 5319747 * eta + 4406349 * eta2)) *
#               PhiDOT5 * r5 * rDOT2 +
#               8 * (349111 + 4 * eta * (-519370 + 33 * eta * (10939 + 42635 * eta))) *
#               PhiDOT4 * r4 * rDOT3 +
#               2j * (1212607 + 3 * eta * (-2012698 - 67827 * eta + 7955628 * eta2)) *
#               PhiDOT3 * r3 * rDOT4 +
#               4 * (201135 + 2 * eta * (-773107 + eta * (1214819 + 1157652 * eta))) *
#               PhiDOT2 * r2 * rDOT5 +
#               5j * (333969 + 2 * eta * (-981471 + 4 * eta * (154039 + 750016 * eta))) *
#               PhiDOT * r * rDOT6 -
#               40 * (13245 + 2 * eta * (-37005 + eta * (14251 + 130160 * eta))) *
#               rDOT7) +
#              2 * total_mass4 *
#              (4 * rDOT *
#               (269279500 * eta3 +
#                2 * (-174108226 +
#                     63063 * S1z * (108j + 5 * (24 + 5 * kappa1) * S1z) +
#                     63063 * S2z * (108j + 5 * (24 + 5 * kappa2) * S2z)) -
#                21 * eta *
#                (103100846 + 1846845 * M_PI2 + 1693692j * S2z -
#                 6006 * (S1z * (-282j + 5 * (-180 + 7 * kappa1) * S1z) -
#                         980 * S1z * S2z +
#                         5 * (-180 + 7 * kappa2) * S2z2)) -
#                2940 * eta2 *
#                (-122855 +
#                 4719 * ((-6 + 7 * kappa1) * S1z2 -
#                         26 * S1z * S2z +
#                         (-6 + 7 * kappa2) * S2z2))) +
#               1j * PhiDOT * r *
#               (-1176172480 * eta3 +
#                8 * (-74084729 +
#                     189189 * S1z * (99j + 5 * (-8 + 133 * kappa1) * S1z) +
#                     189189 * S2z * (99j + 5 * (-8 + 133 * kappa2) * S2z)) +
#                35280 * eta2 *
#                (56255 +
#                 429 * ((-22 + 65 * kappa1) * S1z2 -
#                        174 * S1z * S2z +
#                        (-22 + 65 * kappa2) * S2z2)) -
#                147 * eta *
#                (-65012788 + 4485195 * M_PI2 + 3943368j * S2z +
#                 10296 * (S1z * (383j + 5 * (-96 + 277 * kappa1) * S1z) -
#                          3220 * S1z * S2z +
#                          5 * (-96 + 277 * kappa2) * S2z2)))) +
#              total_mass3 * r *
#              (-12j * PhiDOT * r * rDOT2 *
#               (-1035895280 * eta3 -
#                2 * (-547993687 +
#                     63063 * S1z * (216j + 545 * kappa1 * S1z) +
#                     63063 * S2z * (216j + 545 * kappa2 * S2z)) +
#                77 * eta *
#                (42451610 + 1511055 * M_PI2 + 1749384j * S2z +
#                 6552 * (S1z * (267j + 25 * (6 + 11 * kappa1) * S1z) -
#                         5 * S1z * S2z +
#                         25 * (6 + 11 * kappa2) * S2z2)) +
#                490 * eta2 *
#                (-5802767 +
#                 5148 * ((-6 + 23 * kappa1) * S1z2 -
#                         58 * S1z * S2z +
#                         (-6 + 23 * kappa2) * S2z2))) +
#               4 * rDOT3 *
#               (-1359334480 * eta3 -
#                4 * (-150254558 +
#                     63063 * S1z * (54j + 145 * kappa1 * S1z) +
#                     63063 * S2z * (54j + 145 * kappa2 * S2z)) +
#                231 * eta *
#                (8490448 + 503685 * M_PI2 + 242424j * S2z +
#                 2184 * (S1z * (111j + 110 * kappa1 * S1z) +
#                         70 * S1z * S2z +
#                         110 * kappa2 * S2z2)) +
#                11760 * eta2 *
#                (-312980 +
#                 429 * ((3 + 25 * kappa1) * S1z2 -
#                        44 * S1z * S2z +
#                        (3 + 25 * kappa2) * S2z2))) +
#               6 * PhiDOT2 * r2 * rDOT *
#               (2368900688 * eta3 +
#                8 * (-812986529 +
#                     63063 * S1z * (297j + 250 * kappa1 * S1z) +
#                     63063 * S2z * (297j + 250 * kappa2 * S2z)) -
#                1176 * eta2 *
#                (2423171 +
#                 4290 * ((-3 + 41 * kappa1) * S1z2 -
#                         88 * S1z * S2z +
#                         (-3 + 41 * kappa2) * S2z2)) +
#                539 * eta *
#                (-24139772 + 647595 * M_PI2 + 120744j * S2z -
#                 936 * (S1z * (-129j + 600 * S1z + 205 * kappa1 * S1z) -
#                        460 * S1z * S2z +
#                        5 * (120 + 41 * kappa2) * S2z2))) +
#               1j * PhiDOT3 * r3 *
#               (-4538040136 * eta3 -
#                88 * (259018351 +
#                      17199 * S1z * (-531j + 110 * kappa1 * S1z) +
#                      17199 * S2z * (-531j + 110 * kappa2 * S2z)) +
#                2352 * eta2 *
#                (7332973 +
#                 12870 * ((5 + 23 * kappa1) * S1z2 -
#                          36 * S1z * S2z +
#                          (5 + 23 * kappa2) * S2z2)) +
#                21 * eta *
#                (49864815 * M_PI2 +
#                 8 * (-88128538 - 2099097j * S2z +
#                      9009 * (S1z * (-233j + 40 * (15 + kappa1) * S1z) +
#                              360 * S1z * S2z +
#                              40 * (15 + kappa2) * S2z2))))))
#         )

#         log_term = (
#             74954880 * delta * total_mass3 *
#             (22j * total_mass * PhiDOT * r +
#              59j * PhiDOT3 * r4 + 8 * total_mass * rDOT +
#              66 * PhiDOT2 * r3 * rDOT +
#              24j * PhiDOT * r2 * rDOT2 -
#              4 * r * rDOT3) *
#             np.log(r / r0)
#         )

#         return (
#             1j * (spin_terms + orbital_terms + log_term)
#             / (2.4216192e7 * np.sqrt(210) * r4)
#         )

#     else:
#         return complex(0, 0)

# def hQC_3_m_3(
#                 total_mass: float, 
#                 eta: float, 
#                 vpnorder: int, 
#                 x: float, 
#                 S1z: float, 
#                 S2z: float, 
#                 params: CommonVars,
#             ) -> complex:
    
#     delta: float = np.sqrt(1.0 - 4.0 * eta)
#     EulerGamma: float = 0.5772156649015329
#     b0: float = 2.0 * total_mass / np.exp(0.5)
#     r0: float = b0

#     x3 = x**3
#     x4 = x**4
#     x4p5 = x4*xp5

#     logx = params.logx
#     logb0 = params.logb0
#     logr0 = params.logr0

#     if vpnorder == 4:
#         return (
#             (3.0 * np.sqrt(0.04285714285714286) * delta * x3 * (
#                 -97.0 + 60.0 * EulerGamma - 30.0j * M_PI + 
#                 60.0 * np.log(2.0) + 60.0 * np.log(3.0) + 
#                 60.0 * logb0 + 90.0 * logx
#             )) / 4.0
#         )

#     elif vpnorder == 6:
#         numerator_6 = (delta * x4 * (
#             188568.0 - 116640.0 * EulerGamma - 70537.0 * eta + 43740.0 * EulerGamma * eta + 
#             58320.0j * M_PI - 21870.0j * eta * M_PI - 
#             116640.0 * np.log(2.0) + 43740.0 * eta * np.log(2.0) - 
#             116640.0 * np.log(3.0) + 43740.0 * eta * np.log(3.0) + 
#             14580.0 * (-8.0 + 3.0 * eta) * logb0 + 
#             21870.0 * (-8.0 + 3.0 * eta) * logx
#         ))
#         return numerator_6 / (216.0 * np.sqrt(210.0))

#     elif vpnorder == 7:
#         # Logarithmic expansion for the 3.5PN (order 7) term
#         # log_b0_term handles the final nested log(b0) block
#         log_b0_term = 15876.0 * logb0 * (
#             -194.0j * delta + 120.0j * delta * EulerGamma + 60.0 * delta * M_PI + 
#             5.0 * (-4.0 - 4.0 * delta + 19.0 * eta + 5.0 * delta * eta) * S1z + 
#             20.0 * S2z - 20.0 * delta * S2z - 95.0 * eta * S2z + 25.0 * delta * eta * S2z + 
#             120.0j * delta * np.log(2.0) + 120.0j * delta * np.log(3.0) + 
#             180.0j * delta * logx
#         )

#         numerator_7 = (x4p5 * (
#             1434564.0j * delta - 2490264.0j * delta * EulerGamma + 952560.0j * delta * EulerGamma**2 - 
#             1245132.0 * delta * M_PI + 952560.0 * delta * EulerGamma * M_PI - 
#             79380.0j * delta * (M_PI**2) + 513324.0 * S1z + 513324.0 * delta * S1z - 
#             317520.0 * EulerGamma * S1z - 317520.0 * delta * EulerGamma * S1z - 
#             2439857.0 * eta * S1z - 643223.0 * delta * eta * S1z + 1508220.0 * EulerGamma * eta * S1z + 
#             396900.0 * delta * EulerGamma * eta * S1z + 158760.0j * M_PI * S1z + 158760.0j * delta * M_PI * S1z - 
#             754110.0j * eta * M_PI * S1z - 198450.0j * delta * eta * M_PI * S1z - 
#             513324.0 * S2z + 513324.0 * delta * S2z + 317520.0 * EulerGamma * S2z - 
#             317520.0 * delta * EulerGamma * S2z + 2439857.0 * eta * S2z - 643223.0 * delta * eta * S2z - 
#             1508220.0 * EulerGamma * eta * S2z + 396900.0 * delta * EulerGamma * eta * S2z - 
#             158760.0j * M_PI * S2z + 158760.0j * delta * M_PI * S2z + 
#             754110.0j * eta * M_PI * S2z - 198450.0j * delta * eta * M_PI * S2z - 
#             2490264.0j * delta * np.log(2.0) + 1905120.0j * delta * EulerGamma * np.log(2.0) + 
#             952560.0 * delta * M_PI * np.log(2.0) - 317520.0 * S1z * np.log(2.0) - 
#             317520.0 * delta * S1z * np.log(2.0) + 1508220.0 * eta * S1z * np.log(2.0) + 
#             396900.0 * delta * eta * S1z * np.log(2.0) + 317520.0 * S2z * np.log(2.0) - 
#             317520.0 * delta * S2z * np.log(2.0) - 1508220.0 * eta * S2z * np.log(2.0) + 
#             396900.0 * delta * eta * S2z * np.log(2.0) + 952560.0j * delta * (np.log(2.0)**2) - 
#             2490264.0j * delta * np.log(3.0) + 1905120.0j * delta * EulerGamma * np.log(3.0) + 
#             952560.0 * delta * M_PI * np.log(3.0) - 317520.0 * S1z * np.log(3.0) - 
#             317520.0 * delta * S1z * np.log(3.0) + 1508220.0 * eta * S1z * np.log(3.0) + 
#             396900.0 * delta * eta * S1z * np.log(3.0) + 317520.0 * S2z * np.log(3.0) - 
#             317520.0 * delta * S2z * np.log(3.0) - 1508220.0 * eta * S2z * np.log(3.0) + 
#             396900.0 * delta * eta * S2z * np.log(3.0) + 1905120.0j * delta * np.log(2.0) * np.log(3.0) + 
#             952560.0j * delta * (np.log(3.0)**2) + 952560.0j * delta * (logb0**2) + 
#             589680.0j * delta * logr0 - 3735396.0j * delta * logx + 
#             2857680.0j * delta * EulerGamma * logx + 1428840.0 * delta * M_PI * logx - 
#             476280.0 * S1z * logx - 476280.0 * delta * S1z * logx + 
#             2262330.0 * eta * S1z * logx + 595350.0 * delta * eta * S1z * logx + 
#             476280.0 * S2z * logx - 476280.0 * delta * S2z * logx - 
#             2262330.0 * eta * S2z * logx + 595350.0 * delta * eta * S2z * logx + 
#             2857680.0j * delta * np.log(2.0) * logx + 
#             2857680.0j * delta * np.log(3.0) * logx + 
#             2143260.0j * delta * (logx**2) + 
#             log_b0_term
#         ))
#         return numerator_7 / (2352.0 * np.sqrt(210.0))

#     else:
#         return complex(0, 0)

def generate_hlm(l: int, m: int, total_mass: float, eta: float, r: float, rDOT: float, Phi: float,
             PhiDOT: float, R: float, vpnorder: int, S1z: float,
             S2z: float, x: float, params: CommonVars) -> complex:
    if vpnorder < 0 or vpnorder > 8:
        raise ValueError(f"Error in hl_{l}_m_{m}: Input PN order parameter should be between [0, 8].")

    else:
        # Calculate the leading amplitude coefficient
        amplitude = (4 * total_mass * eta * np.sqrt(M_PI / 5.0)) / R

        GO_func = H_GO_LM_FUNCS.get((l, m))
        QC_func = H_QC_LM_FUNCS.get((l, m))
        
        # Sum the Generalized Orbital (GO) and Quasi-Circular (QC) terms
        waveform_modes = (
                        GO_func(total_mass, eta, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params) +
                        QC_func(total_mass, eta, vpnorder, x, S1z, S2z, params)
                        )
        if m < 0: waveform_modes = waveform_modes.conjugate()
        
        # cpolar(r, theta) in C is equivalent to rect(r, theta) in Python
        phase_factor = rect(1.0, -m * Phi)
        
        return amplitude * waveform_modes * phase_factor

def hlmGOresult(l: int, m: int, total_mass: float, eta: float, r: float, rDOT: float, Phi: float,
                PhiDOT: float, R: float, vpnorder: int, S1z: float, S2z: float, x: float
                ) -> complex:

    if not (0 <= vpnorder <= 8):
        raise ValueError("PN order must be between 0 and 8")

    if not (2 <= l <= 8):
        raise ValueError("l must be between 2 and 8")

    if not (-l <= m <= l):
        raise ValueError(f"m must be between {-l} and {l}")

    # Modes with m = 0 are zero in your C code
    if m == 0:
        return 0j
    
    b0 = 2 * total_mass / np.exp(0.5)
    r0 = b0
    params = CommonVars(
                xp5=np.sqrt(x),
                logx=np.log(x),
                b0 = b0,
                r0 = r0,
                logb0=np.log(b0),
                logr0=np.log(r0),
                delta=np.sqrt(1 - 4*eta),
            )
    hlm = 0j

    for pno in range(vpnorder, -1, -1):
        hlm += generate_hlm(l,m,total_mass, eta, r, rDOT, Phi, PhiDOT, R, pno, S1z, S2z, x, params)

    return hlm