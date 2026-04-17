import math
import cmath
from typing import Any

def Complex(real: float, imag: float) -> complex:
    """Helper function to mimic the C-style Complex constructor."""
    return complex(real, imag)

# All the hGO_l_m functions contain the 3PN non-spinning general orbit (GO) terms
# from Mishra et al. arXiv:1501.07096 and 3.5PN quasi-circular spinning
# corrections as given in Henry et al, arXiv:2209.00374v2 

# The (2,2) mode has newly computed 4PN non-spinning quasi-circular piece from
# arXiv:2304.11185 

# H22

def hGO_2_m_2(mass: float, Nu: float, r: float, rDOT: float,
              PhiDOT: float, vpnorder: int, S1z: float, S2z: float,
              x: float, params: Any) -> complex:
    
    # Note: In Python, `params` should be an object (like a dataclass or SimpleNamespace) 
    # so that attributes can be accessed via dot notation (e.g., params.PhiDOT2).
    
    # For black holes kappa and lambda is 1
    kappa1 = 1.0 
    kappa2 = 1.0
    lambda1 = 1.0
    lambda2 = 1.0
    
    delta = math.sqrt(1 - 4 * Nu)
    
    combination_a = (PhiDOT * r + Complex(0, 1) * rDOT)
    combination_a3 = combination_a * combination_a * combination_a
    combination_a4 = combination_a3 * combination_a
    combination_a5 = combination_a4 * combination_a

    combination_b = (PhiDOT * r - Complex(0, 1) * rDOT)
    combination_b2 = combination_b * combination_b
    combination_b3 = combination_b2 * combination_b

    # Mapping M_PI2 to pi^2 (standard in PN GW theory for this coefficient)
    M_PI = math.pi
    M_PI2 = math.pi ** 2

    if vpnorder == 0:
        return (mass / r + params.PhiDOT2 * params.r2 +
                Complex(0, 2) * PhiDOT * r * rDOT - params.rDOT2)

    elif vpnorder == 2:
        return ((21 * params.Mtot2 * (-10 + Nu) -
                 27 * (-1 + 3 * Nu) * params.r2 * combination_b * combination_a3 +
                 mass * r *
                     ((11 + 156 * Nu) * params.PhiDOT2 * params.r2 +
                      Complex(0, 10) * (5 + 27 * Nu) * PhiDOT * r * rDOT -
                      3 * (15 + 32 * Nu) * params.rDOT2)) /
                (42. * params.r2))

    elif vpnorder == 3:
        return ((params.Mtot2 * (Complex(0, -1) * rDOT *
                                 ((3 + 3 * delta - 8 * Nu) * S1z +
                                  (3 - 3 * delta - 8 * Nu) * S2z) +
                             PhiDOT * r *
                                 ((-3 - 3 * delta + 5 * Nu) * S1z +
                                  (-3 + 3 * delta + 5 * Nu) * S2z))) /
                (3. * params.r2)) 
        # (<--This is the general orbit term)
        # (This is the quasi-circular limit of the general orbit term-->)
        # - ((-4 * ((1 + delta - Nu) * S1z + S2z - (delta + Nu) * S2z) * params.x2p5) / 3.) +
        # ((-4 * (S1z + delta * S1z + S2z - delta * S2z - Nu * (S1z + S2z)) * params.x2p5) / 3.) 
        # (<--This is Quentins quasi-circular term)

    elif vpnorder == 4:
        return ((6 * params.Mtot3 * (3028 + 1267 * Nu + 158 * params.eta2) +
                 9 * (83 - 589 * Nu + 1111 * params.eta2) * params.r3 *
                     combination_b2 * combination_a4 +
                 params.Mtot2 * r *
                     ((-11891 - 36575 * Nu + 13133 * params.eta2) * params.PhiDOT2 *
                          params.r2 +
                      Complex(0, 8) * (-773 - 3767 * Nu + 2852 * params.eta2) *
                          PhiDOT * r * rDOT -
                      6 * (-619 + 2789 * Nu + 934 * params.eta2) * params.rDOT2) -
                 3 * mass * params.r2 *
                     (2 * (-835 - 19 * Nu + 2995 * params.eta2) * params.PhiDOT4 *
                          params.r4 +
                      Complex(0, 6) * (-433 - 721 * Nu + 1703 * params.eta2) *
                          params.PhiDOT3 * params.r3 * rDOT +
                      6 * (-33 + 1014 * Nu + 232 * params.eta2) * params.PhiDOT2 *
                          params.r2 * params.rDOT2 +
                      Complex(0, 4) * (-863 + 1462 * Nu + 2954 * params.eta2) *
                          PhiDOT * r * params.rDOT3 -
                      3 * (-557 + 664 * Nu + 1712 * params.eta2) * params.rDOT4)) /
                    (1512. * params.r3) +
                (3 * params.Mtot3 *
                 (S1z * (4 * Nu * S2z + (1 + delta - 2 * Nu) * S1z * kappa1) -
                  (-1 + delta + 2 * Nu) * params.S2z2 * kappa2)) /
                    (4. * params.r3))
                # This is where circular limit terms added and subtracted
                # - ((kappa1 * (1 + delta - 2 * Nu) * params.S1z2 + S2z * (4 * Nu * S1z - kappa2 * (-1 + delta + 2 * Nu) * S2z)) * params.x3) +
                # ((kappa1 * (1 + delta - 2 * Nu) * params.S1z2 + S2z * (4 * Nu * S1z - kappa2 * (-1 + delta + 2 * Nu) * S2z)) * params.x3)

    elif vpnorder == 5:
        return ((params.Mtot2 * Nu *
                 (2 * mass * (Complex(0, -702) * PhiDOT * r + rDOT) +
                  3 * r *
                      (Complex(0, -316) * params.PhiDOT3 * params.r3 -
                       847 * params.PhiDOT2 * params.r2 * rDOT +
                       Complex(0, 184) * PhiDOT * r * params.rDOT2 -
                       122 * params.rDOT3))) /
                (105. * params.r3) 
                # Henry et al. QC spin terms
                # + ((2*(56*delta*Nu*(-S1z + S2z) + 101*Nu*(S1z + S2z) + 132*params.eta2*(S1z + S2z) - 80*(S1z + delta*S1z + S2z - delta*S2z))*params.x3p5)/63.)
                +
                # Henry et al. ecc spin terms
                ((params.Mtot2 *
                 (mass *
                      ((238 + delta * (238 - 141 * Nu) + Nu * (-181 + 474 * Nu)) *
                           PhiDOT * r * S1z +
                       Complex(0, 8) *
                           (55 + delta * (55 - 19 * Nu) +
                            2 * Nu * (-50 + 43 * Nu)) *
                           rDOT * S1z +
                       (238 + delta * (-238 + 141 * Nu) + Nu * (-181 + 474 * Nu)) *
                           PhiDOT * r * S2z +
                       Complex(0, 8) *
                           (55 + delta * (-55 + 19 * Nu) +
                            2 * Nu * (-50 + 43 * Nu)) *
                           rDOT * S2z) -
                  r * (PhiDOT * r * params.rDOT2 *
                           (-((18 * (1 + delta) + 5 * (-63 + 55 * delta) * Nu +
                               188 * params.eta2) *
                              S1z) +
                            (18 * (-1 + delta) + 5 * (63 + 55 * delta) * Nu -
                             188 * params.eta2) *
                                S2z) -
                       Complex(0, 2) * params.rDOT3 *
                           ((-27 * (1 + delta) + 6 * (5 + 7 * delta) * Nu -
                             4 * params.eta2) *
                                S1z +
                            (-27 + 27 * delta + 30 * Nu - 42 * delta * Nu -
                             4 * params.eta2) *
                                S2z) +
                       Complex(0, 2) * params.PhiDOT2 * params.r2 * rDOT *
                           ((51 + 88 * Nu * (-3 + 5 * Nu) +
                             delta * (51 + 62 * Nu)) *
                                S1z +
                            (51 + 88 * Nu * (-3 + 5 * Nu) -
                             delta * (51 + 62 * Nu)) *
                                S2z) +
                       params.PhiDOT3 * params.r3 *
                           ((120 * (1 + delta) + (-483 + 83 * delta) * Nu +
                             234 * params.eta2) *
                                S1z +
                            (120 + 3 * Nu * (-161 + 78 * Nu) -
                             delta * (120 + 83 * Nu)) *
                                S2z)))) /
                (84. * params.r3)))

    elif vpnorder == 6:
        return ((4 * params.Mtot4 *
                 (-8203424 + 2180250 * params.eta2 + 592600 * params.eta3 +
                  15 * Nu * (-5503804 + 142065 * M_PI2)) -
             2700 * (-507 + 6101 * Nu - 25050 * params.eta2 + 34525 * params.eta3) *
                 params.r4 * combination_b3 * combination_a5 +
             params.Mtot3 * r *
                 (params.PhiDOT2 *
                      (337510808 - 198882000 * params.eta2 +
                       56294600 * params.eta3 +
                       Nu * (183074880 - 6392925 * M_PI2)) *
                      params.r2 +
                  Complex(0, 110) * PhiDOT *
                      (-5498800 - 785120 * params.eta2 + 909200 * params.eta3 +
                       3 * Nu * (-1849216 + 38745 * M_PI2)) *
                      r * rDOT +
                  2 *
                      (51172744 - 94929000 * params.eta2 - 5092400 * params.eta3 +
                       45 * Nu * (2794864 + 142065 * M_PI2)) *
                      params.rDOT2) -
             20 * params.Mtot2 * params.r2 *
                 ((-986439 + 1873255 * Nu - 9961400 * params.eta2 +
                   6704345 * params.eta3) *
                      params.PhiDOT4 * params.r4 +
                  Complex(0, 4) *
                      (-273687 - 978610 * Nu - 4599055 * params.eta2 +
                       2783005 * params.eta3) *
                      params.PhiDOT3 * params.r3 * rDOT +
                  (-181719 + 19395325 * Nu + 8237980 * params.eta2 +
                   2612735 * params.eta3) *
                      params.PhiDOT2 * params.r2 * params.rDOT2 +
                  Complex(0, 8) *
                      (-234312 + 1541140 * Nu + 1230325 * params.eta2 +
                       1828625 * params.eta3) *
                      PhiDOT * r * params.rDOT3 -
                  3 *
                      (-370268 + 1085140 * Nu + 2004715 * params.eta2 +
                       1810425 * params.eta3) *
                      params.rDOT4) +
             300 * mass * params.r3 *
                 (4 *
                      (12203 - 36427 * Nu - 27334 * params.eta2 +
                       149187 * params.eta3) *
                      params.PhiDOT6 * params.r6 +
                  Complex(0, 2) *
                      (44093 - 68279 * Nu - 295346 * params.eta2 +
                       541693 * params.eta3) *
                      params.PhiDOT5 * params.r5 * rDOT +
                  2 *
                      (27432 - 202474 * Nu + 247505 * params.eta2 +
                       394771 * params.eta3) *
                      params.PhiDOT4 * params.r4 * params.rDOT2 +
                  Complex(0, 2) *
                      (97069 - 383990 * Nu - 8741 * params.eta2 +
                       1264800 * params.eta3) *
                      params.PhiDOT3 * params.r3 * params.rDOT3 +
                  (-42811 + 53992 * Nu + 309136 * params.eta2 -
                   470840 * params.eta3) *
                      params.PhiDOT2 * params.r2 * params.rDOT4 +
                  Complex(0, 2) *
                      (51699 - 252256 * Nu + 131150 * params.eta2 +
                       681160 * params.eta3) *
                      PhiDOT * r * params.rDOT5 -
                  3 *
                      (16743 - 75104 * Nu + 26920 * params.eta2 +
                       207200 * params.eta3) *
                      params.rDOT6)) /
                (3.3264e6 * params.r4)
            # Henry et al. QC spin terms
            # + (((4*(1 + delta)*(-7 + 9*kappa1) - 7*(9 + 17*delta)*Nu - 9*(15 + 7*delta)*kappa1*Nu + 12*(7 - 17*kappa1)*params.eta2)*params.S1z2 + 2*S1z*(Complex(0,-42)*(1 + delta - 2*Nu) - 84*(1 + delta - Nu)*M_PI + Nu*(-271 + 288*Nu)*S2z) + S2z*(12*(7 - 17*kappa2)*params.eta2*S2z + 4*(-1 + delta)*(Complex(0,21) + 42*M_PI + 7*S2z - 9*kappa2*S2z) + Nu*(168*(Complex(0,1) + M_PI) + 7*delta*(17 + 9*kappa2)*S2z - 9*(7 + 15*kappa2)*S2z)))*params.x4)/63.
            +
            # Henry et al. ecc spin terms
            (-0.005952380952380952 *
                  (params.Mtot3 *
                   (2 * mass *
                        (S1z * (Complex(0, 14) * (1 + delta - 2 * Nu) +
                                42 * (1 + delta - 2 * Nu) * Nu * S1z +
                                kappa1 *
                                    (438 + delta * (438 + 103 * Nu) +
                                     Nu * (-773 + 108 * Nu)) *
                                    S1z) +
                         2 *
                             (Complex(0, -7) * (-1 + delta + 2 * Nu) +
                              (995 - 192 * Nu) * Nu * S1z) *
                             S2z -
                         (42 * Nu * (-1 + delta + 2 * Nu) +
                          kappa2 * (-438 + (773 - 108 * Nu) * Nu +
                                    delta * (438 + 103 * Nu))) *
                             params.S2z2) +
                    r * (params.rDOT2 *
                             (Complex(0, 56) * (1 + delta - 2 * Nu) * S1z +
                              (-56 * (1 + delta - 2 * Nu) * Nu +
                               kappa1 * (291 * (1 + delta) -
                                         2 * (445 + 154 * delta) * Nu +
                                         24 * params.eta2)) *
                                  params.S1z2 +
                              4 * Nu * (-3 + 44 * Nu) * S1z * S2z +
                              S2z * (Complex(0, -56) * (-1 + delta + 2 * Nu) +
                                     56 * Nu * (-1 + delta + 2 * Nu) * S2z +
                                     kappa2 *
                                         (291 - 890 * Nu + 24 * params.eta2 +
                                          delta * (-291 + 308 * Nu)) *
                                         S2z)) +
                         params.PhiDOT2 * params.r2 *
                             (Complex(0, 196) * (1 + delta - 2 * Nu) * S1z +
                              (56 * Nu * (7 + 7 * delta + Nu) +
                               kappa1 * (-153 * (1 + delta) -
                                         2 * (62 + 215 * delta) * Nu +
                                         804 * params.eta2)) *
                                  params.S1z2 +
                              8 * (60 - 187 * Nu) * Nu * S1z * S2z +
                              S2z * (Complex(0, -196) * (-1 + delta + 2 * Nu) +
                                     56 * Nu * (7 - 7 * delta + Nu) * S2z +
                                     kappa2 *
                                         (-153 + 4 * Nu * (-31 + 201 * Nu) +
                                          delta * (153 + 430 * Nu)) *
                                         S2z)) +
                         2 * PhiDOT * r * rDOT *
                             (Complex(0, -1) *
                                  (117 * (1 + delta) * kappa1 -
                                   434 * (1 + delta) * Nu +
                                   (-23 + 211 * delta) * kappa1 * Nu +
                                   2 * (14 - 195 * kappa1) * params.eta2) *
                                  params.S1z2 +
                              4 * S1z *
                                  (35 * (1 + delta - 2 * Nu) +
                                   Complex(0, 1) * (9 - 209 * Nu) * Nu * S2z) +
                              S2z * (-140 * (-1 + delta + 2 * Nu) +
                                     Complex(0, 1) *
                                         (-14 * Nu * (-31 + 31 * delta + 2 * Nu) +
                                          kappa2 * (117 * (-1 + delta) +
                                                    (23 + 211 * delta) * Nu +
                                                    390 * params.eta2)) *
                                         S2z))))) /
                  params.r4))
            # + Henry et al. QC spinning hereditary terms
            # (((-8 * M_PI * ((1 + delta - Nu) * S1z + S2z - (delta + Nu) * S2z) * params.x4) / 3.))

    elif vpnorder == 7:
        return (
            # Henry et al QC spin terms
            # ((3318*params.eta3*(S1z + S2z) + Nu*(-504*((7 + delta)*kappa1 - 3*(3 + delta)*lambda1)*params.S1z3 - 1008*params.S1z2*(3*kappa1*M_PI - 3*(1 + delta)*S2z + 2*(1 + delta)*kappa1*S2z) + S1z*(17387 + 20761*delta + 1008*S2z*(6*M_PI + (-1 + delta)*(-3 + 2*kappa2)*S2z)) + S2z*(17387 - 20761*delta + 504*S2z*(-6*kappa2*M_PI + (-7 + delta)*kappa2*S2z - 3*(-3 + delta)*lambda2*S2z))) + 2*(2809*(1 + delta)*S1z + 756*(1 + delta)*kappa1*M_PI*params.S1z2 + 756*(1 + delta)*(kappa1 - lambda1)*params.S1z3 - (-1 + delta)*S2z*(2809 + 756*S2z*(-(lambda2*S2z) + kappa2*(M_PI + S2z)))) - 2*params.eta2*(708*delta*(-S1z + S2z) + (S1z + S2z)*(4427 + 1008*(kappa1*params.S1z2 + S2z*(-2*S1z + kappa2*S2z)))))*params.x4p5)/756.
            # +
            # Henry et al. ecc+spin terms
            ((params.Mtot2 *
                 (-3 * mass * r *
                      (Complex(0, -16) * params.rDOT3 *
                           (Complex(0, 12) * Nu * (-16703 + 4427 * Nu) +
                            35 * Nu *
                                (4578 + Nu * (-4288 + 765 * Nu) +
                                 delta * (3748 + 802 * Nu)) *
                                S1z +
                            35 *
                                (delta * (942 - 2 * Nu * (1874 + 401 * Nu)) +
                                 Nu * (4578 + Nu * (-4288 + 765 * Nu))) *
                                S2z -
                            32970 * (S1z + delta * S1z + S2z)) +
                       4 * PhiDOT * r * params.rDOT2 *
                           (-338520 * params.eta3 * (S1z + S2z) +
                            48930 * (S1z + delta * S1z + S2z - delta * S2z) +
                            Nu * (Complex(0, 3420696) -
                                  35 * (14154 + 21167 * delta) * S1z -
                                  495390 * S2z + 740845 * delta * S2z) +
                            params.eta2 * (Complex(0, -612336) +
                                           245 * (3566 - 1565 * delta) * S1z +
                                           245 * (3566 + 1565 * delta) * S2z)) +
                       params.PhiDOT3 * params.r3 *
                           (2515380 * params.eta3 * (S1z + S2z) -
                            5 * Nu *
                                (Complex(0, 1859936) +
                                 7 * (82329 + 37061 * delta) * S1z +
                                 7 * (82329 - 37061 * delta) * S2z) -
                            128100 * (S1z + delta * S1z + S2z - delta * S2z) +
                            4 * params.eta2 *
                                (Complex(0, 381348) +
                                 35 * (-18505 + 1777 * delta) * S1z -
                                 35 * (18505 + 1777 * delta) * S2z)) +
                       Complex(0, 8) * params.PhiDOT2 * params.r2 * rDOT *
                           (779100 * params.eta3 * (S1z + S2z) +
                            5 * Nu *
                                (Complex(0, 828806) +
                                 7 * (4839 + 5971 * delta) * S1z +
                                 7 * (4839 - 5971 * delta) * S2z) -
                            62475 * (S1z + delta * S1z + S2z - delta * S2z) +
                            params.eta2 * (Complex(0, -976002) +
                                           35 * (-29599 + 3109 * delta) * S1z -
                                           35 * (29599 + 3109 * delta) * S2z))) +
                  3 * params.r2 *
                      (Complex(0, 4) * params.PhiDOT2 * params.r2 * params.rDOT3 *
                           (Complex(0, 2) * (65451 - 350563 * Nu) * Nu +
                            105 * Nu *
                                (6020 + delta * (3114 + 411 * Nu) +
                                 Nu * (-4513 + 9136 * Nu)) *
                                S1z +
                            105 *
                                (-3 * delta * (-408 + Nu * (1038 + 137 * Nu)) +
                                 Nu * (6020 + Nu * (-4513 + 9136 * Nu))) *
                                S2z -
                            128520 * (S1z + delta * S1z + S2z)) +
                       PhiDOT * r * params.rDOT4 *
                           (Complex(0, -128) * Nu * (2487 + 18334 * Nu) -
                            105 * Nu *
                                (8689 + 8 * Nu * (-687 + 1402 * Nu) +
                                 delta * (2143 + 8212 * Nu)) *
                                S1z +
                            105 *
                                (Nu * (-8689 + 8 * (687 - 1402 * Nu) * Nu) +
                                 delta * (-1470 + Nu * (2143 + 8212 * Nu))) *
                                S2z +
                            154350 * (S1z + delta * S1z + S2z)) -
                       Complex(0, 6) * params.rDOT5 *
                           (-11760 * params.eta3 * (S1z + S2z) +
                            4 * params.eta2 *
                                (Complex(0, 57338) +
                                 35 * (569 + 301 * delta) * S1z +
                                 35 * (569 - 301 * delta) * S2z) -
                            16 * Nu *
                                (Complex(0, -2517) + 35 * (72 + 71 * delta) * S1z +
                                 35 * (72 - 71 * delta) * S2z) +
                            8715 * (S1z + delta * S1z + S2z - delta * S2z)) +
                       params.PhiDOT5 * params.r5 *
                           (2263380 * params.eta3 * (S1z + S2z) +
                            3 * Nu *
                                (Complex(0, 653432) +
                                 35 * (13283 + 7839 * delta) * S1z +
                                 35 * (13283 - 7839 * delta) * S2z) -
                            219240 * (S1z + delta * S1z + S2z - delta * S2z) +
                            4 * params.eta2 *
                                (Complex(0, -291268) +
                                 105 * (-7669 + 165 * delta) * S1z -
                                 105 * (7669 + 165 * delta) * S2z)) +
                       Complex(0, 14) * params.PhiDOT4 * params.r4 * rDOT *
                           (385170 * params.eta3 * (S1z + S2z) +
                            15 * Nu *
                                (Complex(0, 16914) + (8633 + 2267 * delta) * S1z +
                                 8633 * S2z - 2267 * delta * S2z) -
                            18630 * (S1z + delta * S1z + S2z - delta * S2z) +
                            2 * params.eta2 *
                                (Complex(0, -67904) +
                                 15 * (-13932 + 679 * delta) * S1z -
                                 15 * (13932 + 679 * delta) * S2z)) +
                       6 * params.PhiDOT3 * params.r3 * params.rDOT2 *
                           (-6720 * params.eta3 * (S1z + S2z) +
                            5 * Nu *
                                (Complex(0, -241704) +
                                 7 * (4377 + 3083 * delta) * S1z +
                                 7 * (4377 - 3083 * delta) * S2z) -
                            36960 * (S1z + delta * S1z + S2z - delta * S2z) +
                            params.eta2 * (Complex(0, -364384) +
                                           35 * (5059 - 4753 * delta) * S1z +
                                           35 * (5059 + 4753 * delta) * S2z))) +
                  14 * params.Mtot2 *
                      (Complex(0, 30) * rDOT *
                           (15280 * params.eta3 * (S1z + S2z) -
                            4 * params.eta2 *
                                (Complex(0, -111) + 3562 * S1z + 682 * delta * S1z +
                                 945 * kappa1 * params.S1z3 +
                                 (3562 - 682 * delta +
                                  945 * (-2 + kappa1) * params.S1z2) *
                                     S2z +
                                 945 * (-2 + kappa2) * S1z * params.S2z2 +
                                 945 * kappa2 * params.S2z3) +
                            Nu * (Complex(0, -27520) +
                                  378 *
                                      (5 * (1 + delta) * kappa1 +
                                       3 * (3 + delta) * lambda1) *
                                      params.S1z3 +
                                  (29749 - 13605 * delta) * S2z -
                                  1512 * (1 + delta) * kappa1 * params.S1z2 * S2z -
                                  378 *
                                      (5 * (-1 + delta) * kappa2 +
                                       3 * (-3 + delta) * lambda2) *
                                      params.S2z3 +
                                  S1z * (29749 + 13605 * delta +
                                         1512 * (-1 + delta) * kappa2 *
                                             params.S2z2)) +
                            2 * (-8009 * (1 + delta) * S1z -
                                 567 * (1 + delta) * lambda1 * params.S1z3 +
                                 (-1 + delta) * S2z *
                                     (8009 + 567 * lambda2 * params.S2z2))) +
                       PhiDOT * r *
                           (285840 * params.eta3 * (S1z + S2z) -
                            30 * params.eta2 *
                                (Complex(0, 11504) + 1890 * kappa1 * params.S1z3 +
                                 23090 * S2z - 1823 * delta * S2z +
                                 1890 * (-2 + kappa1) * params.S1z2 * S2z +
                                 1890 * kappa2 * params.S2z3 +
                                 S1z * (23090 + 1823 * delta +
                                        1890 * (-2 + kappa2) * params.S2z2)) +
                            30 * (689 * (1 + delta) * S1z -
                                  1134 * (1 + delta) * lambda1 * params.S1z3 +
                                  (-1 + delta) * S2z *
                                      (-689 + 1134 * lambda2 * params.S2z2)) +
                            2 * Nu *
                                (Complex(0, 415432) - 66840 * S2z +
                                 15 * (8 * (-557 + 1342 * delta) * S1z +
                                       189 *
                                           (5 * (1 + delta) * kappa1 +
                                            6 * (3 + delta) * lambda1) *
                                           params.S1z3 -
                                       10736 * delta * S2z -
                                       2457 * (1 + delta) * kappa1 * params.S1z2 *
                                           S2z +
                                       2457 * (-1 + delta) * kappa2 * S1z *
                                           params.S2z2 -
                                       189 *
                                           (5 * (-1 + delta) * kappa2 +
                                            6 * (-3 + delta) * lambda2) *
                                           params.S2z3)))))) /
                (317520. * params.r4))
            # + Henry et al. QC spinning hereditary terms
            # (2 * M_PI * (kappa1 * (1 + delta - 2 * Nu) * params.S1z2 + S2z * (4 * Nu * S1z - kappa2 * (-1 + delta + 2 * Nu) * S2z)) * params.x4p5)
        )

    else:
        return complex(0, 0)
    


# hQC_l_m() functions contains only the non-spinning hereditary terms at
# particular PN order.

def hQC_2_m_2(mass: float, Nu: float, vpnorder: int, x: float, S1z: float, S2z: float,
              params: Any) -> complex:
    
    EulerGamma = 0.5772156649015329
    x0 = 0.4375079683656479
    b0 = 2 * mass / math.exp(0.5)
    r0 = b0
    kappa1 = 1.0 # for black holes kappa and lambda is 1
    kappa2 = 1.0
    delta = math.sqrt(1 - 4 * Nu)
    
    M_PI = math.pi
    M_PI2 = math.pi ** 2

    # keeping only the hereditary terms
    # if vpnorder == 0:
    #     return (2 * x)

    # elif vpnorder == 2:
    #     return ((-5.095238095238095 + (55 * Nu) / 21.) * params.x2)

    if vpnorder == 3:
        # return (4 * M_PI * params.x2p5)
        return (complex(0, 0.6666666666666666) * params.x2p5 * (-11 + 12 * EulerGamma - complex(0, 6) * M_PI + 24 * math.log(2) + 12 * math.log(b0) + 18 * math.log(x)))

    # elif vpnorder == 4:
    #     return ((-2.874338624338624 - (1069 * Nu) / 108. + (2047 * params.eta2) / 756.) * params.x3)

    elif vpnorder == 5:
        # return ((complex(0,-48)*Nu - (214*M_PI)/21. + (68*Nu*M_PI)/21.)*params.x3p5)
        # return ((2 * (-107 + 34 * Nu) * M_PI * params.x3p5) / 21.) # This is old implementation
        return (complex(0, 0.015873015873015872) * (-107 + 34 * Nu) * params.x3p5 * (-11 + 12 * EulerGamma - complex(0, 6) * M_PI + 24 * math.log(2) + 12 * math.log(b0) + 18 * math.log(x)))

    elif vpnorder == 6:
        # return ((params.x4 * (-27392 * EulerGamma + M_PI * (complex(0, 13696) + 35 * (64 + 41 * Nu) * M_PI) - 13696 * math.log(16 * x))) / 1680.) # This is the old implementation.

        return (params.x4 * (-45.216145124716554 + (456 * EulerGamma) / 35. - 16 * EulerGamma * EulerGamma - 
                complex(0, 6.514285714285714) * M_PI + complex(0, 16) * EulerGamma * M_PI + (4 * M_PI2) / 3. + 
                complex(0, 4.888888888888889) * S1z - complex(0, 5.333333333333333) * EulerGamma * S1z - 
                complex(0, 4.888888888888889) * Nu * S1z + complex(0, 5.333333333333333) * EulerGamma * Nu * S1z - 
                (8 * M_PI * S1z) / 3. + (8 * Nu * M_PI * S1z) / 3. + complex(0, 4.888888888888889) * S2z - 
                complex(0, 5.333333333333333) * EulerGamma * S2z - complex(0, 4.888888888888889) * Nu * S2z + 
                complex(0, 5.333333333333333) * EulerGamma * Nu * S2z - (8 * M_PI * S2z) / 3. + (8 * Nu * M_PI * S2z) / 3. + 
                (912 * math.log(2)) / 35. - 64 * EulerGamma * math.log(2) + complex(0, 32) * M_PI * math.log(2) - 
                complex(0, 10.666666666666666) * S1z * math.log(2) + complex(0, 10.666666666666666) * Nu * S1z * math.log(2) - 
                complex(0, 10.666666666666666) * S2z * math.log(2) + complex(0, 10.666666666666666) * Nu * S2z * math.log(2) - 
                64 * math.log(2) * math.log(2) + (88 * math.log(b0)) / 3. - 32 * EulerGamma * math.log(b0) + 
                complex(0, 16) * M_PI * math.log(b0) - complex(0, 5.333333333333333) * S1z * math.log(b0) + 
                complex(0, 5.333333333333333) * Nu * S1z * math.log(b0) - complex(0, 5.333333333333333) * S2z * math.log(b0) + 
                complex(0, 5.333333333333333) * Nu * S2z * math.log(b0) - 64 * math.log(2) * math.log(b0) - 
                16 * math.log(b0) * math.log(b0) - (1712 * math.log(r0)) / 105. + (684 * math.log(x)) / 35. - 
                48 * EulerGamma * math.log(x) + complex(0, 24) * M_PI * math.log(x) - complex(0, 8) * S1z * math.log(x) + 
                complex(0, 8) * Nu * S1z * math.log(x) - complex(0, 8) * S2z * math.log(x) + complex(0, 8) * Nu * S2z * math.log(x) - 
                96 * math.log(2) * math.log(x) - 48 * math.log(b0) * math.log(x) - 36 * math.log(x) * math.log(x) - 
                complex(0, 0.4444444444444444) * delta * (S1z - S2z) * (-11 + 12 * EulerGamma - complex(0, 6) * M_PI + 
                24 * math.log(2) + 12 * math.log(b0) + 18 * math.log(x))))

    elif vpnorder == 7:
        # return (((-2173 - 4990 * Nu + 1120 * params.eta2) * M_PI * params.x4p5) / 378.) # This is the old implementation.
        return (complex(0, 0.0004409171075837742) * params.x4p5 * (23903 - 26076 * EulerGamma + 57452 * Nu - 
                59880 * EulerGamma * Nu - 18728 * params.eta2 + 13440 * EulerGamma * params.eta2 + 
                complex(0, 13038) * M_PI + complex(0, 29940) * Nu * M_PI - complex(0, 6720) * params.eta2 * M_PI - 
                8316 * kappa1 * params.S1z2 - 8316 * delta * kappa1 * params.S1z2 + 9072 * EulerGamma * kappa1 * params.S1z2 + 
                9072 * delta * EulerGamma * kappa1 * params.S1z2 + 16632 * kappa1 * Nu * params.S1z2 - 
                18144 * EulerGamma * kappa1 * Nu * params.S1z2 - complex(0, 4536) * kappa1 * M_PI * params.S1z2 - 
                complex(0, 4536) * delta * kappa1 * M_PI * params.S1z2 + complex(0, 9072) * kappa1 * Nu * M_PI * params.S1z2 - 
                33264 * Nu * S1z * S2z + 36288 * EulerGamma * Nu * S1z * S2z - complex(0, 18144) * Nu * M_PI * S1z * S2z - 
                8316 * kappa2 * params.S2z2 + 8316 * delta * kappa2 * params.S2z2 + 9072 * EulerGamma * kappa2 * params.S2z2 - 
                9072 * delta * EulerGamma * kappa2 * params.S2z2 + 16632 * kappa2 * Nu * params.S2z2 - 
                18144 * EulerGamma * kappa2 * Nu * params.S2z2 - complex(0, 4536) * kappa2 * M_PI * params.S2z2 + 
                complex(0, 4536) * delta * kappa2 * M_PI * params.S2z2 + complex(0, 9072) * kappa2 * Nu * M_PI * params.S2z2 - 
                52152 * math.log(2) - 119760 * Nu * math.log(2) + 26880 * params.eta2 * math.log(2) + 
                18144 * kappa1 * params.S1z2 * math.log(2) + 18144 * delta * kappa1 * params.S1z2 * math.log(2) - 
                36288 * kappa1 * Nu * params.S1z2 * math.log(2) + 72576 * Nu * S1z * S2z * math.log(2) + 
                18144 * kappa2 * params.S2z2 * math.log(2) - 18144 * delta * kappa2 * params.S2z2 * math.log(2) - 
                36288 * kappa2 * Nu * params.S2z2 * math.log(2) + 12 * (-2173 - 4990 * Nu + 1120 * params.eta2 + 
                756 * kappa1 * (1 + delta - 2 * Nu) * params.S1z2 + 3024 * Nu * S1z * S2z - 
                756 * kappa2 * (-1 + delta + 2 * Nu) * params.S2z2) * math.log(b0) + 18 * (-2173 - 4990 * Nu + 
                1120 * params.eta2 + 756 * kappa1 * (1 + delta - 2 * Nu) * params.S1z2 + 3024 * Nu * S1z * S2z - 
                756 * kappa2 * (-1 + delta + 2 * Nu) * params.S2z2) * math.log(x)))

    # 4PN non-spinning quasi-circular (2,2) mode has been obtained from Blanchet
    # et al. arXiv:2304.11185. This is written in terms of phi. Please refer to file shared by Quentin.

    elif vpnorder == 8:
        return (-1.3109423540262542e-11 * (params.x5 * (29059430400 * EulerGamma * EulerGamma * (-107 + 13 * Nu) + 
                12 * (628830397253 + 1854914893791 * Nu + 421984442880 * math.log(2)) + 415134720 * EulerGamma * (6099 - 40277 * Nu + complex(0, 70) * (107 - 13 * Nu) * M_PI + 280 * (-107 + 13 * Nu) * math.log(2)) + 
                35 * (-28 * params.eta2 * (5385456111 + 5 * Nu * (-163158374 + 26251249 * Nu)) + 
                135135 * (54784 + 5 * Nu * (1951 + 6560 * Nu)) * M_PI2 - 955450349568 * Nu * math.log(2) + 
                3321077760 * (-107 + 13 * Nu) * math.log(2) * math.log(2) - complex(0, 5930496) * M_PI * (6099 - 74773 * Nu + 280 * (-107 + 13 * Nu) * math.log(2))) - 5700491596800 * math.log(mass) + 
                6966444925440 * math.log(x) + 69189120 * (420 * (-107 + 13 * Nu) * math.log(b0) * math.log(b0) + 
                420 * (-107 + 13 * Nu) * math.log(mass) * math.log(mass) - 14 * math.log(mass) * (-11803 * Nu + 
                60 * EulerGamma * (-107 + 13 * Nu) + complex(0, 30) * (107 - 13 * Nu) * M_PI + 
                120 * (-107 + 13 * Nu) * math.log(2) + 90 * (-107 + 13 * Nu) * math.log(x)) + 14 * math.log(b0) * (5885 - 6420 * EulerGamma - 11803 * Nu + 780 * EulerGamma * Nu + complex(0, 3210) * M_PI - 
                complex(0, 390) * Nu * M_PI + 120 * (-107 + 13 * Nu) * math.log(2) + (6420 - 780 * Nu) * math.log(mass) + 
                90 * (-107 + 13 * Nu) * math.log(x)) + 5 * math.log(x) * (-74149 * Nu + 252 * EulerGamma * (-107 + 13 * Nu) - 
                complex(0, 126) * (-107 + 13 * Nu) * M_PI + 504 * (-107 + 13 * Nu) * math.log(2) + 
                189 * (-107 + 13 * Nu) * math.log(x)) + 84672 * Nu * math.log(x0)))))
        
        # return ((params.x5 * (276756480 * EulerGamma * (11449 + 19105 * Nu) - 12 * (846557506853 + 1008017482431 * Nu) + 35 * (28 * params.eta2 * (5385456111 + 5 * Nu * (-163158374 + 26251249 * Nu)) - complex(0, 3953664) * (11449 + 109657 * Nu) * M_PI - 135135 * (54784 + 5 * Nu * (1951 + 6560 * Nu)) * M_PI2) + 138378240 * (11449 + 19105 * Nu) * math.log(16 * x))) / 7.62810048e10)

    else:
        return complex(0, 0)
    

def hl_2_m_2(mass: float, Nu: float, r: float, rDOT: float, Phi: float,
             PhiDOT: float, R: float, vpnorder: int, S1z: float,
             S2z: float, x: float, params: Any) -> complex:

    if vpnorder < 0 or vpnorder > 8:
        raise ValueError("Error in hl_2_m_2: Input PN order parameter should be between [0, 8].")

    else:
        # Calculate the leading amplitude coefficient
        amplitude = (4 * mass * Nu * math.sqrt(math.pi / 5.0)) / R
        
        # Sum the Generalized Orbital (GO) and Quasi-Circular (QC) terms
        waveform_modes = (
            hGO_2_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params) +
            hQC_2_m_2(mass, Nu, vpnorder, x, S1z, S2z, params)
        )
        
        # cpolar(r, theta) in C is equivalent to cmath.rect(r, theta) in Python
        phase_factor = cmath.rect(1.0, -2 * Phi)
        
        return amplitude * waveform_modes * phase_factor
    
def hl_2_m_min2(mass: float, Nu: float, r: float, rDOT: float, Phi: float,
                PhiDOT: float, R: float, vpnorder: int, S1z: float,
                S2z: float, x: float, params: Any) -> complex:

    if vpnorder < 0 or vpnorder > 8:
        raise ValueError("Error in hl_2_m_min2: Input PN order parameter should be between [0, 8].")

    else:
        # Calculate the leading amplitude coefficient
        amplitude = (4 * mass * Nu * math.sqrt(math.pi / 5.0)) / R
        
        # Sum the GO and QC terms, then apply the complex conjugate for the m = -2 mode
        waveform_modes = (
            hGO_2_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params) +
            hQC_2_m_2(mass, Nu, vpnorder, x, S1z, S2z, params)
        ).conjugate()
        
        # Calculate the phase factor with positive 2 * Phi
        phase_factor = cmath.rect(1.0, 2 * Phi)
        
        return amplitude * waveform_modes * phase_factor