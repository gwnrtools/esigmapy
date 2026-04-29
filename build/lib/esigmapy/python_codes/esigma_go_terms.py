"""
esigma_go_terms.py
====================
Python translation of ESIGMA_GOTerms.c
"""

import numpy as np
from dataclasses import dataclass

# Mapping M_PI2 to pi^2 (standard in PN GW theory for this coefficient)
M_PI = np.pi
M_PI2 = M_PI ** 2

def Complex(real: float, imag: float) -> complex:
    """Helper function to mimic the C-style Complex constructor."""
    return complex(real, imag)
def rect(r: float, phi: float) -> complex:
    """Convert polar coordinates to rectangular form."""
    return r * np.exp(1j * phi)

@dataclass
class KeplerVars:
    # Basic orbital elements
    pn_order: int
    eta: float
    Mtot: float
    x: float
    e: float
    l: float
    r: float
    rDOT: float
    PhiDOT: float
    S1z: float
    S2z: float

    # Powers of x
    x2: float = 0.0;   x2p5: float = 0.0; x3: float = 0.0;   x3p5: float = 0.0
    x4: float = 0.0;   x4p5: float = 0.0; x5: float = 0.0;   x6: float = 0.0
    x7: float = 0.0;   x8: float = 0.0;   x9: float = 0.0

    # Powers of e
    e2: float = 0.0;   e3: float = 0.0;   e4: float = 0.0;   e5: float = 0.0
    e6: float = 0.0;   e7: float = 0.0;   e8: float = 0.0;   e9: float = 0.0

    # Powers of r
    r2: float = 0.0;   r3: float = 0.0;   r4: float = 0.0;   r5: float = 0.0
    r6: float = 0.0;   r7: float = 0.0;   r8: float = 0.0;   r9: float = 0.0

    # Powers of rDOT
    rDOT2: float = 0.0; rDOT3: float = 0.0; rDOT4: float = 0.0; rDOT5: float = 0.0
    rDOT6: float = 0.0; rDOT7: float = 0.0; rDOT8: float = 0.0; rDOT9: float = 0.0

    # Powers of PhiDOT
    PhiDOT2: float = 0.0; PhiDOT3: float = 0.0; PhiDOT4: float = 0.0; PhiDOT5: float = 0.0
    PhiDOT6: float = 0.0; PhiDOT7: float = 0.0; PhiDOT8: float = 0.0; PhiDOT9: float = 0.0

    # Powers of S1z
    S1z2: float = 0.0;  S1z3: float = 0.0;  S1z4: float = 0.0;  S1z5: float = 0.0
    S1z6: float = 0.0;  S1z7: float = 0.0;  S1z8: float = 0.0;  S1z9: float = 0.0
    S1z10: float = 0.0; S1z11: float = 0.0; S1z12: float = 0.0; S1z13: float = 0.0
    S1z14: float = 0.0; S1z15: float = 0.0; S1z16: float = 0.0

    # Powers of S2z
    S2z2: float = 0.0;  S2z3: float = 0.0;  S2z4: float = 0.0
    S2z5: float = 0.0;  S2z6: float = 0.0

    # Powers of eta
    eta2: float = 0.0;  eta3: float = 0.0;  eta4: float = 0.0
    eta5: float = 0.0;  eta6: float = 0.0

    # Powers of Mtot
    Mtot2: float = 0.0; Mtot3: float = 0.0; Mtot4: float = 0.0
    Mtot5: float = 0.0; Mtot6: float = 0.0

    def compute_powers(self):
        """Compute all power fields from the base values."""
        self.x2 = self.x**2;     self.x2p5 = self.x**2.5; self.x3 = self.x**3
        self.x3p5 = self.x**3.5; self.x4 = self.x**4;     self.x4p5 = self.x**4.5
        self.x5 = self.x**5;     self.x6 = self.x**6;     self.x7 = self.x**7
        self.x8 = self.x**8;     self.x9 = self.x**9

        self.e2 = self.e**2; self.e3 = self.e**3; self.e4 = self.e**4; self.e5 = self.e**5
        self.e6 = self.e**6; self.e7 = self.e**7; self.e8 = self.e**8; self.e9 = self.e**9

        self.r2 = self.r**2; self.r3 = self.r**3; self.r4 = self.r**4; self.r5 = self.r**5
        self.r6 = self.r**6; self.r7 = self.r**7; self.r8 = self.r**8; self.r9 = self.r**9

        self.rDOT2 = self.rDOT**2; self.rDOT3 = self.rDOT**3; self.rDOT4 = self.rDOT**4
        self.rDOT5 = self.rDOT**5; self.rDOT6 = self.rDOT**6; self.rDOT7 = self.rDOT**7
        self.rDOT8 = self.rDOT**8; self.rDOT9 = self.rDOT**9

        self.PhiDOT2 = self.PhiDOT**2; self.PhiDOT3 = self.PhiDOT**3; self.PhiDOT4 = self.PhiDOT**4
        self.PhiDOT5 = self.PhiDOT**5; self.PhiDOT6 = self.PhiDOT**6; self.PhiDOT7 = self.PhiDOT**7
        self.PhiDOT8 = self.PhiDOT**8; self.PhiDOT9 = self.PhiDOT**9

        self.S1z2 = self.S1z**2;  self.S1z3 = self.S1z**3;  self.S1z4 = self.S1z**4
        self.S1z5 = self.S1z**5;  self.S1z6 = self.S1z**6;  self.S1z7 = self.S1z**7
        self.S1z8 = self.S1z**8;  self.S1z9 = self.S1z**9;  self.S1z10 = self.S1z**10
        self.S1z11 = self.S1z**11; self.S1z12 = self.S1z**12; self.S1z13 = self.S1z**13
        self.S1z14 = self.S1z**14; self.S1z15 = self.S1z**15; self.S1z16 = self.S1z**16

        self.S2z2 = self.S2z**2; self.S2z3 = self.S2z**3; self.S2z4 = self.S2z**4
        self.S2z5 = self.S2z**5; self.S2z6 = self.S2z**6

        self.eta2 = self.eta**2; self.eta3 = self.eta**3; self.eta4 = self.eta**4
        self.eta5 = self.eta**5; self.eta6 = self.eta**6

        self.Mtot2 = self.Mtot**2; self.Mtot3 = self.Mtot**3; self.Mtot4 = self.Mtot**4
        self.Mtot5 = self.Mtot**5; self.Mtot6 = self.Mtot**6

# All the hGO_l_m functions contain the 3PN non-spinning general orbit (GO) terms
# from Mishra et al. arXiv:1501.07096 and 3.5PN quasi-circular spinning
# corrections as given in Henry et al, arXiv:2209.00374v2 

# The (2,2) mode has newly computed 4PN non-spinning quasi-circular piece from
# arXiv:2304.11185 

# H22

def hGO_2_m_2(mass: float, Nu: float, r: float, rDOT: float,
              PhiDOT: float, vpnorder: int, S1z: float, S2z: float,
              x: float, params: KeplerVars) -> complex:
    
    # Note: In Python, `params` should be an object (like a dataclass or SimpleNamespace) 
    # so that attributes can be accessed via dot notation (e.g., params.PhiDOT2).
    
    # For black holes kappa and lambda is 1
    kappa1 = 1.0 
    kappa2 = 1.0
    lambda1 = 1.0
    lambda2 = 1.0
    
    delta = np.sqrt(1 - 4 * Nu)
    
    combination_a = (PhiDOT * r + Complex(0, 1) * rDOT)
    combination_a3 = combination_a * combination_a * combination_a
    combination_a4 = combination_a3 * combination_a
    combination_a5 = combination_a4 * combination_a

    combination_b = (PhiDOT * r - Complex(0, 1) * rDOT)
    combination_b2 = combination_b * combination_b
    combination_b3 = combination_b2 * combination_b

    

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
    
# hQC_l_m() functions contain only the non-spinning hereditary terms at
# particular PN order.

def hQC_2_m_2(mass: float, Nu: float, vpnorder: int, x: float, S1z: float, S2z: float,
              params: KeplerVars) -> complex:
    
    EulerGamma = 0.5772156649015329
    x0 = 0.4375079683656479
    b0 = 2 * mass / np.exp(0.5)
    r0 = b0
    kappa1 = 1.0 # for black holes kappa and lambda is 1
    kappa2 = 1.0
    delta = np.sqrt(1 - 4 * Nu)

    # keeping only the hereditary terms
    # if vpnorder == 0:
    #     return (2 * x)

    # elif vpnorder == 2:
    #     return ((-5.095238095238095 + (55 * Nu) / 21.) * params.x2)

    if vpnorder == 3:
        # return (4 * M_PI * params.x2p5)
        return (complex(0, 0.6666666666666666) * params.x2p5 * (-11 + 12 * EulerGamma - complex(0, 6) * M_PI + 24 * np.log(2) + 12 * np.log(b0) + 18 * np.log(x)))

    # elif vpnorder == 4:
    #     return ((-2.874338624338624 - (1069 * Nu) / 108. + (2047 * params.eta2) / 756.) * params.x3)

    elif vpnorder == 5:
        # return ((complex(0,-48)*Nu - (214*M_PI)/21. + (68*Nu*M_PI)/21.)*params.x3p5)
        # return ((2 * (-107 + 34 * Nu) * M_PI * params.x3p5) / 21.) # This is old implementation
        return (complex(0, 0.015873015873015872) * (-107 + 34 * Nu) * params.x3p5 * (-11 + 12 * EulerGamma - complex(0, 6) * M_PI + 24 * np.log(2) + 12 * np.log(b0) + 18 * np.log(x)))

    elif vpnorder == 6:
        # return ((params.x4 * (-27392 * EulerGamma + M_PI * (complex(0, 13696) + 35 * (64 + 41 * Nu) * M_PI) - 13696 * np.log(16 * x))) / 1680.) # This is the old implementation.

        return (params.x4 * (-45.216145124716554 + (456 * EulerGamma) / 35. - 16 * EulerGamma * EulerGamma - 
                complex(0, 6.514285714285714) * M_PI + complex(0, 16) * EulerGamma * M_PI + (4 * M_PI2) / 3. + 
                complex(0, 4.888888888888889) * S1z - complex(0, 5.333333333333333) * EulerGamma * S1z - 
                complex(0, 4.888888888888889) * Nu * S1z + complex(0, 5.333333333333333) * EulerGamma * Nu * S1z - 
                (8 * M_PI * S1z) / 3. + (8 * Nu * M_PI * S1z) / 3. + complex(0, 4.888888888888889) * S2z - 
                complex(0, 5.333333333333333) * EulerGamma * S2z - complex(0, 4.888888888888889) * Nu * S2z + 
                complex(0, 5.333333333333333) * EulerGamma * Nu * S2z - (8 * M_PI * S2z) / 3. + (8 * Nu * M_PI * S2z) / 3. + 
                (912 * np.log(2)) / 35. - 64 * EulerGamma * np.log(2) + complex(0, 32) * M_PI * np.log(2) - 
                complex(0, 10.666666666666666) * S1z * np.log(2) + complex(0, 10.666666666666666) * Nu * S1z * np.log(2) - 
                complex(0, 10.666666666666666) * S2z * np.log(2) + complex(0, 10.666666666666666) * Nu * S2z * np.log(2) - 
                64 * np.log(2) * np.log(2) + (88 * np.log(b0)) / 3. - 32 * EulerGamma * np.log(b0) + 
                complex(0, 16) * M_PI * np.log(b0) - complex(0, 5.333333333333333) * S1z * np.log(b0) + 
                complex(0, 5.333333333333333) * Nu * S1z * np.log(b0) - complex(0, 5.333333333333333) * S2z * np.log(b0) + 
                complex(0, 5.333333333333333) * Nu * S2z * np.log(b0) - 64 * np.log(2) * np.log(b0) - 
                16 * np.log(b0) * np.log(b0) - (1712 * np.log(r0)) / 105. + (684 * np.log(x)) / 35. - 
                48 * EulerGamma * np.log(x) + complex(0, 24) * M_PI * np.log(x) - complex(0, 8) * S1z * np.log(x) + 
                complex(0, 8) * Nu * S1z * np.log(x) - complex(0, 8) * S2z * np.log(x) + complex(0, 8) * Nu * S2z * np.log(x) - 
                96 * np.log(2) * np.log(x) - 48 * np.log(b0) * np.log(x) - 36 * np.log(x) * np.log(x) - 
                complex(0, 0.4444444444444444) * delta * (S1z - S2z) * (-11 + 12 * EulerGamma - complex(0, 6) * M_PI + 
                24 * np.log(2) + 12 * np.log(b0) + 18 * np.log(x))))

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
                52152 * np.log(2) - 119760 * Nu * np.log(2) + 26880 * params.eta2 * np.log(2) + 
                18144 * kappa1 * params.S1z2 * np.log(2) + 18144 * delta * kappa1 * params.S1z2 * np.log(2) - 
                36288 * kappa1 * Nu * params.S1z2 * np.log(2) + 72576 * Nu * S1z * S2z * np.log(2) + 
                18144 * kappa2 * params.S2z2 * np.log(2) - 18144 * delta * kappa2 * params.S2z2 * np.log(2) - 
                36288 * kappa2 * Nu * params.S2z2 * np.log(2) + 12 * (-2173 - 4990 * Nu + 1120 * params.eta2 + 
                756 * kappa1 * (1 + delta - 2 * Nu) * params.S1z2 + 3024 * Nu * S1z * S2z - 
                756 * kappa2 * (-1 + delta + 2 * Nu) * params.S2z2) * np.log(b0) + 18 * (-2173 - 4990 * Nu + 
                1120 * params.eta2 + 756 * kappa1 * (1 + delta - 2 * Nu) * params.S1z2 + 3024 * Nu * S1z * S2z - 
                756 * kappa2 * (-1 + delta + 2 * Nu) * params.S2z2) * np.log(x)))

    # 4PN non-spinning quasi-circular (2,2) mode has been obtained from Blanchet
    # et al. arXiv:2304.11185. This is written in terms of phi. Please refer to file shared by Quentin.

    elif vpnorder == 8:
        return (-1.3109423540262542e-11 * (params.x5 * (29059430400 * EulerGamma * EulerGamma * (-107 + 13 * Nu) + 
                12 * (628830397253 + 1854914893791 * Nu + 421984442880 * np.log(2)) + 415134720 * EulerGamma * (6099 - 40277 * Nu + complex(0, 70) * (107 - 13 * Nu) * M_PI + 280 * (-107 + 13 * Nu) * np.log(2)) + 
                35 * (-28 * params.eta2 * (5385456111 + 5 * Nu * (-163158374 + 26251249 * Nu)) + 
                135135 * (54784 + 5 * Nu * (1951 + 6560 * Nu)) * M_PI2 - 955450349568 * Nu * np.log(2) + 
                3321077760 * (-107 + 13 * Nu) * np.log(2) * np.log(2) - complex(0, 5930496) * M_PI * (6099 - 74773 * Nu + 280 * (-107 + 13 * Nu) * np.log(2))) - 5700491596800 * np.log(mass) + 
                6966444925440 * np.log(x) + 69189120 * (420 * (-107 + 13 * Nu) * np.log(b0) * np.log(b0) + 
                420 * (-107 + 13 * Nu) * np.log(mass) * np.log(mass) - 14 * np.log(mass) * (-11803 * Nu + 
                60 * EulerGamma * (-107 + 13 * Nu) + complex(0, 30) * (107 - 13 * Nu) * M_PI + 
                120 * (-107 + 13 * Nu) * np.log(2) + 90 * (-107 + 13 * Nu) * np.log(x)) + 14 * np.log(b0) * (5885 - 6420 * EulerGamma - 11803 * Nu + 780 * EulerGamma * Nu + complex(0, 3210) * M_PI - 
                complex(0, 390) * Nu * M_PI + 120 * (-107 + 13 * Nu) * np.log(2) + (6420 - 780 * Nu) * np.log(mass) + 
                90 * (-107 + 13 * Nu) * np.log(x)) + 5 * np.log(x) * (-74149 * Nu + 252 * EulerGamma * (-107 + 13 * Nu) - 
                complex(0, 126) * (-107 + 13 * Nu) * M_PI + 504 * (-107 + 13 * Nu) * np.log(2) + 
                189 * (-107 + 13 * Nu) * np.log(x)) + 84672 * Nu * np.log(x0)))))
        
        # return ((params.x5 * (276756480 * EulerGamma * (11449 + 19105 * Nu) - 12 * (846557506853 + 1008017482431 * Nu) + 35 * (28 * params.eta2 * (5385456111 + 5 * Nu * (-163158374 + 26251249 * Nu)) - complex(0, 3953664) * (11449 + 109657 * Nu) * M_PI - 135135 * (54784 + 5 * Nu * (1951 + 6560 * Nu)) * M_PI2) + 138378240 * (11449 + 19105 * Nu) * np.log(16 * x))) / 7.62810048e10)

    else:
        return complex(0, 0)
    
def hl_2_m_2(mass: float, Nu: float, r: float, rDOT: float, Phi: float,
             PhiDOT: float, R: float, vpnorder: int, S1z: float,
             S2z: float, x: float, params: KeplerVars) -> complex:

    if vpnorder < 0 or vpnorder > 8:
        raise ValueError("Error in hl_2_m_2: Input PN order parameter should be between [0, 8].")

    else:
        # Calculate the leading amplitude coefficient
        amplitude = (4 * mass * Nu * np.sqrt(M_PI / 5.0)) / R
        
        # Sum the Generalized Orbital (GO) and Quasi-Circular (QC) terms
        waveform_modes = (
            hGO_2_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params) +
            hQC_2_m_2(mass, Nu, vpnorder, x, S1z, S2z, params)
        )
        
        # cpolar(r, theta) in C is equivalent to rect(r, theta) in Python
        phase_factor = rect(1.0, -2 * Phi)
        
        return amplitude * waveform_modes * phase_factor
    
def hl_2_m_min2(mass: float, Nu: float, r: float, rDOT: float, Phi: float,
                PhiDOT: float, R: float, vpnorder: int, S1z: float,
                S2z: float, x: float, params: KeplerVars) -> complex:

    if vpnorder < 0 or vpnorder > 8:
        raise ValueError("Error in hl_2_m_min2: Input PN order parameter should be between [0, 8].")

    else:
        # Calculate the leading amplitude coefficient
        amplitude = (4 * mass * Nu * np.sqrt(M_PI / 5.0)) / R
        
        # Sum the GO and QC terms, then apply the complex conjugate for the m = -2 mode
        waveform_modes = (
            hGO_2_m_2(mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params) +
            hQC_2_m_2(mass, Nu, vpnorder, x, S1z, S2z, params)
        ).conjugate()
        
        # Calculate the phase factor with positive 2 * Phi
        phase_factor = rect(1.0, 2 * Phi)
        
        return amplitude * waveform_modes * phase_factor

# H21

def hGO_2_m_1(
                mass: float,
                Nu: float,
                r: float,
                rDOT: float,
                PhiDOT: float,
                vpnorder: int,
                S1z: float,
                S2z: float,
                x: float,
                params: KeplerVars,
            ) -> complex:
    
    delta = np.sqrt(1 - 4 * Nu)
    kappa1 = 1.0
    kappa2 = 1.0
    r0 = 2 * mass / np.exp(0.5)

    if vpnorder == 1:
        return 0.6666666666666666j * delta * mass * PhiDOT

    elif vpnorder == 2:
        return (
            (-0.5j * params.Mtot2 *
             ((1 + delta) * S1z + (-1 + delta) * S2z))
            / params.r2
        )

    elif vpnorder == 3:
        return (
            (0.023809523809523808j * delta * mass * PhiDOT *
             (4 * mass * (-9 + 11 * Nu) +
              r * ((19 - 24 * Nu) * params.PhiDOT2 * params.r2 +
                   2j * (83 + 2 * Nu) * PhiDOT * r * rDOT +
                   2 * (-33 + 10 * Nu) * params.rDOT2)))
            / r
        )

    elif vpnorder == 4:
        return (
            (0.011904761904761904j * params.Mtot2 *
             (2 * mass *
              ((77 + 59 * Nu + 11 * delta * (7 + Nu)) * S1z +
               (-77 - 59 * Nu + 11 * delta * (7 + Nu)) * S2z) +
              r * (-2j * PhiDOT * r * rDOT *
                   (147 * (1 + delta) * S1z + (-83 + 13 * delta) * Nu * S1z +
                    147 * (-1 + delta) * S2z + (83 + 13 * delta) * Nu * S2z) +
                   params.rDOT2 *
                   ((105 * (1 + delta) - 4 * (13 + 15 * delta) * Nu) * S1z +
                    (-105 + 15 * delta * (7 - 4 * Nu) + 52 * Nu) * S2z) +
                   4 * params.PhiDOT2 * params.r2 *
                   ((-21 - 21 * delta + 66 * Nu + 4 * delta * Nu) * S1z +
                    (21 - 21 * delta - 66 * Nu + 4 * delta * Nu) * S2z))))
            / params.r3
        )

    elif vpnorder == 5:
        term1 = (
            (0.0013227513227513227j * delta * mass * PhiDOT *
             (10 * params.Mtot2 * (31 - 205 * Nu + 111 * params.eta2) -
              2 * mass * r *
              ((-197 + 5 * Nu + 660 * params.eta2) * params.PhiDOT2 * params.r2 +
               1j * (-3167 - 5278 * Nu + 201 * params.eta2) * PhiDOT * r * rDOT +
               8 * (202 + 587 * Nu - 177 * params.eta2) * params.rDOT2) +
              3 * params.r2 *
              ((152 - 692 * Nu + 333 * params.eta2) * params.PhiDOT4 * params.r4 +
               2j * (308 - 1607 * Nu + 111 * params.eta2) * params.PhiDOT3 * params.r3 * rDOT -
               3 * (75 - 560 * Nu + 68 * params.eta2) * params.PhiDOT2 * params.r2 * params.rDOT2 -
               2j * (-265 + 526 * Nu + 18 * params.eta2) * PhiDOT * r * params.rDOT3 +
               (-241 + 550 * Nu - 264 * params.eta2) * params.rDOT4)))
            / params.r2
        )

        # Henry et al. ecc spin terms
        term2 = (
            (params.Mtot3 *
             (-4 * rDOT *
              (kappa1 * (1 + delta - 2 * Nu) * params.S1z2 +
               kappa2 * (-1 + delta + 2 * Nu) * params.S2z2) -
              1j * PhiDOT * r *
              ((-((1 + delta) * (9 + kappa1)) +
                2 * (9 + (4 + 3 * delta) * kappa1) * Nu) * params.S1z2 -
               12 * delta * Nu * S1z * S2z +
               (9 + kappa2 - 2 * (9 + 4 * kappa2) * Nu +
                delta * (-9 - kappa2 + 6 * kappa2 * Nu)) * params.S2z2)))
            / (6.0 * params.r3)
        )

        return term1 + term2

    elif vpnorder == 6:
        term1 = (
            (delta * params.Mtot2 * Nu * PhiDOT *
             (mass * (195 * PhiDOT * r - 946j * rDOT) +
              9 * r *
              (270 * params.PhiDOT3 * params.r3 -
               483j * params.PhiDOT2 * params.r2 * rDOT -
               580 * PhiDOT * r * params.rDOT2 +
               42j * params.rDOT3)))
            / (315.0 * params.r2)
        )

        # Henry et al. ecc spin terms
        term2 = (
            0.00033068783068783067j * params.Mtot2 *
            (3 * params.r2 *
             (4j * PhiDOT * r * params.rDOT3 *
              ((-315 * (1 + delta) + 2 * (251 + 463 * delta) * Nu +
                4 * (15 + delta) * params.eta2) * S1z +
               (315 - 315 * delta - 502 * Nu + 926 * delta * Nu +
                4 * (-15 + delta) * params.eta2) * S2z) +
              12 * params.PhiDOT2 * params.r2 * params.rDOT2 *
              ((189 * (1 + delta) - 2 * (521 + 293 * delta) * Nu +
                7 * (55 + 23 * delta) * params.eta2) * S1z +
               (189 * (-1 + delta) + 2 * (521 - 293 * delta) * Nu +
                7 * (-55 + 23 * delta) * params.eta2) * S2z) +
              params.rDOT4 *
              ((567 * (1 + delta) - 16 * (77 + 64 * delta) * Nu +
                8 * (177 + 173 * delta) * params.eta2) * S1z +
               (567 * (-1 + delta) + 16 * (77 - 64 * delta) * Nu +
                8 * (-177 + 173 * delta) * params.eta2) * S2z) -
              4j * params.PhiDOT3 * params.r3 * rDOT *
              ((936 * (1 + delta) - 5 * (979 + 215 * delta) * Nu +
                2 * (1353 + 293 * delta) * params.eta2) * S1z +
               (936 * (-1 + delta) + 5 * (979 - 215 * delta) * Nu +
                2 * (-1353 + 293 * delta) * params.eta2) * S2z) +
              4 * params.PhiDOT4 * params.r4 *
              ((-252 * (1 + delta) + (1315 + 857 * delta) * Nu +
                4 * (-285 + 43 * delta) * params.eta2) * S1z +
               (252 + 5 * Nu * (-263 + 228 * Nu) +
                delta * (-252 + Nu * (857 + 172 * Nu))) * S2z)) -
             2 * mass * r *
             (-1j * PhiDOT * r * rDOT *
              ((2043 * (1 + delta) + (37 + 2597 * delta) * Nu +
                (10635 + 139 * delta) * params.eta2) * S1z +
               (2043 * (-1 + delta) + (-37 + 2597 * delta) * Nu +
                (-10635 + 139 * delta) * params.eta2) * S2z) +
              params.PhiDOT2 * params.r2 *
              ((-765 - Nu * (667 + 7773 * Nu) +
                delta * (-765 + 7 * Nu * (-533 + 245 * Nu))) * S1z +
               (765 + Nu * (667 + 7773 * Nu) +
                delta * (-765 + 7 * Nu * (-533 + 245 * Nu))) * S2z) +
              4 * params.rDOT2 *
              ((-234 * (1 + delta) - 4 * (560 + 901 * delta) * Nu +
                (483 + 1111 * delta) * params.eta2) * S1z +
               (234 + 7 * (320 - 69 * Nu) * Nu +
                delta * (-234 + Nu * (-3604 + 1111 * Nu))) * S2z)) +
             2 * params.Mtot2 *
             (1134 * kappa1 * (-1 - delta + (3 + delta) * Nu) * params.S1z3 +
              1134 * (1 + delta) * (-2 + kappa1) * Nu * params.S1z2 * S2z +
              S1z *
              (-5661 - 5661 * delta - 17156 * Nu - 9172 * delta * Nu +
               231 * params.eta2 + 775 * delta * params.eta2 +
               1134 * (-1 + delta) * (-2 + kappa2) * Nu * params.S2z2) +
              S2z * (5661 - 5661 * delta + 17156 * Nu - 9172 * delta * Nu -
                     231 * params.eta2 + 775 * delta * params.eta2 +
                     1134 * kappa2 * (1 - delta + (-3 + delta) * Nu) *
                     params.S2z2)))
            ) / params.r4
        

        return term1 + term2

    elif vpnorder == 7:

        spin_terms = (
            -23760 * params.Mtot3 *
            (params.PhiDOT3 * params.r4 *
             (-12j * (-35 + 107 * Nu) * S1z +
              5 * (126 - 462 * Nu + 80 * params.eta2 +
                   kappa1 * (-2 - 79 * Nu + 153 * params.eta2)) * params.S1z2 +
              S2z * (-420j + 1284j * Nu +
                     10 * (-63 + kappa2) * S2z +
                     5 * (462 + 79 * kappa2) * Nu * S2z -
                     5 * (80 + 153 * kappa2) * params.eta2 * S2z)) -
             1j * params.PhiDOT2 * params.r3 * rDOT *
             (-24j * (-35 + 202 * Nu) * S1z +
              5 * (-14 * Nu * (3 + 58 * Nu) +
                   kappa1 * (-409 + 1096 * Nu + 6 * params.eta2)) * params.S1z2 +
              S2z * (-840j + 2045 * kappa2 * S2z +
                     10 * (406 - 3 * kappa2) * params.eta2 * S2z +
                     Nu * (4848j + 210 * S2z - 5480 * kappa2 * S2z))) -
             4j * r * params.rDOT3 *
             (3j * (26 + 57 * Nu) * S1z +
              5 * (4 * params.eta2 +
                   kappa1 * (-7 + 49 * Nu + 15 * params.eta2)) * params.S1z2 -
              S2z * (78j - 35 * kappa2 * S2z +
                     5 * (4 + 15 * kappa2) * params.eta2 * S2z +
                     Nu * (171j + 245 * kappa2 * S2z))) +
             2j * mass * rDOT *
             (-4j * (-78 + 769 * Nu) * S1z +
              5 * ((14 - 59 * Nu) * Nu +
                   kappa1 * (-70 - 77 * Nu + 18 * params.eta2)) * params.S1z2 +
              S2z * (-312j + 350 * kappa2 * S2z +
                     5 * (59 - 18 * kappa2) * params.eta2 * S2z +
                     Nu * (3076j - 70 * S2z + 385 * kappa2 * S2z))) +
             PhiDOT * params.r2 * params.rDOT2 *
             (6j * (-62 + 219 * Nu) * S1z -
              5 * (189 - 756 * Nu + 388 * params.eta2 +
                   kappa1 * (98 - 266 * Nu + 276 * params.eta2)) * params.S1z2 +
              S2z * (372j + 945 * S2z + 490 * kappa2 * S2z +
                     20 * (97 + 69 * kappa2) * params.eta2 * S2z -
                     2 * Nu * (657j + 35 * (54 + 19 * kappa2) * S2z))) +
             mass * PhiDOT * r *
             (8j * (-61 + 480 * Nu) * S1z +
              5 * (-392 + 448 * Nu + 474 * params.eta2 +
                   kappa1 * (-11 + 150 * Nu + 58 * params.eta2)) * params.S1z2 -
              S2z * (-488j - 5 * (392 + 11 * kappa2) * S2z +
                     10 * (237 + 29 * kappa2) * params.eta2 * S2z +
                     10 * Nu * (384j + (224 + 75 * kappa2) * S2z))))
        )

        orbital_terms = (
            delta * mass *
            (-240 * mass * PhiDOT * params.r3 *
             ((197936 - 139360 * Nu - 367105 * params.eta2 +
               253245 * params.eta3) * params.PhiDOT4 * params.r4 +
              1j * (279236 - 483940 * Nu - 2817805 * params.eta2 +
                    459180 * params.eta3) * params.PhiDOT3 * params.r3 * rDOT -
              6 * (38627 + 89295 * Nu - 492740 * params.eta2 +
                   75975 * params.eta3) * params.PhiDOT2 * params.r2 * params.rDOT2 -
              1j * (-731008 + 2287930 * Nu + 981060 * params.eta2 +
                    10275 * params.eta3) * PhiDOT * r * params.rDOT3 +
              (-327667 + 436705 * Nu + 659790 * params.eta2 -
               438255 * params.eta3) * params.rDOT4) +
             900 * PhiDOT * params.r4 *
             (2 * (-2594 + 27609 * Nu - 74032 * params.eta2 +
                   25974 * params.eta3) * params.PhiDOT6 * params.r6 +
              4j * (-5730 + 58833 * Nu - 137842 * params.eta2 +
                    17123 * params.eta3) * params.PhiDOT5 * params.r5 * rDOT +
              2 * (-114 - 41622 * Nu + 147569 * params.eta2 +
                   4196 * params.eta3) * params.PhiDOT4 * params.r4 * params.rDOT2 +
              4j * (-9554 + 70788 * Nu - 156227 * params.eta2 +
                    5810 * params.eta3) * params.PhiDOT3 * params.r3 * params.rDOT3 +
              (17619 - 138450 * Nu + 322600 * params.eta2 -
               80816 * params.eta3) * params.PhiDOT2 * params.r2 * params.rDOT4 -
              2j * (8793 - 52230 * Nu + 69340 * params.eta2 +
                    2536 * params.eta3) * PhiDOT * r * params.rDOT5 +
              2 * (3957 - 24534 * Nu + 42584 * params.eta2 -
                   20800 * params.eta3) * params.rDOT6) -
             2 * params.Mtot3 *
             (-23760j * rDOT *
              (5 * (Nu * (-14 + 31 * Nu) + 7 * kappa1 * (10 + 31 * Nu)) * params.S1z2 +
               2 * S1z * (-156j + 155 * params.eta2 * S2z +
                          2 * Nu * (613j + 390 * S2z)) +
               S2z * (-312j + 350 * kappa2 * S2z +
                      155 * params.eta2 * S2z +
                      Nu * (2452j + 35 * (-2 + 31 * kappa2) * S2z))) +
              PhiDOT * r *
              (8946400 * params.eta3 -
               8 * (6991786 + 724680j * S1z +
                    7425 * (392 + 11 * kappa1) * params.S1z2 +
                    724680j * S2z +
                    7425 * (392 + 11 * kappa2) * params.S2z2) -
               3600 * params.eta2 *
               (-628 + 33 * (-19 + 92 * kappa1) * params.S1z2 -
                7326 * S1z * S2z +
                33 * (-19 + 92 * kappa2) * params.S2z2) +
               15 * Nu *
               (994455 * M_PI2 +
                8 * (-2249485 +
                     7920 * (-21 + 8 * kappa1) * params.S1z2 +
                     283536j * S2z +
                     7920 * (-21 + 8 * kappa2) * params.S2z2 -
                     1584 * S1z * (-179j + 170 * S2z))))) +
             3 * params.Mtot2 * r *
             (31680j * params.rDOT3 *
              (5 * (4 * params.eta2 + 7 * kappa1 * (-1 + 5 * Nu)) * params.S1z2 +
               S1z * (78j + 40 * params.eta2 * S2z +
                      Nu * (327j + 420 * S2z)) +
               S2z * (78j - 35 * kappa2 * S2z +
                      20 * params.eta2 * S2z +
                      Nu * (327j + 175 * kappa2 * S2z))) -
              22 * PhiDOT * r * params.rDOT2 *
              (2553200 * params.eta3 -
               24 * (268267 + 5580j * S1z +
                     525 * (27 + 14 * kappa1) * params.S1z2 +
                     5580j * S2z +
                     525 * (27 + 14 * kappa2) * params.S2z2) -
               200 * params.eta2 *
               (39445 + 72 * (-4 + 21 * kappa1) * params.S1z2 -
                3600 * S1z * S2z +
                72 * (-4 + 21 * kappa2) * params.S2z2) +
               25 * Nu *
               (23247 * M_PI2 +
                8 * (-69259 + 1026j * S1z +
                     126 * (27 + 5 * kappa1) * params.S1z2 +
                     1026j * S2z +
                     126 * (27 + 5 * kappa2) * params.S2z2))) +
              params.PhiDOT3 * params.r3 *
              (10071200 * params.eta3 +
               96 * (-421183 - 34650j * S1z +
                     825 * (-63 + kappa1) * params.S1z2 -
                     34650j * S2z +
                     825 * (-63 + kappa2) * params.S2z2) -
               400 * params.eta2 *
               (64177 + 792 * (-5 + 6 * kappa1) * params.S1z2 -
                17424 * S1z * S2z +
                792 * (-5 + 6 * kappa2) * params.S2z2) +
               15 * Nu *
               (426195 * M_PI2 +
                8 * (-509635 +
                     330 * (210 + 83 * kappa1) * params.S1z2 +
                     29304j * S2z +
                     330 * (210 + 83 * kappa2) * params.S2z2 -
                     792 * S1z * (-37j + 70 * S2z)))) -
              2j * params.PhiDOT2 * params.r2 * rDOT *
              (-8330400 * params.eta3 +
               8 * (-2810116 - 415800j * S1z +
                    1012275 * kappa1 * params.S1z2 -
                    415800j * S2z +
                    1012275 * kappa2 * params.S2z2) +
               4800 * params.eta2 *
               (13411 + 33 * (19 + 12 * kappa1) * params.S1z2 +
                462 * S1z * S2z +
                33 * (19 + 12 * kappa2) * params.S2z2) +
               5 * Nu *
               (1278585 * M_PI2 -
                8 * (5139685 +
                     990 * (-21 + 139 * kappa1) * params.S1z2 -
                     313632j * S2z +
                     990 * (-21 + 139 * kappa2) * params.S2z2 -
                     3564 * S1z * (88j + 185 * S2z))))))
        )

        log_term = (
            -13559040 * delta * params.Mtot3 * PhiDOT * r *
            (2 * mass - 3 * params.PhiDOT2 * params.r3 +
             6j * PhiDOT * params.r2 * rDOT +
             6 * r * params.rDOT2) *
            np.log(r / r0)
        )

        return (
            -1.0020843354176688e-7j *
            (spin_terms + orbital_terms + log_term)
        ) / params.r4

    else:
        return 0.0 + 0.0j
    
def hQC_2_m_1(
                mass: float, 
                Nu: float, 
                vpnorder: int, 
                x: float, 
                S1z: float, 
                S2z: float, 
                params: KeplerVars
            ) -> complex:
    
    delta: float = np.sqrt(1.0 - 4.0 * Nu)
    EulerGamma: float = 0.5772156649015329
    b0: float = 2.0 * mass / np.exp(0.5)
    r0: float = b0

    # Common log terms to simplify the messy 4PN-7PN expressions
    log_term_base = (-7.0 + 6.0 * EulerGamma - 3j * M_PI + np.log(64.0) + 6.0 * np.log(b0) + 9.0 * np.log(x))

    if vpnorder == 4:
        return (-2.0 * delta * params.x3 * log_term_base) / 9.0

    elif vpnorder == 5:
        spin_factor = ((1.0 + delta) * S1z + (-1.0 + delta) * S2z)
        return (spin_factor * params.x3p5 * log_term_base) / 6.0

    elif vpnorder == 6:
        return -0.007936507936507936 * (delta * (-17.0 + 6.0 * Nu) * params.x4 * log_term_base)

    elif vpnorder == 7:
        # Heavily nested 3.5PN (order 7) term
        term1 = (params.x4p5 * (
            -98 * S1z + 84 * EulerGamma * S1z + 3017 * Nu * S1z - 2586 * EulerGamma * Nu * S1z - 
            42j * M_PI * S1z + 1293j * Nu * M_PI * S1z + 98 * S2z - 84 * EulerGamma * S2z - 
            3017 * Nu * S2z + 2586 * EulerGamma * Nu * S2z + 42j * M_PI * S2z - 1293j * Nu * M_PI * S2z + 
            84 * S1z * np.log(2) - 2586 * Nu * S1z * np.log(2) - 84 * S2z * np.log(2) + 2586 * Nu * S2z * np.log(2) + 
            84 * S1z * np.log(b0) - 2586 * Nu * S1z * np.log(b0) - 84 * S2z * np.log(b0) + 2586 * Nu * S2z * np.log(b0) + 
            126 * S1z * np.log(x) - 3879 * Nu * S1z * np.log(x) - 126 * S2z * np.log(x) + 3879 * Nu * S2z * np.log(x) + 
            252 * delta * (
                1.5875434618291762j + 1.7523809523809524j * EulerGamma - 1.3333333333333333j * EulerGamma**2 + 
                (92 * M_PI) / 105.0 - (4 * EulerGamma * M_PI) / 3.0 + 0.1111111111111111j * M_PI**2 - 
                (7 * S1z) / 18.0 + (EulerGamma * S1z) / 3.0 + (29 * Nu * S1z) / 12.0 - (29 * EulerGamma * Nu * S1z) / 14.0 - 
                0.16666666666666666j * M_PI * S1z + 1.0357142857142858j * Nu * M_PI * S1z - 
                (7 * S2z) / 18.0 + (EulerGamma * S2z) / 3.0 + (29 * Nu * S2z) / 12.0 - (29 * EulerGamma * Nu * S2z) / 14.0 - 
                0.16666666666666666j * M_PI * S2z + 1.0357142857142858j * Nu * M_PI * S2z + 
                1.7523809523809524j * np.log(2) - 2.6666666666666665j * EulerGamma * np.log(2) - 
                (4 * M_PI * np.log(2)) / 3.0 + (S1z * np.log(2)) / 3.0 - (29 * Nu * S1z * np.log(2)) / 14.0 + 
                (S2z * np.log(2)) / 3.0 - (29 * Nu * S2z * np.log(2)) / 14.0 - 1.3333333333333333j * np.log(2)**2 - 
                1.3333333333333333j * np.log(b0)**2 - 1.3587301587301588j * np.log(r0) + 
                (np.log(b0) * (392j - 336j * EulerGamma - 168 * M_PI + 42 * S1z - 261 * Nu * S1z + 42 * S2z - 261 * Nu * S2z - 336j * np.log(2) - 504j * np.log(x))) / 126.0 + 
                2.6285714285714286j * np.log(x) - 4j * EulerGamma * np.log(x) - 2 * M_PI * np.log(x) + 
                (S1z * np.log(x)) / 2.0 - (87 * Nu * S1z * np.log(x)) / 28.0 + (S2z * np.log(x)) / 2.0 - 
                (87 * Nu * S2z * np.log(x)) / 28.0 - 4j * np.log(2) * np.log(x) - 3j * np.log(x)**2
            )
        )) / 252.0
        return term1

    else:
        return 0.0 + 0.0j

def hl_2_m_1(
                mass: float, 
                Nu: float, 
                r: float, 
                rDOT: float, 
                Phi: float, 
                PhiDOT: float, 
                R: float, 
                vpnorder: int, 
                S1z: float, 
                S2z: float, 
                x: float, 
                params: KeplerVars
            ) -> complex:
    
    if vpnorder > 8:
        raise ValueError(
            "Error in hl_2_m_1: Input PN order parameter should be between [0, 8]."
        )
    
    # Calculate the pre-factor: (4 * M * eta * sqrt(pi/5)) / R
    pre_factor: float = (4.0 * mass * Nu * np.sqrt(M_PI / 5.0)) / R
    
    # Evaluate the General Orbit (GO) and Quasi-Circular (QC) components
    go_component: complex = hGO_2_m_1(
        mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params
    )
    qc_component: complex = hQC_2_m_1(
        mass, Nu, vpnorder, x, S1z, S2z, params
    )
    
    phase_factor: complex = np.exp(-1j * Phi)
    
    # Combine terms
    return pre_factor * (go_component + qc_component) * phase_factor

def hl_2_m_min1(
                    mass: float, 
                    Nu: float, 
                    r: float, 
                    rDOT: float, 
                    Phi: float, 
                    PhiDOT: float, 
                    R: float, 
                    vpnorder: int, 
                    S1z: float, 
                    S2z: float, 
                    x: float, 
                    params: KeplerVars
                ) -> complex:

    if vpnorder > 8:
        raise ValueError(
            "Error in hl_2_m_min1: Input PN order parameter should be between [0, 8]."
        )

    # Pre-factor: (4 * M * eta * sqrt(pi/5)) / R
    pre_factor: float = (4.0 * mass * Nu * np.sqrt(M_PI / 5.0)) / R

    # Evaluate the General Orbit (GO) and Quasi-Circular (QC) components
    # These are the same functions used for the m=1 mode
    go_component: complex = hGO_2_m_1(
        mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params
    )
    qc_component: complex = hQC_2_m_1(
        mass, Nu, vpnorder, x, S1z, S2z, params
    )

    # In C: conj(hGO + hQC)
    # In Python: use the .conjugate() method or np.conj()
    amplitude_term: complex = (go_component + qc_component).conjugate()

    # In C: cpolar(1, 1 * Phi) -> exp(i * Phi)
    # Note the sign change from the m=1 mode
    phase_factor: complex = np.exp(1j * Phi)

    return pre_factor * amplitude_term * phase_factor

# H33

def hGO_3_m_3(
                mass: float,
                Nu: float,
                r: float,
                rDOT: float,
                PhiDOT: float,
                vpnorder: int,
                S1z: float,
                S2z: float,
                x: float,
                params: KeplerVars,
            ) -> complex:
    delta = np.sqrt(1 - 4 * Nu)
    kappa1 = 1.0
    kappa2 = 1.0
    r0 = 2 * mass / np.exp(0.5)

    combination_a  = 1j * PhiDOT * r - rDOT
    combination_a3 = combination_a ** 3

    combination_b  = PhiDOT * r + 1j * rDOT
    combination_b2 = combination_b ** 2
    combination_b4 = combination_b2 ** 2
    combination_b6 = combination_b ** 6

    combination_c  = PhiDOT * r - 1j * rDOT
    combination_c2 = combination_c ** 2

    combination_d  = -1j * PhiDOT * r + rDOT
    combination_d2 = combination_d ** 2
    combination_d5 = combination_d ** 5

    combination_e  = 1j * PhiDOT * r + rDOT
    combination_e3 = combination_e ** 3

    if vpnorder == 1:
        return (
            (np.sqrt(0.11904761904761904) * delta *
             (2 * r * combination_a3 +
              mass * (-7j * PhiDOT * r + 4 * rDOT)))
            / (2.0 * r)
        )

    elif vpnorder == 3:
        return (
            (np.sqrt(0.11904761904761904) * delta *
             (6 * (-5 + 19 * Nu) * params.r2 * combination_b4 *
              (1j * PhiDOT * r + rDOT) +
              2 * params.Mtot2 *
              (-3j * (-101 + 43 * Nu) * PhiDOT * r +
               (-109 + 86 * Nu) * rDOT) +
              3 * mass * r *
              (-12j * (1 + 4 * Nu) * params.PhiDOT3 * params.r3 +
               6 * (14 + 31 * Nu) * params.PhiDOT2 * params.r2 * rDOT +
               3j * (33 + 62 * Nu) * PhiDOT * r * params.rDOT2 -
               4 * (8 + 17 * Nu) * params.rDOT3)))
            / (36.0 * params.r2)
        )

    elif vpnorder == 4:
        return (
            (-0.125j * np.sqrt(0.11904761904761904) * params.Mtot2 *
             (4 * mass * (-1 + 5 * Nu) * ((1 + delta) * S1z + (-1 + delta) * S2z) +
              r * (2 * params.rDOT2 *
                   (6 * (1 + delta) * S1z - 5 * (5 + 3 * delta) * Nu * S1z +
                    (-6 + delta * (6 - 15 * Nu) + 25 * Nu) * S2z) +
                   params.PhiDOT2 * params.r2 *
                   (-24 * (1 + delta) * S1z + (119 + 33 * delta) * Nu * S1z +
                    (24 - 119 * Nu + 3 * delta * (-8 + 11 * Nu)) * S2z) +
                   2j * PhiDOT * r * rDOT *
                   (-18 * (1 + delta) * S1z + (77 + 39 * delta) * Nu * S1z +
                    (18 - 77 * Nu + 3 * delta * (-6 + 13 * Nu)) * S2z))))
            / params.r3
        )

    elif vpnorder == 5:
        term1 = (
            delta *
            (30 * (183 - 1579 * Nu + 3387 * params.eta2) * params.r3 *
             combination_c2 * combination_d5 +
             10 * params.Mtot3 *
             (-1j * (26473 - 27451 * Nu + 9921 * params.eta2) * PhiDOT * r +
              4 * (623 - 732 * Nu + 1913 * params.eta2) * rDOT) +
             2 * params.Mtot2 * r *
             (-11j * (-5353 - 13493 * Nu + 4671 * params.eta2) * params.PhiDOT3 * params.r3 +
              (-75243 - 142713 * Nu + 192821 * params.eta2) * params.PhiDOT2 * params.r2 * rDOT +
              220j * (-256 + 781 * Nu + 840 * params.eta2) * PhiDOT * r * params.rDOT2 -
              10 * (-756 + 8238 * Nu + 7357 * params.eta2) * params.rDOT3) +
             3 * mass * params.r2 *
             (2j * (-7633 + 9137 * Nu + 28911 * params.eta2) * params.PhiDOT5 * params.r5 -
              4 * (-8149 + 1576 * Nu + 43533 * params.eta2) * params.PhiDOT4 * params.r4 * rDOT -
              2j * (-9297 - 19517 * Nu + 64839 * params.eta2) * params.PhiDOT3 * params.r3 * params.rDOT2 -
              32 * (-1288 + 3667 * Nu + 4056 * params.eta2) * params.PhiDOT2 * params.r2 * params.rDOT3 -
              5j * (-9851 + 17954 * Nu + 40968 * params.eta2) * PhiDOT * r * params.rDOT4 +
              20 * (-771 + 1126 * Nu + 3616 * params.eta2) * params.rDOT5))
            / (1584.0 * np.sqrt(210) * params.r3)
        )

        # Henry et al. ecc spin terms
        term2 = (
            0.125j * np.sqrt(1.0714285714285714) * params.Mtot3 *
            (7 * PhiDOT * r + 2j * rDOT) *
            (kappa1 * (-1 - delta + 2 * (2 + delta) * Nu) * params.S1z2 +
             S2z * (-4 * delta * Nu * S1z +
                    kappa2 * (1 - delta + 2 * (-2 + delta) * Nu) * S2z))
            / params.r3
        )

        return term1 + term2

    elif vpnorder == 6:
        term1 = (
            -(delta * params.Mtot2 * Nu *
              (668 * params.Mtot2 +
               2 * mass * r *
               (4081 * params.PhiDOT2 * params.r2 +
                297j * PhiDOT * r * rDOT - 452 * params.rDOT2) +
               5 * params.r2 *
               (1329 * params.PhiDOT4 * params.r4 -
                2926j * params.PhiDOT3 * params.r3 * rDOT -
                384 * params.PhiDOT2 * params.r2 * params.rDOT2 -
                408j * PhiDOT * r * params.rDOT3 +
                200 * params.rDOT4)))
            / (36.0 * np.sqrt(210) * params.r4)
        )

        # Henry et al. ecc spin terms
        term2 = (
            -0.006944444444444444j * params.Mtot2 *
            (10 * params.Mtot2 *
             ((252 * (1 + delta) - (1277 + 1279 * delta) * Nu +
               8 * (12 + 47 * delta) * params.eta2) * S1z +
              (252 * (-1 + delta) + (1277 - 1279 * delta) * Nu +
               8 * (-12 + 47 * delta) * params.eta2) * S2z) +
             2 * mass * r *
             (2 * params.PhiDOT2 * params.r2 *
              ((1320 * (1 + delta) - 2 * (4469 + 211 * delta) * Nu +
                (8709 + 2777 * delta) * params.eta2) * S1z +
               (1320 * (-1 + delta) + 8938 * Nu - 422 * delta * Nu +
                (-8709 + 2777 * delta) * params.eta2) * S2z) +
              3j * PhiDOT * r * rDOT *
              ((2000 * (1 + delta) - (9147 + 3173 * delta) * Nu +
                (8911 + 5273 * delta) * params.eta2) * S1z +
               (2000 * (-1 + delta) + (9147 - 3173 * delta) * Nu +
                (-8911 + 5273 * delta) * params.eta2) * S2z) +
              10 * params.rDOT2 *
              ((-105 * (1 + delta) + (541 + 77 * delta) * Nu -
                2 * (462 + 247 * delta) * params.eta2) * S1z +
               (105 + Nu * (-541 + 924 * Nu) +
                delta * (-105 + (77 - 494 * Nu) * Nu)) * S2z)) -
             3 * params.r2 *
             (-3 * params.PhiDOT2 * params.r2 * params.rDOT2 *
              ((480 * (1 + delta) - (1711 + 1889 * delta) * Nu +
                2 * (-1161 + 757 * delta) * params.eta2) * S1z +
               (480 * (-1 + delta) + (1711 - 1889 * delta) * Nu +
                2 * (1161 + 757 * delta) * params.eta2) * S2z) +
              2 * params.PhiDOT4 * params.r4 *
              ((350 * (1 + delta) - 4 * (404 + 461 * delta) * Nu +
                (883 + 769 * delta) * params.eta2) * S1z +
               (350 * (-1 + delta) + 4 * (404 - 461 * delta) * Nu +
                (-883 + 769 * delta) * params.eta2) * S2z) +
              2j * params.PhiDOT3 * params.r3 * rDOT *
              ((660 * (1 + delta) - (4061 + 2899 * delta) * Nu +
                (2643 + 4789 * delta) * params.eta2) * S1z +
               (660 * (-1 + delta) + (4061 - 2899 * delta) * Nu +
                (-2643 + 4789 * delta) * params.eta2) * S2z) +
              10 * params.rDOT4 *
              ((-30 * (1 + delta) + (187 + 101 * delta) * Nu -
                2 * (159 + 61 * delta) * params.eta2) * S1z +
               (30 + Nu * (-187 + 318 * Nu) +
                delta * (-30 + (101 - 122 * Nu) * Nu)) * S2z) +
              2j * PhiDOT * r * params.rDOT3 *
              ((90 + Nu * (-1321 + 5118 * Nu) +
                delta * (90 + Nu * (-319 + 714 * Nu))) * S1z +
               (-90 + (1321 - 5118 * Nu) * Nu +
                delta * (90 + Nu * (-319 + 714 * Nu))) * S2z)))
            / (np.sqrt(210) * params.r4)
        )

        return term1 + term2

    elif vpnorder == 7:
        M_PI2 = M_PI ** 2

        spin_terms = (
            504504 * params.Mtot3 *
            (2 * mass *
             (rDOT *
              (S1z * (108 - 498 * Nu +
                      5j * (-24 - 5 * kappa1 + 3 * (76 + kappa1) * Nu +
                            4 * (-111 + 13 * kappa1) * params.eta2) * S1z) +
               6 * (-18 + 83 * Nu) * S2z -
               5j * (-24 - 5 * kappa2 + 3 * (76 + kappa2) * Nu +
                     4 * (-111 + 13 * kappa2) * params.eta2) * params.S2z2) +
              PhiDOT * r *
              (S1z * (-3j * (-99 + 581 * Nu) +
                      5 * (-24 + 399 * kappa1 + 48 * (7 - 19 * Nu) * Nu +
                           kappa1 * Nu * (-1629 + 188 * Nu)) * S1z) +
               3j * (-99 + 581 * Nu) * S2z -
               5 * (-24 + 399 * kappa2 + 48 * (7 - 19 * Nu) * Nu +
                    kappa2 * Nu * (-1629 + 188 * Nu)) * params.S2z2)) +
             r * (3 * PhiDOT * r * params.rDOT2 *
                  (S1z * (216j + 545 * kappa1 * S1z +
                          40 * (45 + 8 * kappa1) * params.eta2 * S1z -
                          30 * Nu * (50j + 20 * S1z + 73 * kappa1 * S1z)) +
                   12j * (-18 + 125 * Nu) * S2z -
                   5 * (109 * kappa2 - 6 * (20 + 73 * kappa2) * Nu +
                        8 * (45 + 8 * kappa2) * params.eta2) * params.S2z2) +
                  2 * params.rDOT3 *
                  (S1z * (-54 + 145j * kappa1 * S1z -
                          30j * Nu * (11j + 2 * Nu * S1z +
                                      kappa1 * (17 + 8 * Nu) * S1z)) +
                   6 * (9 - 55 * Nu) * S2z +
                   5j * (12 * params.eta2 +
                         kappa2 * (-29 + 6 * Nu * (17 + 8 * Nu))) * params.S2z2) +
                  6 * params.PhiDOT2 * params.r2 * rDOT *
                  (S1z * (297 - 465 * Nu +
                          5j * (6 * (20 - 87 * Nu) * Nu +
                                kappa1 * (-50 + 3 * Nu * (47 + 76 * Nu))) * S1z) +
                   3 * (-99 + 155 * Nu) * S2z -
                   5j * (6 * (20 - 87 * Nu) * Nu +
                         kappa2 * (-50 + 3 * Nu * (47 + 76 * Nu))) * params.S2z2) +
                  params.PhiDOT3 * params.r3 *
                  (3j * (531 - 1295 * Nu) * S1z +
                   10 * (-33 * kappa1 + 6 * (30 + 13 * kappa1) * Nu +
                         4 * (-96 + 67 * kappa1) * params.eta2) * params.S1z2 +
                   S2z * (-1593j + 330 * kappa2 * S2z +
                          5 * Nu * (777j -
                                    4 * (90 + 39 * kappa2 - 192 * Nu +
                                         134 * kappa2 * Nu) * S2z)))))
        )

        orbital_terms = (
            delta *
            (-17640 * (-4083 + Nu * (58311 + Nu * (-269240 + 405617 * Nu))) *
             params.r4 * combination_b6 * combination_e3 +
             168 * params.Mtot2 * params.r2 *
             (1j * (-7508635 + 7 * Nu * (-1318438 + Nu * (-10231834 + 9667755 * Nu))) *
              params.PhiDOT5 * params.r5 +
              7 * (1235591 + Nu * (884445 + (23935218 - 26913443 * Nu) * Nu)) *
              params.PhiDOT4 * params.r4 * rDOT -
              1j * (8961149 + 7 * Nu * (-31755709 + Nu * (-11134798 + 22187331 * Nu))) *
              params.PhiDOT3 * params.r3 * params.rDOT2 -
              (-36806435 + 7 * Nu * (33178545 + Nu * (24565078 + 22873537 * Nu))) *
              params.PhiDOT2 * params.r2 * params.rDOT3 -
              5j * (-7761899 + 7 * Nu * (2892563 + 5998602 * Nu + 7493619 * params.eta2)) *
              PhiDOT * r * params.rDOT4 +
              5 * (-2422057 + 7 * Nu * (501045 + Nu * (2033141 + 2771816 * Nu))) *
              params.rDOT5) +
             1764 * mass * params.r3 *
             (-2j * (239087 + Nu * (-1206515 + Nu * (422631 + 3979375 * Nu))) *
              params.PhiDOT7 * params.r7 +
              2 * (621284 + Nu * (-2279907 + 2 * Nu * (-1180187 + 5876531 * Nu))) *
              params.PhiDOT6 * params.r6 * rDOT +
              2j * (39270 + Nu * (1235486 - 5319747 * Nu + 4406349 * params.eta2)) *
              params.PhiDOT5 * params.r5 * params.rDOT2 +
              8 * (349111 + 4 * Nu * (-519370 + 33 * Nu * (10939 + 42635 * Nu))) *
              params.PhiDOT4 * params.r4 * params.rDOT3 +
              2j * (1212607 + 3 * Nu * (-2012698 - 67827 * Nu + 7955628 * params.eta2)) *
              params.PhiDOT3 * params.r3 * params.rDOT4 +
              4 * (201135 + 2 * Nu * (-773107 + Nu * (1214819 + 1157652 * Nu))) *
              params.PhiDOT2 * params.r2 * params.rDOT5 +
              5j * (333969 + 2 * Nu * (-981471 + 4 * Nu * (154039 + 750016 * Nu))) *
              PhiDOT * r * params.rDOT6 -
              40 * (13245 + 2 * Nu * (-37005 + Nu * (14251 + 130160 * Nu))) *
              params.rDOT7) +
             2 * params.Mtot4 *
             (4 * rDOT *
              (269279500 * params.eta3 +
               2 * (-174108226 +
                    63063 * S1z * (108j + 5 * (24 + 5 * kappa1) * S1z) +
                    63063 * S2z * (108j + 5 * (24 + 5 * kappa2) * S2z)) -
               21 * Nu *
               (103100846 + 1846845 * M_PI2 + 1693692j * S2z -
                6006 * (S1z * (-282j + 5 * (-180 + 7 * kappa1) * S1z) -
                        980 * S1z * S2z +
                        5 * (-180 + 7 * kappa2) * params.S2z2)) -
               2940 * params.eta2 *
               (-122855 +
                4719 * ((-6 + 7 * kappa1) * params.S1z2 -
                        26 * S1z * S2z +
                        (-6 + 7 * kappa2) * params.S2z2))) +
              1j * PhiDOT * r *
              (-1176172480 * params.eta3 +
               8 * (-74084729 +
                    189189 * S1z * (99j + 5 * (-8 + 133 * kappa1) * S1z) +
                    189189 * S2z * (99j + 5 * (-8 + 133 * kappa2) * S2z)) +
               35280 * params.eta2 *
               (56255 +
                429 * ((-22 + 65 * kappa1) * params.S1z2 -
                       174 * S1z * S2z +
                       (-22 + 65 * kappa2) * params.S2z2)) -
               147 * Nu *
               (-65012788 + 4485195 * M_PI2 + 3943368j * S2z +
                10296 * (S1z * (383j + 5 * (-96 + 277 * kappa1) * S1z) -
                         3220 * S1z * S2z +
                         5 * (-96 + 277 * kappa2) * params.S2z2)))) +
             params.Mtot3 * r *
             (-12j * PhiDOT * r * params.rDOT2 *
              (-1035895280 * params.eta3 -
               2 * (-547993687 +
                    63063 * S1z * (216j + 545 * kappa1 * S1z) +
                    63063 * S2z * (216j + 545 * kappa2 * S2z)) +
               77 * Nu *
               (42451610 + 1511055 * M_PI2 + 1749384j * S2z +
                6552 * (S1z * (267j + 25 * (6 + 11 * kappa1) * S1z) -
                        5 * S1z * S2z +
                        25 * (6 + 11 * kappa2) * params.S2z2)) +
               490 * params.eta2 *
               (-5802767 +
                5148 * ((-6 + 23 * kappa1) * params.S1z2 -
                        58 * S1z * S2z +
                        (-6 + 23 * kappa2) * params.S2z2))) +
              4 * params.rDOT3 *
              (-1359334480 * params.eta3 -
               4 * (-150254558 +
                    63063 * S1z * (54j + 145 * kappa1 * S1z) +
                    63063 * S2z * (54j + 145 * kappa2 * S2z)) +
               231 * Nu *
               (8490448 + 503685 * M_PI2 + 242424j * S2z +
                2184 * (S1z * (111j + 110 * kappa1 * S1z) +
                        70 * S1z * S2z +
                        110 * kappa2 * params.S2z2)) +
               11760 * params.eta2 *
               (-312980 +
                429 * ((3 + 25 * kappa1) * params.S1z2 -
                       44 * S1z * S2z +
                       (3 + 25 * kappa2) * params.S2z2))) +
              6 * params.PhiDOT2 * params.r2 * rDOT *
              (2368900688 * params.eta3 +
               8 * (-812986529 +
                    63063 * S1z * (297j + 250 * kappa1 * S1z) +
                    63063 * S2z * (297j + 250 * kappa2 * S2z)) -
               1176 * params.eta2 *
               (2423171 +
                4290 * ((-3 + 41 * kappa1) * params.S1z2 -
                        88 * S1z * S2z +
                        (-3 + 41 * kappa2) * params.S2z2)) +
               539 * Nu *
               (-24139772 + 647595 * M_PI2 + 120744j * S2z -
                936 * (S1z * (-129j + 600 * S1z + 205 * kappa1 * S1z) -
                       460 * S1z * S2z +
                       5 * (120 + 41 * kappa2) * params.S2z2))) +
              1j * params.PhiDOT3 * params.r3 *
              (-4538040136 * params.eta3 -
               88 * (259018351 +
                     17199 * S1z * (-531j + 110 * kappa1 * S1z) +
                     17199 * S2z * (-531j + 110 * kappa2 * S2z)) +
               2352 * params.eta2 *
               (7332973 +
                12870 * ((5 + 23 * kappa1) * params.S1z2 -
                         36 * S1z * S2z +
                         (5 + 23 * kappa2) * params.S2z2)) +
               21 * Nu *
               (49864815 * M_PI2 +
                8 * (-88128538 - 2099097j * S2z +
                     9009 * (S1z * (-233j + 40 * (15 + kappa1) * S1z) +
                             360 * S1z * S2z +
                             40 * (15 + kappa2) * params.S2z2))))))
        )

        log_term = (
            74954880 * delta * params.Mtot3 *
            (22j * mass * PhiDOT * r +
             59j * params.PhiDOT3 * params.r4 + 8 * mass * rDOT +
             66 * params.PhiDOT2 * params.r3 * rDOT +
             24j * PhiDOT * params.r2 * params.rDOT2 -
             4 * r * params.rDOT3) *
            np.log(r / r0)
        )

        return (
            1j * (spin_terms + orbital_terms + log_term)
            / (2.4216192e7 * np.sqrt(210) * params.r4)
        )

    else:
        return 0.0 + 0.0j
    
def hQC_3_m_3(
                mass: float, 
                Nu: float, 
                vpnorder: int, 
                x: float, 
                S1z: float, 
                S2z: float, 
                params: KeplerVars
            ) -> complex:
    
    delta: float = np.sqrt(1.0 - 4.0 * Nu)
    EulerGamma: float = 0.5772156649015329
    b0: float = 2.0 * mass / np.exp(0.5)
    r0: float = b0

    if vpnorder == 4:
        return (
            (3.0 * np.sqrt(0.04285714285714286) * delta * params.x3 * (
                -97.0 + 60.0 * EulerGamma - 30.0j * M_PI + 
                60.0 * np.log(2.0) + 60.0 * np.log(3.0) + 
                60.0 * np.log(b0) + 90.0 * np.log(x)
            )) / 4.0
        )

    elif vpnorder == 6:
        numerator_6 = (delta * params.x4 * (
            188568.0 - 116640.0 * EulerGamma - 70537.0 * Nu + 43740.0 * EulerGamma * Nu + 
            58320.0j * M_PI - 21870.0j * Nu * M_PI - 
            116640.0 * np.log(2.0) + 43740.0 * Nu * np.log(2.0) - 
            116640.0 * np.log(3.0) + 43740.0 * Nu * np.log(3.0) + 
            14580.0 * (-8.0 + 3.0 * Nu) * np.log(b0) + 
            21870.0 * (-8.0 + 3.0 * Nu) * np.log(x)
        ))
        return numerator_6 / (216.0 * np.sqrt(210.0))

    elif vpnorder == 7:
        # Logarithmic expansion for the 3.5PN (order 7) term
        # log_b0_term handles the final nested log(b0) block
        log_b0_term = 15876.0 * np.log(b0) * (
            -194.0j * delta + 120.0j * delta * EulerGamma + 60.0 * delta * M_PI + 
            5.0 * (-4.0 - 4.0 * delta + 19.0 * Nu + 5.0 * delta * Nu) * S1z + 
            20.0 * S2z - 20.0 * delta * S2z - 95.0 * Nu * S2z + 25.0 * delta * Nu * S2z + 
            120.0j * delta * np.log(2.0) + 120.0j * delta * np.log(3.0) + 
            180.0j * delta * np.log(x)
        )

        numerator_7 = (params.x4p5 * (
            1434564.0j * delta - 2490264.0j * delta * EulerGamma + 952560.0j * delta * EulerGamma**2 - 
            1245132.0 * delta * M_PI + 952560.0 * delta * EulerGamma * M_PI - 
            79380.0j * delta * (M_PI**2) + 513324.0 * S1z + 513324.0 * delta * S1z - 
            317520.0 * EulerGamma * S1z - 317520.0 * delta * EulerGamma * S1z - 
            2439857.0 * Nu * S1z - 643223.0 * delta * Nu * S1z + 1508220.0 * EulerGamma * Nu * S1z + 
            396900.0 * delta * EulerGamma * Nu * S1z + 158760.0j * M_PI * S1z + 158760.0j * delta * M_PI * S1z - 
            754110.0j * Nu * M_PI * S1z - 198450.0j * delta * Nu * M_PI * S1z - 
            513324.0 * S2z + 513324.0 * delta * S2z + 317520.0 * EulerGamma * S2z - 
            317520.0 * delta * EulerGamma * S2z + 2439857.0 * Nu * S2z - 643223.0 * delta * Nu * S2z - 
            1508220.0 * EulerGamma * Nu * S2z + 396900.0 * delta * EulerGamma * Nu * S2z - 
            158760.0j * M_PI * S2z + 158760.0j * delta * M_PI * S2z + 
            754110.0j * Nu * M_PI * S2z - 198450.0j * delta * Nu * M_PI * S2z - 
            2490264.0j * delta * np.log(2.0) + 1905120.0j * delta * EulerGamma * np.log(2.0) + 
            952560.0 * delta * M_PI * np.log(2.0) - 317520.0 * S1z * np.log(2.0) - 
            317520.0 * delta * S1z * np.log(2.0) + 1508220.0 * Nu * S1z * np.log(2.0) + 
            396900.0 * delta * Nu * S1z * np.log(2.0) + 317520.0 * S2z * np.log(2.0) - 
            317520.0 * delta * S2z * np.log(2.0) - 1508220.0 * Nu * S2z * np.log(2.0) + 
            396900.0 * delta * Nu * S2z * np.log(2.0) + 952560.0j * delta * (np.log(2.0)**2) - 
            2490264.0j * delta * np.log(3.0) + 1905120.0j * delta * EulerGamma * np.log(3.0) + 
            952560.0 * delta * M_PI * np.log(3.0) - 317520.0 * S1z * np.log(3.0) - 
            317520.0 * delta * S1z * np.log(3.0) + 1508220.0 * Nu * S1z * np.log(3.0) + 
            396900.0 * delta * Nu * S1z * np.log(3.0) + 317520.0 * S2z * np.log(3.0) - 
            317520.0 * delta * S2z * np.log(3.0) - 1508220.0 * Nu * S2z * np.log(3.0) + 
            396900.0 * delta * Nu * S2z * np.log(3.0) + 1905120.0j * delta * np.log(2.0) * np.log(3.0) + 
            952560.0j * delta * (np.log(3.0)**2) + 952560.0j * delta * (np.log(b0)**2) + 
            589680.0j * delta * np.log(r0) - 3735396.0j * delta * np.log(x) + 
            2857680.0j * delta * EulerGamma * np.log(x) + 1428840.0 * delta * M_PI * np.log(x) - 
            476280.0 * S1z * np.log(x) - 476280.0 * delta * S1z * np.log(x) + 
            2262330.0 * Nu * S1z * np.log(x) + 595350.0 * delta * Nu * S1z * np.log(x) + 
            476280.0 * S2z * np.log(x) - 476280.0 * delta * S2z * np.log(x) - 
            2262330.0 * Nu * S2z * np.log(x) + 595350.0 * delta * Nu * S2z * np.log(x) + 
            2857680.0j * delta * np.log(2.0) * np.log(x) + 
            2857680.0j * delta * np.log(3.0) * np.log(x) + 
            2143260.0j * delta * (np.log(x)**2) + 
            log_b0_term
        ))
        return numerator_7 / (2352.0 * np.sqrt(210.0))

    else:
        return 0.0 + 0.0j
    
def hl_3_m_3(
                mass: float, 
                Nu: float, 
                r: float, 
                rDOT: float, 
                Phi: float, 
                PhiDOT: float, 
                R: float, 
                vpnorder: int, 
                S1z: float, 
                S2z: float, 
                x: float, 
                params: KeplerVars
            ) -> complex:

    if vpnorder > 8:
        raise ValueError(
            "Error in hl_3_m_3: Input PN order parameter should be between [0, 8]."
        )

    # Pre-factor: (4 * M * eta * sqrt(pi/5)) / R
    # Note: This specific pre-factor is identical across several modes in the C source
    pre_factor: float = (4.0 * mass * Nu * np.sqrt(np.pi / 5.0)) / R

    # Evaluate the General Orbit (GO) and Quasi-Circular (QC) components
    go_component: complex = hGO_3_m_3(
        mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params
    )
    qc_component: complex = hQC_3_m_3(
        mass, Nu, vpnorder, x, S1z, S2z, params
    )

    # cpolar(1, -3 * Phi) -> exp(-i * 3 * Phi)
    # The higher 'm' index means the phase evolves 3 times faster than the orbit
    phase_factor: complex = np.exp(-3j * Phi)

    return pre_factor * (go_component + qc_component) * phase_factor

def hl_3_m_min3(
                mass: float, 
                Nu: float, 
                r: float, 
                rDOT: float, 
                Phi: float, 
                PhiDOT: float, 
                R: float, 
                vpnorder: int, 
                S1z: float, 
                S2z: float, 
                x: float, 
                params: KeplerVars
            ) -> complex:

    if vpnorder > 8:
        raise ValueError(
            "Error in hl_3_m_3: Input PN order parameter should be between [0, 8]."
        )

    # Pre-factor: (4 * M * eta * sqrt(pi/5)) / R
    # Note: This specific pre-factor is identical across several modes in the C source
    pre_factor: float = (4.0 * mass * Nu * np.sqrt(np.pi / 5.0)) / R

    # Evaluate the General Orbit (GO) and Quasi-Circular (QC) components
    # We assume hGO_3_m_3 and hQC_3_m_3 are defined elsewhere in your package
    go_component: complex = hGO_3_m_3(
        mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params
    )
    qc_component: complex = hQC_3_m_3(
        mass, Nu, vpnorder, x, S1z, S2z, params
    )
    # In C: conj(hGO + hQC)
    # In Python: use the .conjugate() method or np.conj()
    amplitude_term: complex = (go_component + qc_component).conjugate()

    # cpolar(1, -3 * Phi) -> exp(-i * 3 * Phi)
    # The higher 'm' index means the phase evolves 3 times faster than the orbit
    phase_factor: complex = np.exp(3j * Phi)

    return pre_factor * amplitude_term * phase_factor

# H32

def hGO_3_m_2(
                mass: float,
                Nu: float,
                r: float,
                rDOT: float,
                PhiDOT: float,
                vpnorder: int,
                S1z: float,
                S2z: float,
                x: float,
                params: KeplerVars,
            ) -> complex:
    delta = np.sqrt(1 - 4 * Nu)
    kappa1 = 1.0
    kappa2 = 1.0

    if vpnorder == 2:
        return (
            -(np.sqrt(0.7142857142857143) * mass * (-1 + 3 * Nu) * PhiDOT *
              (4 * PhiDOT * r + 1j * rDOT))
            / 6.0
        )

    elif vpnorder == 3:
        return (
            (np.sqrt(0.7142857142857143) * params.Mtot2 * Nu *
             (4 * PhiDOT * r + 1j * rDOT) * (S1z + S2z))
            / (3.0 * params.r2)
        )

    elif vpnorder == 4:
        return (
            -(mass * PhiDOT *
              (2 * mass *
               ((167 - 925 * Nu + 1615 * params.eta2) * PhiDOT * r +
                5j * (-82 + 239 * Nu + 55 * params.eta2) * rDOT) -
               3 * r *
               (2 * (-13 - 25 * Nu + 355 * params.eta2) * params.PhiDOT3 * params.r3 -
                60j * (-8 + 25 * Nu + params.eta2) * params.PhiDOT2 * params.r2 * rDOT +
                12 * (-23 + 70 * Nu + 10 * params.eta2) * PhiDOT * r * params.rDOT2 +
                5j * (-13 + 38 * Nu + 10 * params.eta2) * params.rDOT3)))
            / (108.0 * np.sqrt(35) * r)
        )

    elif vpnorder == 5:
        term1 = (
            (params.Mtot2 * Nu * PhiDOT *
             (7j * mass +
              r * (49j * params.PhiDOT2 * params.r2 +
                   90 * PhiDOT * r * rDOT - 6j * params.rDOT2)))
            / (4.0 * np.sqrt(35) * params.r2)
        )

        # Henry et al. ecc spin terms
        term2 = (
            (np.sqrt(0.7142857142857143) * params.Mtot2 *
             (2j * mass * rDOT *
              ((-12 + Nu * (97 + 4 * Nu) + delta * (-12 + 5 * Nu)) * S1z +
               (-12 + delta * (12 - 5 * Nu) + Nu * (97 + 4 * Nu)) * S2z) +
              4 * mass * PhiDOT * r *
              (-((12 + delta * (12 - 23 * Nu) + Nu * (53 + 8 * Nu)) * S1z) -
               (12 + Nu * (53 + 8 * Nu) + delta * (-12 + 23 * Nu)) * S2z) -
              3 * r *
              (16 * params.eta2 * PhiDOT * r *
               (4 * params.PhiDOT2 * params.r2 -
                2j * PhiDOT * r * rDOT + params.rDOT2) *
               (S1z + S2z) +
               30j * params.PhiDOT2 * params.r2 * rDOT *
               (S1z + delta * S1z + S2z - delta * S2z) +
               Nu * (-4 * params.PhiDOT3 * params.r3 *
                     ((5 + 17 * delta) * S1z + (5 - 17 * delta) * S2z) -
                     1j * params.PhiDOT2 * params.r2 * rDOT *
                     ((189 + 17 * delta) * S1z + (189 - 17 * delta) * S2z) +
                     20 * PhiDOT * r * params.rDOT2 *
                     (-((-3 + delta) * S1z) + (3 + delta) * S2z) -
                     4j * params.rDOT3 *
                     ((-4 + delta) * S1z - (4 + delta) * S2z)))))
            / (72.0 * params.r3)
        )

        return term1 + term2

    elif vpnorder == 6:
        term1 = (
            -(mass * PhiDOT *
              (4 * params.Mtot2 *
               (2 * (5377 + 6438 * Nu - 79866 * params.eta2 + 37348 * params.eta3) *
                PhiDOT * r -
                5j * (-4115 + 18399 * Nu - 20276 * params.eta2 + 7 * params.eta3) *
                rDOT) -
               4 * mass * r *
               ((4599 - 15737 * Nu + 36259 * params.eta2 + 108563 * params.eta3) *
                params.PhiDOT3 * params.r3 -
                1j * (-34053 + 59698 * Nu + 192949 * params.eta2 + 16193 * params.eta3) *
                params.PhiDOT2 * params.r2 * rDOT +
                (-59058 + 77983 * Nu + 322468 * params.eta2 - 4264 * params.eta3) *
                PhiDOT * r * params.rDOT2 +
                5j * (-3387 + 8518 * Nu + 8968 * params.eta2 + 884 * params.eta3) *
                params.rDOT3) +
               3 * params.r2 *
               (4 * (-710 + 3892 * Nu - 10655 * params.eta2 + 24000 * params.eta3) *
                params.PhiDOT5 * params.r5 +
                11j * (-1484 + 11693 * Nu - 25006 * params.eta2 + 428 * params.eta3) *
                params.PhiDOT4 * params.r4 * rDOT +
                4 * (4161 - 25618 * Nu + 29489 * params.eta2 + 22078 * params.eta3) *
                params.PhiDOT3 * params.r3 * params.rDOT2 +
                44j * (-151 + 1067 * Nu - 2419 * params.eta2 + 57 * params.eta3) *
                params.PhiDOT2 * params.r2 * params.rDOT3 +
                4 * (2041 - 11680 * Nu + 19334 * params.eta2 + 3368 * params.eta3) *
                PhiDOT * r * params.rDOT4 +
                5j * (477 - 2624 * Nu + 3862 * params.eta2 + 1160 * params.eta3) *
                params.rDOT5)))
            / (4752.0 * np.sqrt(35) * params.r2)
        )

        # Henry et al. ecc spin terms
        term2 = (
            (np.sqrt(0.7142857142857143) * params.Mtot3 *
             (2 * mass *
              (2j * (1 + delta - 2 * Nu) * S1z +
               Nu * ((1 + delta) * (-6 + kappa1) + 12 * Nu) * params.S1z2 +
               8 * Nu * (-1 + 3 * Nu) * S1z * S2z +
               S2z * (-2j * (-1 + delta + 2 * Nu) +
                      Nu * (-6 - delta * (-6 + kappa2) + kappa2 + 12 * Nu) * S2z)) +
              r * (4 * params.rDOT2 *
                   (2j * (1 + delta - 2 * Nu) * S1z +
                    Nu * ((1 + delta) * (-2 + kappa1) + 4 * Nu) * params.S1z2 +
                    8 * params.eta2 * S1z * S2z +
                    S2z * (-2j * (-1 + delta + 2 * Nu) +
                           Nu * (-2 - delta * (-2 + kappa2) + kappa2 + 4 * Nu) * S2z)) +
                   2 * params.PhiDOT2 * params.r2 *
                   (14j * (1 + delta - 2 * Nu) * S1z +
                    (6 * (1 + delta) * kappa1 -
                     (26 + 23 * kappa1 + delta * (26 + 11 * kappa1)) * Nu +
                     4 * (1 + 9 * kappa1) * params.eta2) * params.S1z2 -
                    64 * params.eta2 * S1z * S2z +
                    S2z * (-14j * (-1 + delta + 2 * Nu) +
                           2 * Nu * (-13 + 13 * delta + 2 * Nu) * S2z +
                           kappa2 * (6 + delta * (-6 + 11 * Nu) +
                                     Nu * (-23 + 36 * Nu)) * S2z)) +
                   PhiDOT * r * rDOT *
                   (40 * (1 + delta - 2 * Nu) * S1z -
                    1j * (8 * Nu * (-2 - 2 * delta + Nu) +
                          kappa1 * (3 + delta * (3 + 11 * Nu) +
                                    Nu * (5 + 18 * Nu))) * params.S1z2 +
                    20j * (-3 + Nu) * Nu * S1z * S2z +
                    S2z * (-40 * (-1 + delta + 2 * Nu) -
                           1j * (8 * Nu * (-2 + 2 * delta + Nu) +
                                 kappa2 * (3 - delta * (3 + 11 * Nu) +
                                           Nu * (5 + 18 * Nu))) * S2z)))))
            / (24.0 * params.r4)
        )

        return term1 + term2

    elif vpnorder == 7:
        return (
            -0.000014029180695847363 *
            (params.Mtot2 *
             (3 * params.r2 *
              (-120 * params.eta3 *
               (3565 * params.PhiDOT5 * params.r5 +
                2321j * params.PhiDOT4 * params.r4 * rDOT +
                8244 * params.PhiDOT3 * params.r3 * params.rDOT2 -
                869j * params.PhiDOT2 * params.r2 * params.rDOT3 -
                56 * PhiDOT * r * params.rDOT4 -
                120j * params.rDOT5) *
               (S1z + S2z) +
               2475 * params.PhiDOT2 * params.r2 *
               (6 * params.PhiDOT3 * params.r3 +
                77j * params.PhiDOT2 * params.r2 * rDOT -
                72 * PhiDOT * r * params.rDOT2 +
                6j * params.rDOT3) *
               (S1z + delta * S1z + S2z - delta * S2z) -
               3 * Nu *
               (22j * params.PhiDOT4 * params.r4 * rDOT *
                (36322j + 5 * (2993 + 3893 * delta) * S1z +
                 5 * (2993 - 3893 * delta) * S2z) -
                25j * params.rDOT5 *
                ((1053 + 443 * delta) * S1z + (1053 - 443 * delta) * S2z) +
                44j * params.PhiDOT2 * params.r2 * params.rDOT3 *
                (-5444j + 5 * (1424 + 849 * delta) * S1z +
                 7120 * S2z - 4245 * delta * S2z) -
                20 * PhiDOT * r * params.rDOT4 *
                (-1782j + (5963 + 2969 * delta) * S1z +
                 5963 * S2z - 2969 * delta * S2z) +
                4 * params.PhiDOT5 * params.r5 *
                (-86889j + 10 * (2063 + 225 * delta) * S1z +
                 20630 * S2z - 2250 * delta * S2z) +
                4 * params.PhiDOT3 * params.r3 * params.rDOT2 *
                (234861j + 40 * (-1824 + 97 * delta) * S1z -
                 40 * (1824 + 97 * delta) * S2z)) +
               2 * params.eta2 *
               (params.PhiDOT5 * params.r5 *
                (-1549757j + 300 * (1448 + 1311 * delta) * S1z +
                 300 * (1448 - 1311 * delta) * S2z) +
                11j * params.PhiDOT4 * params.r4 * rDOT *
                (329548j + 15 * (4113 + 1411 * delta) * S1z +
                 61695 * S2z - 21165 * delta * S2z) +
                22j * params.PhiDOT2 * params.r2 * params.rDOT3 *
                (-23971j + 15 * (3829 + 243 * delta) * S1z +
                 57435 * S2z - 3645 * delta * S2z) +
                150j * params.rDOT5 *
                ((-503 + 92 * delta) * S1z - (503 + 92 * delta) * S2z) +
                10 * PhiDOT * r * params.rDOT4 *
                (4565j + 6 * (-6327 + 991 * delta) * S1z -
                 6 * (6327 + 991 * delta) * S2z) +
                21 * params.PhiDOT3 * params.r3 * params.rDOT2 *
                (161403j + 10 * (-1897 + 1471 * delta) * S1z -
                 10 * (1897 + 1471 * delta) * S2z))) -
              6 * mass * r *
              (60 * params.eta3 *
               (2417 * params.PhiDOT3 * params.r3 +
                7258j * params.PhiDOT2 * params.r2 * rDOT -
                4381 * PhiDOT * r * params.rDOT2 +
                480j * params.rDOT3) *
               (S1z + S2z) -
               165 *
               (1161 * params.PhiDOT3 * params.r3 -
                536j * params.PhiDOT2 * params.r2 * rDOT -
                2412 * PhiDOT * r * params.rDOT2 -
                270j * params.rDOT3) *
               (S1z + delta * S1z + S2z - delta * S2z) +
               2 * Nu *
               (1j * params.PhiDOT2 * params.r2 * rDOT *
                (1015784j + 5 * (30849 + 88721 * delta) * S1z +
                 5 * (30849 - 88721 * delta) * S2z) +
                params.PhiDOT3 * params.r3 *
                (-173371j + 5 * (61569 + 10789 * delta) * S1z +
                 307845 * S2z - 53945 * delta * S2z) -
                5 * PhiDOT * r * params.rDOT2 *
                (-115368j + (177417 + 52307 * delta) * S1z +
                 177417 * S2z - 52307 * delta * S2z) +
                100j * params.rDOT3 *
                ((-1545 + 181 * delta) * S1z - (1545 + 181 * delta) * S2z)) +
               params.eta2 *
               (20 * PhiDOT * r * params.rDOT2 *
                (-11187j - 48074 * S1z + 11057 * delta * S1z -
                 48074 * S2z - 11057 * delta * S2z) +
                725j * params.rDOT3 *
                (-73 * S1z + 31 * delta * S1z - 73 * S2z - 31 * delta * S2z) +
                params.PhiDOT3 * params.r3 *
                (603141j - 543040 * S1z + 404620 * delta * S1z -
                 20 * (27152 + 20231 * delta) * S2z) +
                1j * params.PhiDOT2 * params.r2 * rDOT *
                (-1798104j - 648485 * S1z + 105755 * delta * S1z -
                 5 * (129697 + 21151 * delta) * S2z))) +
              10 * params.Mtot2 *
              (24 * params.eta3 *
               (6981 * PhiDOT * r + 1600j * rDOT) *
               (S1z + S2z) -
               66 * (2027 * PhiDOT * r + 380j * rDOT) *
               (S1z + delta * S1z + S2z - delta * S2z) +
               30j * Nu * rDOT *
               (297 * (1 + delta) * kappa1 * params.S1z3 +
                297 * (1 + delta) * kappa1 * params.S1z2 * S2z +
                S1z * (17261 - 1641 * delta -
                       297 * (-1 + delta) * kappa2 * params.S2z2) +
                S2z * (17261 + 1641 * delta -
                       297 * (-1 + delta) * kappa2 * params.S2z2)) +
               8 * Nu * PhiDOT * r *
               (-7315j -
                4455 * (1 + delta) * kappa1 * params.S1z3 +
                3 * (3881 - 7757 * delta) * S2z -
                4455 * (1 + delta) * kappa1 * params.S1z2 * S2z +
                4455 * (-1 + delta) * kappa2 * params.S2z3 +
                3 * S1z *
                (3881 + 7757 * delta +
                 1485 * (-1 + delta) * kappa2 * params.S2z2)) +
               3 * params.eta2 *
               (-5j * rDOT *
                (S1z * (18793 + 223 * delta + 1188 * kappa1 * params.S1z2) +
                 (18793 - 223 * delta + 1188 * (-2 + kappa1) * params.S1z2) * S2z +
                 1188 * (-2 + kappa2) * S1z * params.S2z2 +
                 1188 * kappa2 * params.S2z3) +
                4 * PhiDOT * r *
                (-4939j - 23359 * S1z - 5563 * delta * S1z +
                 5940 * kappa1 * params.S1z3 +
                 (-23359 + 5563 * delta +
                  5940 * (-2 + kappa1) * params.S1z2) * S2z +
                 5940 * (-2 + kappa2) * S1z * params.S2z2 +
                 5940 * kappa2 * params.S2z3)))))
            / (np.sqrt(35) * params.r4)
        )

    else:
        return 0.0 + 0.0j
    
def hQC_3_m_2(
                mass: float, 
                Nu: float, 
                vpnorder: int, 
                x: float, 
                S1z: float, 
                S2z: float, 
                params: KeplerVars
            ) -> complex:
    """
    Translates the hQC_3_m_2 mode from C to Python using NumPy.
    'params' contains pre-computed powers of x (x3p5, x4, x4p5) and eta (eta2).
    """
    # Using np.exp for real-valued constants is fine
    b0: float = 2.0 * mass / np.exp(0.5)
    EulerGamma: float = 0.5772156649015329
    
    # Pre-calculating the shared log term using np.log
    common_log_term: complex = (
        -10.0 + 6.0 * EulerGamma - 3.0j * np.pi + 
        np.log(4096.0) + 6.0 * np.log(b0) + 9.0 * np.log(x)
    )

    if vpnorder == 5:
        return (
            -0.4444444444444444j * np.sqrt(0.7142857142857143) * (-1.0 + 3.0 * Nu) * params.x3p5 * common_log_term
        )

    elif vpnorder == 6:
        return (
            0.8888888888888888j * np.sqrt(0.7142857142857143) * Nu * (S1z + S2z) * params.x4 * common_log_term
        )

    elif vpnorder == 7:
        numerator = (
            -0.024691358024691357j * (193.0 - 680.0 * Nu + 230.0 * params.eta2) * params.x4p5 * common_log_term
        )
        return numerator / np.sqrt(35.0)

    else:
        return 0.0 + 0.0j
    
def hl_3_m_2(
                mass: float, 
                Nu: float, 
                r: float, 
                rDOT: float, 
                Phi: float, 
                PhiDOT: float, 
                R: float, 
                vpnorder: int, 
                S1z: float, 
                S2z: float, 
                x: float, 
                params: KeplerVars
            ) -> complex:

    if vpnorder > 8:
        raise ValueError(
            "Error in hl_3_m_2: Input PN order parameter should be between [0, 8]."
        )

    # Pre-factor: (4 * M * Nu * sqrt(pi/5)) / R
    pre_factor: float = (4.0 * mass * Nu * np.sqrt(np.pi / 5.0)) / R

    # Evaluate the General Orbit (GO) and Quasi-Circular (QC) components
    # hGO_3_m_2 and hQC_3_m_2 are assumed to be defined in your package
    go_component: complex = hGO_3_m_2(
        mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params
    )
    qc_component: complex = hQC_3_m_2(
        mass, Nu, vpnorder, x, S1z, S2z, params
    )

    # cpolar(1, -2 * Phi) -> exp(-i * 2 * Phi)
    phase_factor: complex = np.exp(-2j * Phi)

    return pre_factor * (go_component + qc_component) * phase_factor

def hl_3_m_min2(
                mass: float, 
                Nu: float, 
                r: float, 
                rDOT: float, 
                Phi: float, 
                PhiDOT: float, 
                R: float, 
                vpnorder: int, 
                S1z: float, 
                S2z: float, 
                x: float, 
                params: KeplerVars
            ) -> complex:

    if vpnorder > 8:
        raise ValueError(
            "Error in hl_3_m_2: Input PN order parameter should be between [0, 8]."
        )

    # Pre-factor: (4 * M * Nu * sqrt(pi/5)) / R
    pre_factor: float = (4.0 * mass * Nu * np.sqrt(np.pi / 5.0)) / R

    # Evaluate the General Orbit (GO) and Quasi-Circular (QC) components
    # hGO_3_m_2 and hQC_3_m_2 are assumed to be defined in your package
    go_component: complex = hGO_3_m_2(
        mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params
    )
    qc_component: complex = hQC_3_m_2(
        mass, Nu, vpnorder, x, S1z, S2z, params
    )
    amplitude_term: complex = (go_component + qc_component).conjugate()

    # cpolar(1, -2 * Phi) -> exp(-i * 2 * Phi)
    phase_factor: complex = np.exp(2j * Phi)

    return pre_factor * amplitude_term * phase_factor

# H31

def hGO_3_m_1(
                mass: float,
                Nu: float,
                r: float,
                rDOT: float,
                PhiDOT: float,
                vpnorder: int,
                S1z: float,
                S2z: float,
                x: float,
                params: KeplerVars,
            ) -> complex:
    delta = np.sqrt(1 - 4 * Nu)
    kappa1 = 1.0
    kappa2 = 1.0
    r0 = 2 * mass / np.exp(0.5)

    combination_a  = PhiDOT * r + 1j * rDOT
    combination_a2 = combination_a ** 2
    combination_a3 = combination_a ** 3
    combination_a4 = combination_a ** 4
    combination_a5 = combination_a ** 5

    combination_b  = PhiDOT * r - 1j * rDOT
    combination_b2 = combination_b ** 2
    combination_b3 = combination_b ** 3
    combination_b4 = combination_b ** 4

    if vpnorder == 1:
        return (
            delta * (mass * (7j * PhiDOT * r - 12 * rDOT) -
                     6j * r * combination_b * combination_a2)
            / (6.0 * np.sqrt(14) * r)
        )

    elif vpnorder == 3:
        return (
            delta *
            (6j * (-5 + 19 * Nu) * params.r2 * combination_b2 * combination_a3 +
             2 * params.Mtot2 *
             (1j * (-101 + 43 * Nu) * PhiDOT * r +
              (109 - 86 * Nu) * rDOT) +
             3 * mass * r *
             (-4j * (-9 + 14 * Nu) * params.PhiDOT3 * params.r3 +
              6 * (2 + 9 * Nu) * params.PhiDOT2 * params.r2 * rDOT -
              1j * (33 + 62 * Nu) * PhiDOT * r * params.rDOT2 +
              4 * (8 + 17 * Nu) * params.rDOT3))
            / (36.0 * np.sqrt(14) * params.r2)
        )

    elif vpnorder == 4:
        return (
            0.041666666666666664j * params.Mtot2 *
            (4 * mass * (-1 + 5 * Nu) * ((1 + delta) * S1z + (-1 + delta) * S2z) -
             r * (2 * params.rDOT2 *
                  (-6 * (1 + delta) * S1z + 5 * (5 + 3 * delta) * Nu * S1z +
                   6 * S2z - 6 * delta * S2z +
                   5 * (-5 + 3 * delta) * Nu * S2z) +
                  params.PhiDOT2 * params.r2 *
                  ((24 + 24 * delta - 87 * Nu + 31 * delta * Nu) * S1z +
                   (-24 + 24 * delta + 87 * Nu + 31 * delta * Nu) * S2z) +
                  2j * PhiDOT * r * rDOT *
                  ((6 + 6 * delta - 31 * Nu + 35 * delta * Nu) * S1z +
                   (-6 + 6 * delta + 31 * Nu + 35 * delta * Nu) * S2z)))
            / (np.sqrt(14) * params.r3)
        )

    elif vpnorder == 5:
        term1 = (
            delta *
            (-18j * (183 - 1579 * Nu + 3387 * params.eta2) * params.r3 *
             combination_b3 * combination_a4 +
             2 * params.Mtot3 *
             (1j * (26473 - 27451 * Nu + 9921 * params.eta2) * PhiDOT * r -
              12 * (623 - 732 * Nu + 1913 * params.eta2) * rDOT) +
             2 * params.Mtot2 * r *
             (-1j * (-8641 - 59189 * Nu + 31959 * params.eta2) * params.PhiDOT3 * params.r3 +
              (-32635 - 29345 * Nu + 29541 * params.eta2) * params.PhiDOT2 * params.r2 * rDOT -
              44j * (-256 + 781 * Nu + 840 * params.eta2) * PhiDOT * r * params.rDOT2 +
              6 * (-756 + 8238 * Nu + 7357 * params.eta2) * params.rDOT3) +
             3 * mass * params.r2 *
             (2j * (-2479 - 4505 * Nu + 16785 * params.eta2) * params.PhiDOT5 * params.r5 +
              4 * (817 + 1220 * Nu - 7449 * params.eta2) * params.PhiDOT4 * params.r4 * rDOT +
              6j * (-1679 + 1469 * Nu + 12233 * params.eta2) * params.PhiDOT3 * params.r3 * params.rDOT2 -
              32 * (-460 + 421 * Nu + 2514 * params.eta2) * params.PhiDOT2 * params.r2 * params.rDOT3 +
              1j * (-9851 + 17954 * Nu + 40968 * params.eta2) * PhiDOT * r * params.rDOT4 -
              12 * (-771 + 1126 * Nu + 3616 * params.eta2) * params.rDOT5))
            / (4752.0 * np.sqrt(14) * params.r3)
        )

        # Henry et al. ecc spin terms
        term2 = (
            params.Mtot3 *
            (kappa1 *
             (-1j * (-13 - 13 * delta + 68 * Nu + 42 * delta * Nu) * PhiDOT * r -
              34 * (1 + delta) * rDOT + 4 * (26 + 9 * delta) * Nu * rDOT) *
             params.S1z2 +
             S2z * (12 * delta * Nu * (7j * PhiDOT * r - 6 * rDOT) * S1z +
                    kappa2 *
                    (-1j * (13 - 13 * delta - 68 * Nu + 42 * delta * Nu) * PhiDOT * r +
                     2 * (17 - 17 * delta - 52 * Nu + 18 * delta * Nu) * rDOT) *
                    S2z))
            / (24.0 * np.sqrt(14) * params.r3)
        )

        return term1 + term2

    elif vpnorder == 6:
        term1 = (
            delta * params.Mtot2 * Nu *
            (668 * params.Mtot2 -
             2 * mass * r *
             (727 * params.PhiDOT2 * params.r2 -
              99j * PhiDOT * r * rDOT + 452 * params.rDOT2) +
             params.r2 *
             (-499 * params.PhiDOT4 * params.r4 +
              1534j * params.PhiDOT3 * params.r3 * rDOT +
              3072 * params.PhiDOT2 * params.r2 * params.rDOT2 -
              680j * PhiDOT * r * params.rDOT3 +
              1000 * params.rDOT4))
            / (180.0 * np.sqrt(14) * params.r4)
        )

        # Henry et al. ecc spin terms
        term2 = (
            0.0023148148148148147j * params.Mtot2 *
            (2 * params.Mtot2 *
             ((252 * (1 + delta) - (1277 + 1279 * delta) * Nu +
               8 * (12 + 47 * delta) * params.eta2) * S1z +
              (252 * (-1 + delta) + (1277 - 1279 * delta) * Nu +
               8 * (-12 + 47 * delta) * params.eta2) * S2z) +
             3 * params.r2 *
             (2 * params.rDOT4 *
              ((30 + Nu * (-187 + 318 * Nu) +
                delta * (30 + Nu * (-101 + 122 * Nu))) * S1z +
               (-30 + (187 - 318 * Nu) * Nu +
                delta * (30 + Nu * (-101 + 122 * Nu))) * S2z) +
              2 * params.PhiDOT4 * params.r4 *
              ((90 - Nu * (28 + 579 * Nu) +
                delta * (90 + Nu * (-800 + 551 * Nu))) * S1z +
               (-90 + Nu * (28 + 579 * Nu) +
                delta * (90 + Nu * (-800 + 551 * Nu))) * S2z) +
              2j * PhiDOT * r * params.rDOT3 *
              ((186 - Nu * (745 + 354 * Nu) +
                delta * (186 + Nu * (-191 + 554 * Nu))) * S1z +
               (-186 + Nu * (745 + 354 * Nu) +
                delta * (186 + Nu * (-191 + 554 * Nu))) * S2z) +
              3 * params.PhiDOT2 * params.r2 * params.rDOT2 *
              ((32 + Nu * (-451 + 230 * Nu) +
                delta * (32 + Nu * (691 + 626 * Nu))) * S1z +
               (-32 + (451 - 230 * Nu) * Nu +
                delta * (32 + Nu * (691 + 626 * Nu))) * S2z) +
              2j * params.PhiDOT3 * params.r3 * rDOT *
              ((-12 + Nu * (-341 + 315 * Nu) +
                delta * (-12 + Nu * (-91 + 1213 * Nu))) * S1z +
               (12 + (341 - 315 * Nu) * Nu +
                delta * (-12 + Nu * (-91 + 1213 * Nu))) * S2z)) -
             2 * mass * r *
             (2 * params.PhiDOT2 * params.r2 *
              ((-312 * (1 + delta) + 2 * (827 - 923 * delta) * Nu +
                5 * (-201 + 131 * delta) * params.eta2) * S1z +
               (-312 * (-1 + delta) - 2 * (827 + 923 * delta) * Nu +
                5 * (201 + 131 * delta) * params.eta2) * S2z) +
              2 * params.rDOT2 *
              ((105 + Nu * (-541 + 924 * Nu) +
                delta * (105 + Nu * (-77 + 494 * Nu))) * S1z +
               (-105 + (541 - 924 * Nu) * Nu +
                delta * (105 + Nu * (-77 + 494 * Nu))) * S2z) +
              1j * PhiDOT * r * rDOT *
              ((1104 - 7 * Nu * (439 + 597 * Nu) +
                delta * (1104 + Nu * (-3071 + 3083 * Nu))) * S1z +
               (-1104 + 7 * Nu * (439 + 597 * Nu) +
                delta * (1104 + Nu * (-3071 + 3083 * Nu))) * S2z)))
            / (np.sqrt(14) * params.r4)
        )

        return term1 + term2

    elif vpnorder == 7:
        M_PI2 = np.pi ** 2

        spin_terms = (
            1513512 * params.Mtot3 *
            (2 * mass *
             (PhiDOT * r *
              (1j * (-97 + 631 * Nu) * S1z +
               5 * (8 + 16 * Nu * (-7 + 15 * Nu) +
                    3 * kappa1 * (-39 + Nu * (149 + 4 * Nu))) * params.S1z2 +
               S2z * (97j - 631j * Nu -
                      5 * (8 + 16 * Nu * (-7 + 15 * Nu) +
                           3 * kappa2 * (-39 + Nu * (149 + 4 * Nu))) * S2z)) +
              rDOT *
              (2 * (-18 + 83 * Nu) * S1z -
               5j * (-4 * (6 + Nu * (-25 + 7 * Nu)) +
                     kappa1 * (155 + Nu * (-373 + 164 * Nu))) * params.S1z2 +
               S2z * (36 - 166 * Nu +
                      5j * (-4 * (6 + Nu * (-25 + 7 * Nu)) +
                            kappa2 * (155 + Nu * (-373 + 164 * Nu))) * S2z))) +
             r * (2 * params.rDOT3 *
                  (S1z * (18 - 110 * Nu -
                          5j * (69 * kappa1 - 214 * kappa1 * Nu +
                                4 * (5 + 4 * kappa1) * params.eta2) * S1z) +
                   2 * (-9 + 55 * Nu) * S2z +
                   5j * (69 * kappa2 - 214 * kappa2 * Nu +
                         4 * (5 + 4 * kappa2) * params.eta2) * params.S2z2) +
                  params.PhiDOT3 * params.r3 *
                  (S1z * (255j - 1403j * Nu +
                          10 * (28 * (3 - 8 * Nu) * Nu +
                                kappa1 * (51 + 2 * Nu * (-91 + 118 * Nu))) * S1z) +
                   1j * (-255 + 1403 * Nu) * S2z -
                   10 * (28 * (3 - 8 * Nu) * Nu +
                         kappa2 * (51 + 2 * Nu * (-91 + 118 * Nu))) * params.S2z2) +
                  2 * params.PhiDOT2 * params.r2 * rDOT *
                  (S1z * (255 - 1079 * Nu +
                          5j * (2 * (60 - 361 * Nu) * Nu +
                                kappa1 * (6 + Nu * (-47 + 164 * Nu))) * S1z) +
                   (-255 + 1079 * Nu) * S2z -
                   5j * (2 * (60 - 361 * Nu) * Nu +
                         kappa2 * (6 + Nu * (-47 + 164 * Nu))) * params.S2z2) +
                  PhiDOT * r * params.rDOT2 *
                  (4j * (-114 + 781 * Nu) * S1z +
                   5 * (-213 * kappa1 - 72 * Nu + 278 * kappa1 * Nu +
                        8 * (7 + 44 * kappa1) * params.eta2) * params.S1z2 +
                   S2z * (456j + 1065 * kappa2 * S2z +
                          2 * Nu * (-1562j - 5 * (-36 + 139 * kappa2 +
                                                   4 * (7 + 44 * kappa2) * Nu) * S2z)))))
        )

        orbital_terms = (
            delta *
            (52920 * (-4083 + Nu * (58311 + Nu * (-269240 + 405617 * Nu))) *
             params.r4 * combination_b4 * combination_a5 +
             840 * params.Mtot2 * params.r2 *
             ((-2555489 + 7 * Nu * (820078 + Nu * (-6623390 + 4948497 * Nu))) *
              params.PhiDOT5 * params.r5 +
              1j * (3537631 + 7 * Nu * (-2817653 + Nu * (-7052042 + 4017147 * Nu))) *
              params.PhiDOT4 * params.r4 * rDOT +
              3 * (-1428997 + 7 * Nu * (-1230747 + Nu * (-237418 + 4061717 * Nu))) *
              params.PhiDOT3 * params.r3 * params.rDOT2 +
              1j * (-5153011 + 7 * Nu * (-2375327 + 9 * Nu * (218846 + 1640185 * Nu))) *
              params.PhiDOT2 * params.r2 * params.rDOT3 +
              (-7761899 + 7 * Nu * (2892563 + 5998602 * Nu + 7493619 * params.eta2)) *
              PhiDOT * r * params.rDOT4 +
              3j * (-2422057 + 7 * Nu * (501045 + Nu * (2033141 + 2771816 * Nu))) *
              params.rDOT5) -
             8820 * mass * params.r3 *
             (2 * (111737 + Nu * (-366573 + Nu * (-618923 + 2278593 * Nu))) *
              params.PhiDOT7 * params.r7 +
              2j * (101844 + Nu * (-273675 - 871630 * Nu + 2069774 * params.eta2)) *
              params.PhiDOT6 * params.r6 * rDOT +
              2 * (341322 + Nu * (-1429938 + Nu * (-1206083 + 7690681 * Nu))) *
              params.PhiDOT5 * params.r5 * params.rDOT2 +
              8j * (90241 + 2 * Nu * (-206022 + Nu * (-62113 + 1003558 * Nu))) *
              params.PhiDOT4 * params.r4 * params.rDOT3 +
              2 * (410547 + Nu * (-2269686 + Nu * (762091 + 8400052 * Nu))) *
              params.PhiDOT3 * params.r3 * params.rDOT4 +
              4j * (217935 + 2 * Nu * (-573699 + 5 * Nu * (18671 + 445748 * Nu))) *
              params.PhiDOT2 * params.r2 * params.rDOT5 +
              (333969 + 2 * Nu * (-981471 + 4 * Nu * (154039 + 750016 * Nu))) *
              PhiDOT * r * params.rDOT6 +
              24j * (13245 + 2 * Nu * (-37005 + Nu * (14251 + 130160 * Nu))) *
              params.rDOT7) +
             2 * params.Mtot4 *
             (-4178597424j * rDOT +
              84j * rDOT *
              (38468500 * params.eta3 +
               648648j * (S1z + S2z) -
               90090 * ((-24 + 155 * kappa1) * params.S1z2 +
                        (-24 + 155 * kappa2) * params.S2z2) -
               420 * params.eta2 *
               (-122855 +
                3003 * ((2 + 11 * kappa1) * params.S1z2 -
                        18 * S1z * S2z +
                        (2 + 11 * kappa2) * params.S2z2)) +
               3 * Nu *
               (-103100846 - 1846845 * M_PI2 -
                564564j * S2z +
                6006 * (S1z * (-94j + 5 * (-52 + 63 * kappa1) * S1z) -
                        20 * S1z * S2z +
                        5 * (-52 + 63 * kappa2) * params.S2z2))) +
              PhiDOT * r *
              (1176172480 * params.eta3 +
               8 * (74084729 -
                    189189 * S1z * (97j + 5 * (-8 + 117 * kappa1) * S1z) -
                    189189 * S2z * (97j + 5 * (-8 + 117 * kappa2) * S2z)) -
               176400 * params.eta2 *
               (11251 +
                429 * ((2 + 13 * kappa1) * params.S1z2 -
                       22 * S1z * S2z +
                       (2 + 13 * kappa2) * params.S2z2)) +
               147 * Nu *
               (-65012788 + 4485195 * M_PI2 +
                4499352j * S2z +
                10296 * (S1z * (437j + 15 * (-32 + 71 * kappa1) * S1z) -
                         3860 * S1z * S2z +
                         15 * (-32 + 71 * kappa2) * params.S2z2)))) -
             3 * params.Mtot3 * r *
             (-4j * params.rDOT3 *
              (601018232 - 1359334480 * params.eta3 -
               756756 * S1z * (6j + 115 * kappa1 * S1z) -
               756756 * S2z * (6j + 115 * kappa2 * S2z) +
               231 * Nu *
               (8490448 + 503685 * M_PI2 +
                80808j * S2z +
                2184 * (S1z * (37j + 190 * kappa1 * S1z) +
                        70 * S1z * S2z +
                        190 * kappa2 * params.S2z2)) +
               58800 * params.eta2 *
               (-62596 +
                429 * ((-1 + 5 * kappa1) * params.S1z2 -
                       12 * S1z * S2z +
                       (-1 + 5 * kappa2) * params.S2z2))) -
              14j * params.PhiDOT2 * params.r2 * rDOT *
              (-229522160 * params.eta3 +
               8 * (48303859 +
                    135135 * S1z * (-17j + 2 * kappa1 * S1z) +
                    135135 * S2z * (-17j + 2 * kappa2 * S2z)) +
               2520 * params.eta2 *
               (100913 +
                286 * ((-31 + 5 * kappa1) * params.S1z2 -
                       72 * S1z * S2z +
                       (-31 + 5 * kappa2) * params.S2z2)) +
               7 * Nu *
               (125038052 + 2374515 * M_PI2 +
                5858424j * S2z -
                10296 * (S1z * (-569j + 25 * (-24 + 7 * kappa1) * S1z) +
                         700 * S1z * S2z +
                         25 * (-24 + 7 * kappa2) * params.S2z2))) +
              4 * PhiDOT * r * params.rDOT2 *
              (-1095987374 + 1035895280 * params.eta3 +
               378378 * S1z * (152j + 355 * kappa1 * S1z) +
               378378 * S2z * (152j + 355 * kappa2 * S2z) -
               490 * params.eta2 *
               (-5802767 +
                5148 * ((2 + 23 * kappa1) * params.S1z2 -
                        42 * S1z * S2z +
                        (2 + 23 * kappa2) * params.S2z2)) -
               77 * Nu *
               (42451610 + 1511055 * M_PI2 +
                3623256j * S2z -
                6552 * (S1z * (-553j + 5 * (18 + 37 * kappa1) * S1z) +
                        965 * S1z * S2z +
                        5 * (18 + 37 * kappa2) * params.S2z2))) +
              7 * params.PhiDOT3 * params.r3 *
              (512893080 * params.eta3 -
               136 * (-2089567 +
                      135135 * S1z * (1j + 2 * kappa1 * S1z) +
                      135135 * S2z * (1j + 2 * kappa2 * S2z)) -
               560 * params.eta2 *
               (2457671 +
                2574 * ((11 + 53 * kappa1) * params.S1z2 -
                        84 * S1z * S2z +
                        (11 + 53 * kappa2) * params.S2z2)) +
               3 * Nu *
               (16621605 * M_PI2 +
                8 * (27468722 +
                     2681679j * S2z +
                     3003 * (S1z * (893j - 840 * S1z + 800 * kappa1 * S1z) -
                             3160 * S1z * S2z +
                             40 * (-21 + 20 * kappa2) * params.S2z2))))))
        )

        log_term = (
            74954880 * delta * params.Mtot3 *
            (mass * (-22j * PhiDOT * r - 24 * rDOT) +
             3 * r *
             (7j * params.PhiDOT3 * params.r3 +
              14 * params.PhiDOT2 * params.r2 * rDOT -
              8j * PhiDOT * r * params.rDOT2 +
              4 * params.rDOT3)) *
            np.log(r / r0)
        )

        return (
            1j * (spin_terms + orbital_terms + log_term)
            / (3.6324288e8 * np.sqrt(14) * params.r4)
        )

    else:
        return 0.0 + 0.0j
    
def hQC_3_m_1(
                mass: float, 
                Nu: float, 
                vpnorder: int, 
                x: float, 
                S1z: float, 
                S2z: float, 
                params: KeplerVars
            ) -> complex:
    
    delta: float = np.sqrt(1.0 - 4.0 * Nu)
    EulerGamma: float = 0.5772156649015329
    b0: float = 2.0 * mass / np.exp(0.5)
    r0: float = b0

    if vpnorder == 4:
        return (
            -0.005555555555555556 * (
                delta * params.x3 * (
                    -97.0 + 60.0 * EulerGamma - 30.0j * np.pi + 
                    60.0 * np.log(2.0) + 60.0 * np.log(b0) + 90.0 * np.log(x)
                )
            ) / np.sqrt(14.0)
        )

    elif vpnorder == 6:
        return (
            (delta * params.x4 * (
                -1552.0 + 960.0 * EulerGamma - 6487.0 * Nu + 420.0 * EulerGamma * Nu - 
                480.0j * np.pi - 210.0j * Nu * np.pi + 960.0 * np.log(2.0) + 
                420.0 * Nu * np.log(2.0) + 60.0 * (16.0 + 7.0 * Nu) * np.log(b0) + 
                90.0 * (16.0 + 7.0 * Nu) * np.log(x)
            )) / (1080.0 * np.sqrt(14.0))
        )

    elif vpnorder == 7:
        numerator_7 = (
            -9.44822373393802e-6 * (
                params.x4p5 * (
                    53132.0j * delta - 92232.0j * delta * EulerGamma + 35280.0j * delta * EulerGamma**2 - 
                    46116.0 * delta * np.pi + 35280.0 * delta * EulerGamma * np.pi - 2940.0j * delta * np.pi**2 + 
                    57036.0 * S1z + 57036.0 * delta * S1z - 35280.0 * EulerGamma * S1z - 35280.0 * delta * EulerGamma * S1z - 
                    114513.0 * Nu * S1z - 143031.0 * delta * Nu * S1z + 97020.0 * EulerGamma * Nu * S1z + 
                    114660.0 * delta * EulerGamma * Nu * S1z + 17640.0j * np.pi * S1z + 17640.0j * delta * np.pi * S1z - 
                    48510.0j * Nu * np.pi * S1z - 57330.0j * delta * Nu * np.pi * S1z - 57036.0 * S2z + 
                    57036.0 * delta * S2z + 35280.0 * EulerGamma * S2z - 35280.0 * delta * EulerGamma * S2z + 
                    114513.0 * Nu * S2z - 143031.0 * delta * Nu * S2z - 97020.0 * EulerGamma * Nu * S2z + 
                    114660.0 * delta * EulerGamma * Nu * S2z - 17640.0j * np.pi * S2z + 17640.0j * delta * np.pi * S2z + 
                    48510.0j * Nu * np.pi * S2z - 57330.0j * delta * Nu * np.pi * S2z - 92232.0j * delta * np.log(2.0) + 
                    70560.0j * delta * EulerGamma * np.log(2.0) + 35280.0 * delta * np.pi * np.log(2.0) - 
                    35280.0 * S1z * np.log(2.0) - 35280.0 * delta * S1z * np.log(2.0) + 97020.0 * Nu * S1z * np.log(2.0) + 
                    114660.0 * delta * Nu * S1z * np.log(2.0) + 35280.0 * S2z * np.log(2.0) - 35280.0 * delta * S2z * np.log(2.0) - 
                    97020.0 * Nu * S2z * np.log(2.0) + 114660.0 * delta * Nu * S2z * np.log(2.0) + 
                    35280.0j * delta * np.log(2.0)**2 - 2490264.0j * delta * np.log(3.0) + 
                    1905120.0j * delta * EulerGamma * np.log(3.0) + 952560.0 * delta * np.pi * np.log(3.0) - 
                    317520.0 * S1z * np.log(3.0) - 317520.0 * delta * S1z * np.log(3.0) + 1508220.0 * Nu * S1z * np.log(3.0) + 
                    396900.0 * delta * Nu * S1z * np.log(3.0) + 317520.0 * S2z * np.log(3.0) - 317520.0 * delta * S2z * np.log(3.0) - 
                    1508220.0 * Nu * S2z * np.log(3.0) + 396900.0 * delta * Nu * S2z * np.log(3.0) + 
                    1905120.0j * delta * np.log(2.0) * np.log(3.0) + 952560.0j * delta * np.log(3.0)**2 + 
                    35280.0j * delta * np.log(b0)**2 + 21840.0j * delta * np.log(r0) - 138348.0j * delta * np.log(x) + 
                    105840.0j * delta * EulerGamma * np.log(x) + 52920.0 * delta * np.pi * np.log(x) - 52920.0 * S1z * np.log(x) - 
                    52920.0 * delta * S1z * np.log(x) + 145530.0 * Nu * S1z * np.log(x) + 171990.0 * delta * Nu * S1z * np.log(x) + 
                    52920.0 * S2z * np.log(x) - 52920.0 * delta * S2z * np.log(x) - 145530.0 * Nu * S2z * np.log(x) + 
                    171990.0 * delta * Nu * S2z * np.log(x) + 105840.0j * delta * np.log(2.0) * np.log(x) + 
                    79380.0j * delta * np.log(x)**2 + 
                    588.0 * np.log(b0) * (
                        -194.0j * delta + 120.0j * delta * EulerGamma + 60.0 * delta * np.pi + 
                        15.0 * (-4.0 - 4.0 * delta + 11.0 * Nu + 13.0 * delta * Nu) * S1z + 
                        60.0 * S2z - 60.0 * delta * S2z - 165.0 * Nu * S2z + 195.0 * delta * Nu * S2z + 
                        120.0j * delta * np.log(2.0) + 180.0j * delta * np.log(x)
                    )
                )
            )
        )
        return numerator_7 / np.sqrt(14.0)

    else:
        return 0.0 + 0.0j

def hl_3_m_1(
                mass: float, 
                Nu: float, 
                r: float, 
                rDOT: float, 
                Phi: float, 
                PhiDOT: float, 
                R: float, 
                vpnorder: int, 
                S1z: float, 
                S2z: float, 
                x: float, 
                params: KeplerVars
            ) -> complex:

    if vpnorder > 8:
        # Maintaining the logic of the C-source error check
        raise ValueError(
            "Error in hl_3_m_1: Input PN order parameter should be between [0, 8]."
        )

    # Pre-factor: (4 * M * Nu * sqrt(pi/5)) / R
    pre_factor: float = (4.0 * mass * Nu * np.sqrt(np.pi / 5.0)) / R

    # Evaluate the General Orbit (GO) and Quasi-Circular (QC) components
    # Assuming hGO_3_m_1 and hQC_3_m_1 are already translated in your module
    go_component: complex = hGO_3_m_1(
        mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params
    )
    qc_component: complex = hQC_3_m_1(
        mass, Nu, vpnorder, x, S1z, S2z, params
    )

    # cpolar(1, -1 * Phi) -> exp(-i * Phi)
    phase_factor: complex = np.exp(-1j * Phi)

    return pre_factor * (go_component + qc_component) * phase_factor

def hl_3_m_min1(
                mass: float, 
                Nu: float, 
                r: float, 
                rDOT: float, 
                Phi: float, 
                PhiDOT: float, 
                R: float, 
                vpnorder: int, 
                S1z: float, 
                S2z: float, 
                x: float, 
                params: KeplerVars
            ) -> complex:

    if vpnorder > 8:
        # Maintaining the logic of the C-source error check
        raise ValueError(
            "Error in hl_3_m_1: Input PN order parameter should be between [0, 8]."
        )

    # Pre-factor: (4 * M * Nu * sqrt(pi/5)) / R
    pre_factor: float = (4.0 * mass * Nu * np.sqrt(np.pi / 5.0)) / R

    # Evaluate the General Orbit (GO) and Quasi-Circular (QC) components
    # Assuming hGO_3_m_1 and hQC_3_m_1 are already translated in your module
    go_component: complex = hGO_3_m_1(
        mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, x, params
    )
    qc_component: complex = hQC_3_m_1(
        mass, Nu, vpnorder, x, S1z, S2z, params
    )
    amplitude_factor = (go_component + qc_component).conjugate()

    # cpolar(1, 1 * Phi) -> exp(i * Phi)
    phase_factor: complex = np.exp(1j * Phi)

    return pre_factor * amplitude_factor * phase_factor

# H44

def hGO_4_m_4(
                mass: float,
                Nu: float,
                r: float,
                rDOT: float,
                PhiDOT: float,
                vpnorder: int,
                S1z: float,
                S2z: float,
                params: KeplerVars,
            ) -> complex:
    delta = np.sqrt((1 - 4 * Nu) + 0j)
    kappa1 = 1.0
    kappa2 = 1.0

    combination_a  = PhiDOT * r + 1j * rDOT
    combination_a4 = combination_a ** 4
    combination_a5 = combination_a ** 5
    combination_a6 = combination_a ** 6

    combination_b  = PhiDOT * r - 1j * rDOT
    combination_b2 = combination_b ** 2

    if vpnorder == 2:
        return (
            np.sqrt(0.7142857142857143) * (-1 + 3 * Nu) *
            (7 * params.Mtot2 + 6 * params.r2 * combination_a4 +
             3 * mass * r *
             (17 * params.PhiDOT2 * params.r2 +
              18j * PhiDOT * r * rDOT - 6 * params.rDOT2))
            / (36.0 * params.r2)
        )

    elif vpnorder == 4:
        return (
            (40 * params.Mtot3 * (314 - 987 * Nu + 195 * params.eta2) -
             60 * (23 - 159 * Nu + 291 * params.eta2) * params.r3 *
             combination_b * combination_a5 +
             params.Mtot2 * r *
             ((53143 - 199660 * Nu + 127500 * params.eta2) * params.PhiDOT2 * params.r2 +
              24j * (967 - 4615 * Nu + 5935 * params.eta2) * PhiDOT * r * rDOT -
              10 * (290 - 2033 * Nu + 4365 * params.eta2) * params.rDOT2) -
             3 * mass * params.r2 *
             ((613 - 920 * Nu + 6420 * params.eta2) * params.PhiDOT4 * params.r4 -
              8j * (-976 + 1745 * Nu + 3150 * params.eta2) * params.PhiDOT3 * params.r3 * rDOT +
              2 * (-6141 + 8980 * Nu + 31500 * params.eta2) * params.PhiDOT2 * params.r2 * params.rDOT2 +
              4j * (-1853 + 1730 * Nu + 13230 * params.eta2) * PhiDOT * r * params.rDOT3 -
              20 * (-83 + 30 * Nu + 762 * params.eta2) * params.rDOT4))
            / (1584.0 * np.sqrt(35) * params.r3)
        )

    elif vpnorder == 5:
        term1 = (
            params.Mtot2 * Nu *
            (6 * mass * (-43j * PhiDOT * r + 9 * rDOT) +
             r * (-734j * params.PhiDOT3 * params.r3 +
                  129 * params.PhiDOT2 * params.r2 * rDOT +
                  156j * PhiDOT * r * params.rDOT2 -
                  26 * params.rDOT3))
            / (24.0 * np.sqrt(35) * params.r3)
        )

        # Henry et al. ecc spin terms
        term2 = (
            params.Mtot2 *
            (-3j * params.PhiDOT2 * params.r3 * rDOT *
             ((-250 + 1221 * Nu - 1512 * params.eta2 +
               delta * (-250 + 849 * Nu)) * S1z +
              (-250 + delta * (250 - 849 * Nu) + 1221 * Nu -
               1512 * params.eta2) * S2z) -
             2 * mass * PhiDOT * r *
             ((-130 + 757 * Nu - 1224 * params.eta2 +
               delta * (-130 + 513 * Nu)) * S1z +
              (-130 + delta * (130 - 513 * Nu) + 757 * Nu -
               1224 * params.eta2) * S2z) -
             2j * mass * rDOT *
             ((-100 + 577 * Nu - 864 * params.eta2 +
               delta * (-100 + 333 * Nu)) * S1z +
              (-100 + delta * (100 - 333 * Nu) + 577 * Nu -
               864 * params.eta2) * S2z) -
             6 * params.PhiDOT3 * params.r4 *
             ((-65 + 263 * Nu - 291 * params.eta2 +
               delta * (-65 + 282 * Nu)) * S1z +
              (-65 + delta * (65 - 282 * Nu) + 263 * Nu -
               291 * params.eta2) * S2z) +
             12 * PhiDOT * params.r2 * params.rDOT2 *
             ((-40 + 201 * Nu - 252 * params.eta2 +
               delta * (-40 + 129 * Nu)) * S1z +
              (-40 + delta * (40 - 129 * Nu) + 201 * Nu -
               252 * params.eta2) * S2z) +
             6j * r * params.rDOT3 *
             ((-20 + 107 * Nu - 144 * params.eta2 +
               delta * (-20 + 63 * Nu)) * S1z +
              (-20 + delta * (20 - 63 * Nu) + 107 * Nu -
               144 * params.eta2) * S2z))
            / (72.0 * np.sqrt(35) * params.r3)
        )

        return term1 + term2

    elif vpnorder == 6:
        term1 = (
            (10 * params.Mtot4 *
             (-4477296 + 12734393 * Nu - 6895 * params.eta2 + 1043805 * params.eta3) +
             3150 * (-367 + 4337 * Nu - 17462 * params.eta2 + 23577 * params.eta3) *
             params.r4 * combination_b2 * combination_a6 +
             2 * params.Mtot3 * r *
             ((-36967579 + 245501977 * Nu - 459916170 * params.eta2 +
               150200680 * params.eta3) * params.PhiDOT2 * params.r2 +
              4j * (7571073 - 10780154 * Nu - 56898800 * params.eta2 +
                    43665510 * params.eta3) * PhiDOT * r * rDOT -
              10 * (1283609 - 5800627 * Nu + 3725295 * params.eta2 +
                    4771935 * params.eta3) * params.rDOT2) -
             params.Mtot2 * params.r2 *
             ((-28258134 + 3245207 * Nu + 144051250 * params.eta2 +
               136991820 * params.eta3) * params.PhiDOT4 * params.r4 -
              24j * (2371982 - 7733376 * Nu - 7948185 * params.eta2 +
                     9074870 * params.eta3) * params.PhiDOT3 * params.r3 * rDOT +
              7 * (6557973 - 50558069 * Nu + 59901380 * params.eta2 +
                   104752320 * params.eta3) * params.PhiDOT2 * params.r2 * params.rDOT2 +
              168j * (52044 - 1084807 * Nu + 1849450 * params.eta2 +
                      4171730 * params.eta3) * PhiDOT * r * params.rDOT3 -
              35 * (1083 - 1246819 * Nu + 2524240 * params.eta2 +
                    5995845 * params.eta3) * params.rDOT4) -
             105 * mass * params.r3 *
             ((116396 - 551405 * Nu + 560658 * params.eta2 +
               293036 * params.eta3) * params.PhiDOT6 * params.r6 +
              2j * (158192 - 670661 * Nu + 177718 * params.eta2 +
                    2163976 * params.eta3) * params.PhiDOT5 * params.r5 * rDOT +
              (-393665 + 1322392 * Nu + 1589680 * params.eta2 -
               8622660 * params.eta3) * params.PhiDOT4 * params.r4 * params.rDOT2 -
              8j * (-23048 + 209397 * Nu - 487057 * params.eta2 +
                    260396 * params.eta3) * params.PhiDOT3 * params.r3 * params.rDOT3 -
              (630647 - 3391000 * Nu + 2501958 * params.eta2 +
               7664096 * params.eta3) * params.PhiDOT2 * params.r2 * params.rDOT4 -
              2j * (218975 - 1037408 * Nu + 148970 * params.eta2 +
                    3699480 * params.eta3) * PhiDOT * r * params.rDOT5 +
              10 * (10233 - 44864 * Nu - 13050 * params.eta2 +
                    203280 * params.eta3) * params.rDOT6))
            / (1.44144e6 * np.sqrt(35) * params.r4)
        )

        # Henry et al. ecc spin terms
        term2 = (
            np.sqrt(0.7142857142857143) * params.Mtot3 * (-1 + 3 * Nu) *
            (12 * mass +
             r * (53 * params.PhiDOT2 * params.r2 +
                  26j * PhiDOT * r * rDOT - 8 * params.rDOT2)) *
            (kappa1 * (1 + delta - 2 * Nu) * params.S1z2 +
             S2z * (4 * Nu * S1z + kappa2 * S2z - delta * kappa2 * S2z -
                    2 * kappa2 * Nu * S2z))
            / (48.0 * params.r4)
        )

        return term1 + term2

    elif vpnorder == 7:
        # Henry et al. ecc+spin terms
        return (
            params.Mtot2 *
            (14 * params.Mtot2 *
             (120 * params.eta3 *
              (24635 * PhiDOT * r + 18657j * rDOT) *
              (S1z + S2z) -
              60 * (10039 * PhiDOT * r + 7706j * rDOT) *
              (S1z + delta * S1z + S2z - delta * S2z) +
              5 * Nu *
              (PhiDOT * r *
               (448616j + 703833 * S1z + 505635 * delta * S1z +
                703833 * S2z - 505635 * delta * S2z) +
               4j * rDOT *
               (4175j + 123114 * S1z + 76938 * delta * S1z +
                123114 * S2z - 76938 * delta * S2z)) -
              6 * params.eta2 *
              (2 * PhiDOT * r *
               (217374j + 549175 * S1z + 61075 * delta * S1z +
                549175 * S2z - 61075 * delta * S2z) +
               5j * rDOT *
               (9861j + 132241 * S1z + 18825 * delta * S1z +
                132241 * S2z - 18825 * delta * S2z))) -
             3 * params.r2 *
             (1680 * params.eta3 *
              (2833 * params.PhiDOT5 * params.r5 +
               18796j * params.PhiDOT4 * params.r4 * rDOT -
               13185 * params.PhiDOT3 * params.r3 * params.rDOT2 +
               1186j * params.PhiDOT2 * params.r2 * params.rDOT3 -
               5863 * PhiDOT * r * params.rDOT4 -
               2100j * params.rDOT5) *
              (S1z + S2z) -
              1050 *
              (454 * params.PhiDOT5 * params.r5 +
               1195j * params.PhiDOT4 * params.r4 * rDOT -
               1950 * params.PhiDOT3 * params.r3 * params.rDOT2 -
               442j * params.PhiDOT2 * params.r2 * params.rDOT3 -
               384 * PhiDOT * r * params.rDOT4 -
               184j * params.rDOT5) *
              (S1z + delta * S1z + S2z - delta * S2z) -
              6 * params.eta2 *
              (2 * params.PhiDOT5 * params.r5 *
               (-2459811j + 35 * (26517 + 10223 * delta) * S1z +
                (928095 - 357805 * delta) * S2z) -
               10j * params.rDOT5 *
               (35291j + 28 * (2183 + 1155 * delta) * S1z +
                (61124 - 32340 * delta) * S2z) +
               1j * params.PhiDOT2 * params.r2 * params.rDOT3 *
               (917901j + 1120 * (-191 + 1616 * delta) * S1z -
                1120 * (191 + 1616 * delta) * S2z) +
               5 * params.PhiDOT3 * params.r3 * params.rDOT2 *
               (85426j + 7 * (-148363 + 4365 * delta) * S1z -
                7 * (148363 + 4365 * delta) * S2z) -
               4 * PhiDOT * r * params.rDOT4 *
               (280067j + 70 * (5844 + 4817 * delta) * S1z -
                70 * (-5844 + 4817 * delta) * S2z) +
               1j * params.PhiDOT4 * params.r4 * rDOT *
               (10375501j + 70 * (65831 + 22871 * delta) * S1z -
                70 * (-65831 + 22871 * delta) * S2z)) +
              Nu * (-40j * params.rDOT5 *
                    (12203j + 42 * (843 + 656 * delta) * S1z -
                     42 * (-843 + 656 * delta) * S2z) +
                    4j * params.PhiDOT2 * params.r2 * params.rDOT3 *
                    (-266071j + 210 * (-2279 + 1010 * delta) * S1z -
                     210 * (2279 + 1010 * delta) * S2z) -
                    16 * PhiDOT * r * params.rDOT4 *
                    (58753j + 105 * (2030 + 1947 * delta) * S1z -
                     105 * (-2030 + 1947 * delta) * S2z) +
                    params.PhiDOT5 * params.r5 *
                    (-8997592j + 105 * (39959 + 15835 * delta) * S1z -
                     105 * (-39959 + 15835 * delta) * S2z) -
                    10 * params.PhiDOT3 * params.r3 * params.rDOT2 *
                    (254228j + 21 * (65293 + 34551 * delta) * S1z -
                     21 * (-65293 + 34551 * delta) * S2z) +
                    2j * params.PhiDOT4 * params.r4 * rDOT *
                    (12351083j + 105 * (43193 + 37913 * delta) * S1z -
                     105 * (-43193 + 37913 * delta) * S2z))) +
             mass * r *
             (2520 * params.eta3 *
              (8036 * params.PhiDOT3 * params.r3 +
               41814j * params.PhiDOT2 * params.r2 * rDOT -
               30537 * PhiDOT * r * params.rDOT2 -
               9064j * params.rDOT3) *
              (S1z + S2z) -
              210 *
              (11849 * params.PhiDOT3 * params.r3 +
               31868j * params.PhiDOT2 * params.r2 * rDOT +
               3572 * PhiDOT * r * params.rDOT2 +
               1508j * params.rDOT3) *
              (S1z + delta * S1z + S2z - delta * S2z) -
              12 * params.eta2 *
              (1j * params.PhiDOT2 * params.r2 * rDOT *
               (-1150397j + 35 * (154765 + 157449 * delta) * S1z +
                5416775 * S2z - 5510715 * delta * S2z) -
               PhiDOT * r * params.rDOT2 *
               (2306552j + 35 * (4931 + 110733 * delta) * S1z +
                172585 * S2z - 3875655 * delta * S2z) +
               params.PhiDOT3 * params.r3 *
               (12461121j + 35 * (18331 + 87381 * delta) * S1z +
                641585 * S2z - 3058335 * delta * S2z) -
               25j * params.rDOT3 *
               (36676j + 35 * (137 + 1053 * delta) * S1z +
                4795 * S2z - 36855 * delta * S2z)) +
              Nu * (-200j * params.rDOT3 *
                    (-2501j + 42 * (-308 + 51 * delta) * S1z -
                     42 * (308 + 51 * delta) * S2z) +
                    2 * PhiDOT * r * params.rDOT2 *
                    (-1951984j - 105 * (-37907 + 14661 * delta) * S1z +
                     105 * (37907 + 14661 * delta) * S2z) +
                    8j * params.PhiDOT2 * params.r2 * rDOT *
                    (-10005028j + 105 * (40991 + 31689 * delta) * S1z -
                     105 * (-40991 + 31689 * delta) * S2z) +
                    params.PhiDOT3 * params.r3 *
                    (88418488j + 105 * (57793 + 266391 * delta) * S1z -
                     105 * (-57793 + 266391 * delta) * S2z))))
            / (332640.0 * np.sqrt(35) * params.r4)
        )

    else:
        return 0.0 + 0.0j
    
def hQC_4_m_4(
                mass: float, 
                Nu: float, 
                vpnorder: int, 
                x: float, 
                S1z: float, 
                S2z: float, 
                params: KeplerVars
            ) -> complex:
    EulerGamma: float = 0.5772156649015329
    b0: float = 2.0 * mass / np.exp(0.5)

    if vpnorder == 5:
        # Pre-factor: Complex(0, 0.07407407407407407) / sqrt(35)
        return (
            (0.07407407407407407j * params.x3p5 * (
                1888.0 - 960.0 * EulerGamma - 5661.0 * Nu + 2880.0 * EulerGamma * Nu + 
                480.0j * np.pi - 1440.0j * Nu * np.pi - 2880.0 * np.log(2.0) + 
                8640.0 * Nu * np.log(2.0) + 960.0 * (-1.0 + 3.0 * Nu) * np.log(b0) + 
                1440.0 * (-1.0 + 3.0 * Nu) * np.log(x)
            )) / np.sqrt(35.0)
        )

    elif vpnorder == 7:
        # Pre-factor: Complex(0, 1.0020843354176687e-6) / sqrt(35)
        return (
            (1.0020843354176687e-6j * params.x4p5 * (
                -752360448.0 + 382556160.0 * EulerGamma + 2620364605.0 * Nu - 
                1333248000.0 * EulerGamma * Nu - 898312500.0 * params.eta2 + 
                458035200.0 * EulerGamma * params.eta2 - 191278080.0j * np.pi + 
                666624000.0j * Nu * np.pi - 229017600.0j * params.eta2 * np.pi + 
                1147668480.0 * np.log(2.0) - 3999744000.0 * Nu * np.log(2.0) + 
                1374105600.0 * params.eta2 * np.log(2.0) + 
                215040.0 * (1779.0 - 6200.0 * Nu + 2130.0 * params.eta2) * np.log(b0) + 
                322560.0 * (1779.0 - 6200.0 * Nu + 2130.0 * params.eta2) * np.log(x)
            )) / np.sqrt(35.0)
        )

    else:
        return 0.0 + 0.0j
    
def hl_4_m_4(
                mass: float, 
                Nu: float, 
                r: float, 
                rDOT: float, 
                Phi: float, 
                PhiDOT: float, 
                R: float, 
                vpnorder: int, 
                S1z: float, 
                S2z: float, 
                x: float, 
                params: KeplerVars
            ) -> complex:

    if vpnorder > 8:
        raise ValueError(
            "Error in hl_4_m_4: Input PN order parameter should be between [0, 8]."
        )

    # Pre-factor: (4 * M * Nu * sqrt(pi/5)) / R
    pre_factor: float = (4.0 * mass * Nu * np.sqrt(np.pi / 5.0)) / R

    # Evaluate the General Orbit (GO) and Quasi-Circular (QC) components
    # hGO_4_m_4 and hQC_4_m_4 are assuming your existing structure
    go_component: complex = hGO_4_m_4(
        mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, params
    )
    qc_component: complex = hQC_4_m_4(
        mass, Nu, vpnorder, x, S1z, S2z, params
    )

    # cpolar(1, -4 * Phi) -> exp(-i * 4 * Phi)
    phase_factor: complex = np.exp(-4j * Phi)

    return pre_factor * (go_component + qc_component) * phase_factor

def hl_4_m_min4(
                mass: float, 
                Nu: float, 
                r: float, 
                rDOT: float, 
                Phi: float, 
                PhiDOT: float, 
                R: float, 
                vpnorder: int, 
                S1z: float, 
                S2z: float, 
                x: float, 
                params: KeplerVars
            ) -> complex:

    if vpnorder > 8:
        raise ValueError(
            "Error in hl_4_m_4: Input PN order parameter should be between [0, 8]."
        )

    # Pre-factor: (4 * M * Nu * sqrt(pi/5)) / R
    pre_factor: float = (4.0 * mass * Nu * np.sqrt(np.pi / 5.0)) / R

    # Evaluate the General Orbit (GO) and Quasi-Circular (QC) components
    # hGO_4_m_4 and hQC_4_m_4 are assuming your existing structure
    go_component: complex = hGO_4_m_4(
        mass, Nu, r, rDOT, PhiDOT, vpnorder, S1z, S2z, params
    )
    qc_component: complex = hQC_4_m_4(
        mass, Nu, vpnorder, x, S1z, S2z, params
    )

    amplitude_factor = (go_component + qc_component).conjugate()

    # cpolar(1, 4 * Phi) -> exp(i * 4 * Phi)
    phase_factor: complex = np.exp(4j * Phi)

    return pre_factor * amplitude_factor * phase_factor

# H43

def hGO_4_m_3(
                mass: float,
                Nu: float,
                r: float,
                rDOT: float,
                PhiDOT: float,
                vpnorder: int,
                S1z: float,
                S2z: float,
                x: float,
                params: KeplerVars,
            ) -> complex:
    delta = np.sqrt(1 - 4 * Nu)
    kappa1 = 1.0
    kappa2 = 1.0

    if vpnorder == 3:
        return (
            0.16666666666666666j * delta * mass * (-1 + 2 * Nu) * PhiDOT *
            (4 * mass +
             r * (23 * params.PhiDOT2 * params.r2 +
                  10j * PhiDOT * r * rDOT - 2 * params.rDOT2))
            / (np.sqrt(70) * r)
        )

    elif vpnorder == 4:
        return (
            -0.041666666666666664j * np.sqrt(0.35714285714285715) *
            params.Mtot2 * Nu *
            (4 * mass +
             r * (23 * params.PhiDOT2 * params.r2 +
                  10j * PhiDOT * r * rDOT - 2 * params.rDOT2)) *
            ((-1 + delta) * S1z + S2z + delta * S2z)
            / params.r3
        )

    elif vpnorder == 5:
        return (
            0.0012626262626262627j * delta * mass * PhiDOT *
            (2 * params.Mtot2 * (972 - 2293 * Nu + 1398 * params.eta2) +
             2 * mass * r *
             ((1788 - 9077 * Nu + 13416 * params.eta2) * params.PhiDOT2 * params.r2 +
              3j * (-2796 + 5299 * Nu + 1622 * params.eta2) * PhiDOT * r * rDOT -
              2 * (-1200 + 2545 * Nu + 162 * params.eta2) * params.rDOT2) -
             3 * params.r2 *
             ((-524 - 489 * Nu + 6392 * params.eta2) * params.PhiDOT4 * params.r4 +
              4j * (796 - 1864 * Nu + 133 * params.eta2) * params.PhiDOT3 * params.r3 * rDOT +
              42 * (-51 + 94 * Nu + 56 * params.eta2) * params.PhiDOT2 * params.r2 * params.rDOT2 +
              4j * (-229 + 366 * Nu + 358 * params.eta2) * PhiDOT * r * params.rDOT3 -
              4 * (-43 + 62 * Nu + 80 * params.eta2) * params.rDOT4))
            / (np.sqrt(70) * params.r2)
        )

    elif vpnorder == 6:
        term1 = (
            delta * params.Mtot2 * Nu * PhiDOT *
            (6 * mass * (181 * PhiDOT * r - 89j * rDOT) +
             r * (4847 * params.PhiDOT3 * params.r3 -
                  7338j * params.PhiDOT2 * params.r2 * rDOT -
                  408 * PhiDOT * r * params.rDOT2 +
                  112j * params.rDOT3))
            / (180.0 * np.sqrt(70) * params.r2)
        )

        # Henry et al. ecc spin terms
        term2 = (
            -0.0006313131313131314j * params.Mtot2 *
            (2 * params.Mtot2 *
             ((-440 + 6801 * Nu - 1428 * params.eta2 +
               delta * (-440 - 3193 * Nu + 300 * params.eta2)) * S1z +
              (440 - 6801 * Nu + 1428 * params.eta2 +
               delta * (-440 - 3193 * Nu + 300 * params.eta2)) * S2z) -
             2 * mass * r *
             (-3j * PhiDOT * r * rDOT *
              (delta * (-1320 + 9093 * Nu + 59 * params.eta2) * S1z -
               5 * (264 - 311 * Nu + 823 * params.eta2) * S1z +
               delta * (-1320 + 9093 * Nu + 59 * params.eta2) * S2z +
               5 * (264 - 311 * Nu + 823 * params.eta2) * S2z) -
              2 * params.rDOT2 *
              ((220 + 1659 * Nu - 1512 * params.eta2 +
                delta * (220 - 3067 * Nu + 240 * params.eta2)) * S1z +
               (-220 - 1659 * Nu + 1512 * params.eta2 +
                delta * (220 - 3067 * Nu + 240 * params.eta2)) * S2z) +
              2 * params.PhiDOT2 * params.r2 *
              ((1826 - 19530 * Nu + 20145 * params.eta2 +
                delta * (1826 + 1534 * Nu + 567 * params.eta2)) * S1z +
               (-1826 + 19530 * Nu - 20145 * params.eta2 +
                delta * (1826 + 1534 * Nu + 567 * params.eta2)) * S2z)) -
             3 * params.r2 *
             (3080j * params.PhiDOT3 * params.r3 * rDOT *
              ((1 + delta) * S1z + (-1 + delta) * S2z) +
              2 * params.eta2 *
              (129 * params.PhiDOT2 * params.r2 * params.rDOT2 *
               (41 * S1z + 23 * delta * S1z - 41 * S2z + 23 * delta * S2z) -
               2 * params.rDOT4 *
               (149 * S1z + 75 * delta * S1z - 149 * S2z + 75 * delta * S2z) +
               2j * PhiDOT * r * params.rDOT3 *
               (925 * S1z + 491 * delta * S1z - 925 * S2z + 491 * delta * S2z) -
               1j * params.PhiDOT3 * params.r3 * rDOT *
               (-2105 * S1z + 753 * delta * S1z + 2105 * S2z + 753 * delta * S2z) +
               params.PhiDOT4 * params.r4 *
               (11617 * S1z + 4847 * delta * S1z - 11617 * S2z + 4847 * delta * S2z)) +
              Nu * (16 * params.PhiDOT4 * params.r4 *
                    (-413 * S1z + 127 * delta * S1z + 413 * S2z + 127 * delta * S2z) -
                    2 * params.rDOT4 *
                    (-301 * S1z + 213 * delta * S1z + 301 * S2z + 213 * delta * S2z) +
                    2j * PhiDOT * r * params.rDOT3 *
                    (-1625 * S1z + 1009 * delta * S1z + 1625 * S2z + 1009 * delta * S2z) +
                    3 * params.PhiDOT2 * params.r2 * params.rDOT2 *
                    (-2587 * S1z + 1267 * delta * S1z + 2587 * S2z + 1267 * delta * S2z) -
                    2j * params.PhiDOT3 * params.r3 * rDOT *
                    (3635 * S1z + 7981 * delta * S1z - 3635 * S2z + 7981 * delta * S2z))))
            / (np.sqrt(70) * params.r4)
        )

        return term1 + term2

    elif vpnorder == 7:
        # Henry et al. ecc+spin terms
        spin_terms = (
            5005 * params.Mtot2 *
            (-24 * PhiDOT * params.r2 * params.rDOT2 *
             (1j * (-11 + 48 * Nu) * S1z +
              6 * (-5 + 4 * kappa1) * params.eta2 * params.S1z2 +
              S2z * (11j - 48j * Nu -
                     6 * (-5 + 4 * kappa2) * params.eta2 * S2z)) +
             4 * r * params.rDOT3 *
             ((-11 + 93 * Nu) * S1z -
              6j * (-5 + 4 * kappa1) * params.eta2 * params.S1z2 +
              S2z * (11 - 93 * Nu +
                     6j * (-5 + 4 * kappa2) * params.eta2 * S2z)) +
             4 * mass * rDOT *
             ((22 - 111 * Nu) * S1z +
              6j * (15 + 8 * kappa1) * params.eta2 * params.S1z2 +
              S2z * (-22 + 111 * Nu -
                     6j * (15 + 8 * kappa2) * params.eta2 * S2z)) +
             6j * params.PhiDOT2 * params.r3 * rDOT *
             (1j * (-121 + 963 * Nu) * S1z +
              6 * (-55 * params.eta2 +
                   kappa1 * (-5 + 30 * Nu + 4 * params.eta2)) * params.S1z2 +
              S2z * (121j + 30 * kappa2 * S2z -
                     6 * (-55 + 4 * kappa2) * params.eta2 * S2z -
                     9 * Nu * (107j + 20 * kappa2 * S2z))) +
             2 * mass * PhiDOT * r *
             (-1j * (-121 + 633 * Nu) * S1z +
              6 * (180 * params.eta2 +
                   kappa1 * (9 - 54 * Nu + 28 * params.eta2)) * params.S1z2 -
              S2z * (121j + 54 * kappa2 * S2z +
                     24 * (45 + 7 * kappa2) * params.eta2 * S2z -
                     3 * Nu * (211j + 108 * kappa2 * S2z))) +
             params.PhiDOT3 * params.r4 *
             (-1j * (-649 + 4767 * Nu) * S1z +
              6 * (820 * params.eta2 +
                   kappa1 * (75 - 450 * Nu + 364 * params.eta2)) * params.S1z2 -
              S2z * (649j + 450 * kappa2 * S2z +
                     24 * (205 + 91 * kappa2) * params.eta2 * S2z -
                     3 * Nu * (1589j + 900 * kappa2 * S2z))))
        )

        orbital_terms = (
            delta *
            (2 * mass * PhiDOT * params.r3 *
             (10 * (234744 - 1010534 * Nu + 1024443 * params.eta2 +
                    5451096 * params.eta3) * params.PhiDOT4 * params.r4 -
              3j * (-2426804 + 1512854 * Nu + 4994115 * params.eta2 +
                    610960 * params.eta3) * params.PhiDOT3 * params.r3 * rDOT +
              (-30341028 + 23936528 * Nu + 89326545 * params.eta2 +
               19329660 * params.eta3) * params.PhiDOT2 * params.r2 * params.rDOT2 +
              21j * (-668008 + 803028 * Nu + 1908955 * params.eta2 +
                     540370 * params.eta3) * PhiDOT * r * params.rDOT3 -
              14 * (-172143 + 155683 * Nu + 680580 * params.eta2 +
                    111840 * params.eta3) * params.rDOT4) -
             105 * PhiDOT * params.r4 *
             ((-8280 + 24681 * Nu - 151973 * params.eta2 +
               624074 * params.eta3) * params.PhiDOT6 * params.r6 +
              2j * (-32208 + 248485 * Nu - 524074 * params.eta2 +
                    24546 * params.eta3) * params.PhiDOT5 * params.r5 * rDOT +
              2 * (48924 - 239802 * Nu + 137447 * params.eta2 +
                   358156 * params.eta3) * params.PhiDOT4 * params.r4 * params.rDOT2 +
              4j * (174 + 24488 * Nu - 102039 * params.eta2 +
                    44882 * params.eta3) * params.PhiDOT3 * params.r3 * params.rDOT3 +
              3 * (10455 - 56490 * Nu + 84504 * params.eta2 +
                   54016 * params.eta3) * params.PhiDOT2 * params.r2 * params.rDOT4 +
              2j * (11175 - 52698 * Nu + 57436 * params.eta2 +
                    60808 * params.eta3) * PhiDOT * r * params.rDOT5 -
              6 * (829 - 3726 * Nu + 3480 * params.eta2 +
                   4640 * params.eta3) * params.rDOT6) +
             params.Mtot2 * r *
             (20020 * params.rDOT3 * (S1z + S2z) *
              (-11 + 71 * Nu + 30j * params.eta2 * (S1z + S2z)) -
              52 * PhiDOT * r * params.rDOT2 *
              (129150 * params.eta3 +
               7 * Nu * (31961 + 8580j * S1z + 8580j * S2z) -
               3j * (32671j + 8470 * S1z + 8470 * S2z) -
               35 * params.eta2 *
               (33313 + 1980 * params.S1z2 + 3960 * S1z * S2z + 1980 * params.S2z2)) +
              params.PhiDOT3 * params.r3 *
              (-10566168 - 70869960 * params.eta3 +
               3248245j * S1z + 2252250 * kappa1 * params.S1z2 +
               3248245j * S2z + 2252250 * kappa2 * params.S2z2 -
               7 * Nu *
               (4818166 + 2480335j * S1z + 1287000 * kappa1 * params.S1z2 +
                2480335j * S2z + 1287000 * kappa2 * params.S2z2) +
               3850 * params.eta2 *
               (38873 + 78 * (7 + 30 * kappa1) * params.S1z2 -
                3588 * S1z * S2z +
                78 * (7 + 30 * kappa2) * params.S2z2)) -
              2j * params.PhiDOT2 * params.r2 * rDOT *
              (6531280 * params.eta3 +
               3 * (5382288 + 605605j * S1z + 150150 * kappa1 * params.S1z2 +
                    605605j * S2z + 150150 * kappa2 * params.S2z2) +
               210 * params.eta2 *
               (322144 + 2145 * (1 + 4 * kappa1) * params.S1z2 -
                12870 * S1z * S2z +
                2145 * (1 + 4 * kappa2) * params.S2z2) -
               21 * Nu *
               (2826484 + 85800 * kappa1 * params.S1z2 +
                515515j * S2z + 85800 * kappa2 * params.S2z2 -
                5005 * S1z * (-103j + 60 * S2z)))) -
             26 * params.Mtot3 *
             (770 * rDOT *
              (-90j * params.eta2 * params.S1z2 +
               S2z * (-22 + 67 * Nu - 90j * params.eta2 * S2z) +
               S1z * (-22 + Nu * (67 + 300j * S2z) - 180j * params.eta2 * S2z)) +
              PhiDOT * r *
              (-38076 + 174720 * params.eta3 -
               46585j * S1z - 20790 * kappa1 * params.S1z2 -
               46585j * S2z - 20790 * kappa2 * params.S2z2 -
               1260 * params.eta2 *
               (158 + 33 * (5 + 2 * kappa1) * params.S1z2 +
                198 * S1z * S2z +
                33 * (5 + 2 * kappa2) * params.S2z2) +
               7 * Nu *
               (6188 + 11880 * kappa1 * params.S1z2 +
                21505j * S2z + 11880 * kappa2 * params.S2z2 +
                55 * S1z * (391j + 744 * S2z)))))
        )

        return (
            -1.3875013875013875e-6j * mass *
            (spin_terms + orbital_terms)
            / (np.sqrt(70) * params.r4)
        )

    else:
        return Complex(0.0, 0.0)