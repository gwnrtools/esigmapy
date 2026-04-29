"""
esigma_pn_inspiral.py
====================
Python translation of ESIGMA_PNInspiral.c

References:
  - T. Damour, A. Gopakumar, and B. Iyer (PRD 70 064028), Eqs. 71a, 71b
  - K. G. Arun, Luc Blanchet, Bala R. Iyer and Siddharta Sinha, arXiv:0908.3854
  - I. Hinder, F. Herrmann, P. Laguna and D. Shoemaker, arXiv:0806.1037, Eqs. A35, A36
  - Huerta et al.

Notes
-----
* REAL8      -> float  (Python float is already double precision)
* M_PI       -> np.pi
* M_PI2      -> np.pi**2  (macro used in original)
* M_PI3      -> np.pi**3
* LAL_GAMMA  -> EULER_GAMMA constant defined below
* GSL_SUCCESS / GSL_FAILURE -> 0 / 1 (returned from ode function)
* XLAL_REAL8_FAIL_NAN -> float('nan')
* static functions -> module-level functions (no class needed)
* The ODE dispatcher (eccentric_x_model_odes) is translated to return
  (dxdt, dedt, dldt, dphidt) directly – the caller should use
  scipy.integrate.ode or solve_ivp with this as the rhs.
"""

import numpy as np

# ------ Constants --------------------------------------------------
EULER_GAMMA = 0.5772156649015329   # LAL_GAMMA / Euler–Mascheroni constant
LOG2        = 0.693147180559945309417232121458
LOG3        = 1.09861228866810969139524523692
LOG4        = 1.38629436111989061883446424292
LOG5        = 1.60943791243410037460075933323
REAL8_FAIL_NAN = float('nan')

# ------ Eccentric enhancement factors -------------------------------

def phi_e(e: float) -> float:
    e2  = e * e
    ef  = 1.0 - e2
    e4  = e2 * e2
    e6  = e4 * e2
    e8  = e6 * e2
    e10 = e8 * e2
    e12 = e10 * e2
    num = (1.0
           + (18970894028.0 * e2)  / 2649026657.0
           + (157473274.0  * e4)  / 30734301.0
           + (48176523.0   * e6)  / 177473701.0
           + (9293260.0    * e8)  / 3542508891.0
           - (5034498.0    * e10) / 7491716851.0
           + (428340.0     * e12) / 9958749469.0)
    den = ef**5
    return num / den


def psi_e(e: float) -> float:
    e2 = e * e
    ef = 1.0 - e2
    e4 = e2 * e2
    e6 = e4 * e2
    num = 1.0 - (185.0 * e2) / 21.0 - (3733.0 * e4) / 99.0 - (1423.0 * e6) / 104.0
    den = ef**6
    return num / den


def zed_e(e: float) -> float:
    e2 = e * e
    ef = 1.0 - e2
    e4 = e2 * e2
    e6 = e4 * e2
    num = 1.0 + (2095.0 * e2) / 143.0 + (1590.0 * e4) / 59.0 + (977.0 * e6) / 113.0
    den = ef**6
    return num / den


def kappa_e(e: float) -> float:
    e2  = e * e
    ef  = 1.0 - e2
    e4  = e2 * e2
    e6  = e4 * e2
    e8  = e6 * e2
    e10 = e8 * e2
    num = (1.0
           + (1497.0 * e2)  / 79.0
           + (7021.0 * e4)  / 143.0
           + (997.0  * e6)  / 98.0
           + (463.0  * e8)  / 51.0
           - (3829.0 * e10) / 120.0)
    den = ef**7
    return num / den


def phi_e_tilde(e: float) -> float:
    e2  = e * e
    ef  = 1.0 - e2
    e4  = e2 * e2
    e6  = e4 * e2
    e8  = e6 * e2
    e10 = e8 * e2
    num = (1.0
           + (413137256.0 * e2)  / 136292703.0
           + (37570495.0  * e4)  / 98143337.0
           - (2640201.0   * e6)  / 993226448.0
           - (4679700.0   * e8)  / 6316712563.0
           - (328675.0    * e10) / 8674876481.0)
    den = ef**3 * np.sqrt(ef)
    return num / den


def psi_e_tilde(e: float) -> float:
    e2  = e * e
    ef  = 1.0 - e2
    e4  = e2 * e2
    e6  = e4 * e2
    e8  = e6 * e2
    e10 = e8 * e2
    num = (1.0
           - (2022.0 * e2)  / 305.0
           - (249.0  * e4)  / 26.0
           - (193.0  * e6)  / 239.0
           + (23.0   * e8)  / 43.0
           - (102.0  * e10) / 463.0)
    den = ef**4 * np.sqrt(ef)
    return num / den


def zed_e_tilde(e: float) -> float:
    e2 = e * e
    ef = 1.0 - e2
    e4 = e2 * e2
    e6 = e4 * e2
    e8 = e6 * e2
    num = (1.0
           + (1563.0 * e2) / 194.0
           + (1142.0 * e4) / 193.0
           + (123.0  * e6) / 281.0
           - (27.0   * e8) / 328.0)
    den = ef**4 * np.sqrt(ef)
    return num / den


def kappa_e_tilde(e: float) -> float:
    e2  = e * e
    ef  = 1.0 - e2
    e4  = e2 * e2
    e6  = e4 * e2
    e8  = e6 * e2
    e10 = e8 * e2
    num = (1.0
           + (1789.0 * e2)  / 167.0
           + (5391.0 * e4)  / 340.0
           + (2150.0 * e6)  / 219.0
           - (1007.0 * e8)  / 320.0
           + (2588.0 * e10) / 189.0)
    den = ef**5
    return num / den


# ---------- Derived eccentricity functions ----------------------------------

def phi_e_rad(e: float) -> float:
    ef   = 1.0 - e * e
    e2   = e * e
    pre  = 192.0 * np.sqrt(ef) / (985.0 * e2)
    return pre * (np.sqrt(ef) * phi_e(e) - phi_e_tilde(e))


def psi_e_rad(e: float) -> float:
    ef   = 1.0 - e * e
    e2   = e * e
    pf1  = 18816.0 / (55691.0 * e2 * np.sqrt(ef))
    pf2  = 16382.0 * np.sqrt(ef) / (55691.0 * e2)
    b1   = np.sqrt(ef) * (1.0 - (11.0 / 7.0) * e2) * phi_e(e) - (1.0 - (3.0 / 7.0) * e2) * phi_e_tilde(e)
    b2   = np.sqrt(ef) * psi_e(e) - psi_e_tilde(e)
    return pf1 * b1 + pf2 * b2


def zed_e_rad(e: float) -> float:
    ef   = 1.0 - e * e
    e2   = e * e
    pf1  = 924.0 / (19067.0 * e2 * np.sqrt(ef))
    pf2  = 12243.0 * np.sqrt(ef) / (76268.0 * e2)
    b1   = -ef * np.sqrt(ef) * phi_e(e) + (1.0 - (5.0 / 11.0) * e2) * phi_e_tilde(e)
    b2   = np.sqrt(ef) * zed_e(e) - zed_e_tilde(e)
    return pf1 * b1 + pf2 * b2


def kappa_e_rad(e: float) -> float:
    ef   = 1.0 - e * e
    e2   = e * e
    den  = 769.0 / 96.0 - 3059665.0 * LOG2 / 700566.0 + 8190315.0 * LOG3 / 1868176.0
    pf   = np.sqrt(ef) / e2
    b1   = np.sqrt(ef) * kappa_e(e) - kappa_e_tilde(e)
    return pf * b1 / den


def f_e(e: float) -> float:
    ef   = 1.0 - e * e
    e2   = e * e
    e4   = e2 * e2
    e6   = e4 * e2
    e8   = e6 * e2
    denom = np.sqrt(ef) * ef**6
    num  = 1.0 + 85.0 * e2 / 6.0 + 5171.0 * e4 / 192.0 + 1751.0 * e6 / 192.0 + 297.0 * e8 / 1024.0
    return num / denom


def capital_f_e(e: float) -> float:
    ef   = 1.0 - e * e
    e2   = e * e
    e4   = e2 * e2
    e6   = e4 * e2
    denom = np.sqrt(ef) * ef**5
    num  = 1.0 + 2782.0 * e2 / 769.0 + 10721.0 * e4 / 6152.0 + 1719.0 * e6 / 24608.0
    return num / denom


def psi_n(e: float) -> float:
    ef   = 1.0 - e * e
    e2   = e * e
    pf1  = 1344.0 * (7.0 - 5.0 * e2) / (17599.0 * ef)
    pf2  = 8191.0 / 17599.0
    return pf1 * phi_e(e) + pf2 * psi_e(e)


def zed_n(e: float) -> float:
    pf1 = 583.0 / 567.0
    pf2 = 16.0  / 567.0
    return pf1 * zed_e(e) - pf2 * phi_e(e)


# ========= x_dot (dx/dt) terms ========================================

## --------- 0PN terms ------------

def x_dot_0pn(e: float, eta: float) -> float:
    """Eq. (A26)"""
    e2   = e * e
    e4   = e2 * e2
    ef   = 1.0 - e2
    num  = 2.0 * eta * (37.0 * e4 + 292.0 * e2 + 96.0)
    den  = 15.0 * ef**3 * np.sqrt(ef)
    e0   = 2.0 * eta * 96.0 / 15.0
    return (num / den) if e else e0

## --------- 1PN terms ------------

def x_dot_1pn(e: float, eta: float) -> float:
    """Eq. (A27)"""
    e2   = e * e
    e4   = e2 * e2
    e6   = e4 * e2
    ef   = 1.0 - e2
    t0   = 11888.0 + 14784.0 * eta
    t2   = e2  * (-87720.0 + 159600.0 * eta)
    t4   = e4  * (-171038.0 + 141708.0 * eta)
    t6   = e6  * (-11717.0  + 8288.0  * eta)
    den  = 420.0 * ef**4 * np.sqrt(ef)
    num  = -eta * (t0 + t2 + t4 + t6)
    e0   = -eta * (2972.0 / 105.0 + 176.0 * eta / 5.0)
    return (num / den) if e else e0

## --------- 1.5PN terms ------------

def x_dot_1_5_pn(e: float, eta: float, m1: float, m2: float,
                  S1z: float, S2z: float) -> float:
    """1.5PN spin-orbit term (Klein et al. arXiv:1801.08542, Eq. C1a)"""
    e2   = e * e
    e4   = e2 * e2
    e6   = e4 * e2
    ef   = 1.0 - e2
    M    = m1 + m2
    if e:
        return -(
            m1 * m2 * (
                (5424 + 27608*e2 + 16694*e4 + 585*e6) * m1*m1 * S1z +
                (5424 + 27608*e2 + 16694*e4 + 585*e6) * m2*m2 * S2z +
                3 * (1200 + 6976*e2 + 4886*e4 + 207*e6) * m1*m2 * (S1z + S2z)
            )
        ) / (45.0 * ef**5 * M**4)
    else:
        pre = 64.0 * eta / 5.0
        return pre * (-1.0/12.0 * (113*m1*m1*S1z + 113*m2*m2*S2z + 75*m1*m2*(S1z+S2z)) / M**2)


def x_dot_hereditary_1_5(e: float, eta: float, x: float) -> float:
    """Eq. (A28)"""
    pre = eta * x**6 * np.sqrt(x)
    if e:
        return pre * 256.0 * np.pi * phi_e(e) / 5.0
    else:
        return pre * 256.0 * np.pi / 5.0

## --------- 2PN terms ------------

def x_dot_2pn(e: float, eta: float) -> float:
    """Eq. (A29)"""

    e2 = e * e
    e4 = e2 * e2
    e6 = e4 * e2
    e8 = e4 * e4

    ef = 1.0 - e2
    eta2 = eta * eta

    # Small-e limit
    e0_lim = eta * (68206.0 / 2835.0 +
                    27322.0 * eta / 315.0 +
                    1888.0 * eta2 / 45.0)

    if e: return e0_lim
    else:
        sqrt_ef = np.sqrt(ef)

        den = 45360.0 * ef**5 * sqrt_ef

        e0_term = -360224.0 + 4514976.0 * eta + 1903104.0 * eta2
        e2_term = e2 * (-92846560.0 + 15464736.0 * eta + 61282032.0 * eta2)
        e4_term = e4 * (783768.0 - 207204264.0 * eta + 166506060.0 * eta2)
        e6_term = e6 * (83424402.0 - 123108426.0 * eta + 64828848.0 * eta2)
        e8_term = e8 * (3523113.0 - 3259980.0 * eta + 1964256.0 * eta2)

        rt_e0_term = 1451520.0 - 580608.0 * eta
        rt_e2_term = e2 * (64532160.0 - 25812864.0 * eta)
        rt_e4_term = e4 * (66316320.0 - 26526528.0 * eta)
        rt_e6_term = e6 * (2646000.0 - 1058400.0 * eta)

        num = eta * (
            e0_term + e2_term + e4_term + e6_term + e8_term
            + sqrt_ef * (rt_e0_term + rt_e2_term + rt_e4_term + rt_e6_term)
        )

        return num / den

def x_dot_2pn_SS(e: float, eta: float, m1: float, m2: float,
                 S1z: float, S2z: float) -> float:
    """2PN spin-spin term (Quentin Henry et al., arXiv:2308.13606)"""
    kappa1 = 1.0
    kappa2 = 1.0
    e2   = e * e
    e4   = e2 * e2
    e6   = e4 * e2
    ef   = 1.0 - e2
    ef_sqrt = np.sqrt(ef)
    ef_55 = ef**4 * ef * ef_sqrt
    M2   = (m1 + m2)**2
    M4   = M2**2
    pre  = 64.0 * eta / 5.0

    if e:
        return (
            m1 * m2 * (
                (48*(1+80*kappa1) + 3*e6*(9+236*kappa1) +
                 8*e2*(57+2692*kappa1) + 2*e4*(207+7472*kappa1)) * m1*m1 * S1z*S1z +
                2*(3792 + 21080*e2 + 14530*e4 + 681*e6) * m1*m2 * S1z*S2z +
                (48*(1+80*kappa2) + 3*e6*(9+236*kappa2) +
                 8*e2*(57+2692*kappa2) + 2*e4*(207+7472*kappa2)) * m2*m2 * S2z*S2z
            )
        ) / (60.0 * ef_55 * M4)
    else:
        return pre * (
            ((1 + 80*kappa1)*m1*m1*S1z*S1z +
             158*m1*m2*S1z*S2z +
             (1 + 80*kappa2)*m2*m2*S2z*S2z)
        ) / (16.0 * M2)

## --------- 2.5PN terms ------------

def x_dot_hereditary_2_5(e: float, eta: float, x: float) -> float:
    """See Huerta et al."""
    ef   = 1.0 - e*e
    ef2  = ef * ef
    e2   = e * e
    pre  = eta * x**7 * np.sqrt(x)

    if e:
        b1 = 256.0 * np.pi * phi_e(e) / ef
        b2 = (-17599.0 * np.pi * psi_n(e) / 35.0
              - 2268.0 * eta * np.pi * zed_n(e) / 5.0
              - 788.0 * e2 * np.pi * phi_e_rad(e) / ef2)
        return pre * (b1 + 2.0*b2/3.0)
    else:
        b1e0 = 256.0 * np.pi
        b2e0 = -17599.0 * np.pi / 35.0 - 2268.0 * eta * np.pi / 5.0
        return pre * (b1e0 + 2.0*b2e0/3.0)
    
def x_dot_2_5pn_SO(e: float, eta: float, m1: float, m2: float,
                    S1z: float, S2z: float) -> float:
    """2.5PN spin-orbit (Quentin Henry et al., arXiv:2308.13606)"""
    e2   = e * e
    e4   = e2 * e2
    e6   = e4 * e2
    e8   = e6 * e2
    ef   = 1.0 - e2
    ef_sqrt = np.sqrt(ef)
    M2   = (m1 + m2)**2
    M4   = M2 * M2
    M6   = M4 * M2
    pre  = 64.0 * eta / 5.0

    if e:
        return (
            -0.0000992063492063492 *
            m1 * m2 * (
                (4008832 + 808515*e8 +
                 896*e4*(126373 + 26748*ef_sqrt) +
                 16*e6*(2581907 + 30576*ef_sqrt) +
                 384*e2*(111403 + 61264*ef_sqrt)) * m1**4 * S1z +
                (4008832 + 808515*e8 +
                 896*e4*(126373 + 26748*ef_sqrt) +
                 16*e6*(2581907 + 30576*ef_sqrt) +
                 384*e2*(111403 + 61264*ef_sqrt)) * m2**4 * S2z +
                (1962112 + 1803645*e8 +
                 32*e6*(2542879 + 26754*ef_sqrt) +
                 896*e4*(202075 + 46809*ef_sqrt) +
                 768*e2*(51189 + 53606*ef_sqrt)) * m1*m1*m2*m2 * (S1z + S2z) +
                m1**3 * m2 * (
                    (3029504 + 1807155*e8 +
                     64*e6*(1314169 + 16562*ef_sqrt) +
                     128*e2*(340325 + 398216*ef_sqrt) +
                     112*e4*(1732273 + 463632*ef_sqrt)) * S1z +
                    3 * (414208 + 281775*e8 +
                         112*e6*(103639 + 728*ef_sqrt) +
                         336*e4*(77531 + 11888*ef_sqrt) +
                         128*e2*(55999 + 30632*ef_sqrt)) * S2z
                ) +
                m1 * m2**3 * (
                    3 * (414208 + 281775*e8 +
                         112*e6*(103639 + 728*ef_sqrt) +
                         336*e4*(77531 + 11888*ef_sqrt) +
                         128*e2*(55999 + 30632*ef_sqrt)) * S1z +
                    (3029504 + 1807155*e8 +
                     64*e6*(1314169 + 16562*ef_sqrt) +
                     128*e2*(340325 + 398216*ef_sqrt) +
                     112*e4*(1732273 + 463632*ef_sqrt)) * S2z
                )
            ) / ((-1 + e2)**6 * M6)
        )
    else:
        return pre * (
            -0.000992063492063492 *
            (31319*m1**4*S1z + 31319*m2**4*S2z +
             15329*m1*m1*m2*m2*(S1z + S2z) +
             4*m1**3*m2*(5917*S1z + 2427*S2z) +
             4*m1*m2**3*(2427*S1z + 5917*S2z)) / M4
        )

## --------- 3PN terms ------------

def x_dot_hereditary_3(e: float, eta: float, x: float) -> float:
    """See Huerta et al."""
    pi2  = np.pi * np.pi
    pre  = eta * x**8

    if e:
        x_3_term = (64.0 *
                    (-116761.0 * kappa_e(e) +
                     (19600.0*pi2 - 59920.0*EULER_GAMMA - 59920.0*LOG4 - 89880.0*np.log(x)) * f_e(e))
                    ) / 18375.0
        return pre * x_3_term
    else:
        x_3_term_e0 = (64.0 *
                  (-59920.0*EULER_GAMMA - 116761.0 + 19600.0*pi2 - 59920.0*LOG4 - 89880.0*np.log(x))
                  ) / 18375.0
        return pre * x_3_term_e0


def x_dot_3pn(e: float, eta: float, x: float) -> float:
    """3PN term (Huerta et al.)"""

    e2 = e * e
    e4 = e2 * e2
    e6 = e4 * e2
    e8 = e4 * e4
    e10 = e8 * e2

    ef = 1.0 - e2
    eta2 = eta * eta
    pi2 = np.pi * np.pi

    # Small-e limit
    e0_lim = eta * (
        426247111.0 / 222750.0
        - 56198689.0 * eta / 17010.0
        + 541.0 * eta2 / 70.0
        - 2242.0 * eta2 * eta / 81.0
        + 1804.0 * eta * pi2 / 15.0
        + 109568.0 * np.log(x) / 525.0
    )

    if abs(e)<1e-12: return e0_lim

    sqrt_ef = np.sqrt(ef)

    den = 598752000.0 * ef**7

    bit_1 = -25.0 * e10 * (
        81.0 * (99269280.0 - 33332681.0 * sqrt_ef)
        + 176.0 * eta * (
            9.0 * (-5132796.0 + 1874543.0 * sqrt_ef)
            + 4.0 * eta * (
                3582684.0 - 2962791.0 * sqrt_ef
                + 2320640.0 * sqrt_ef * eta
            )
        )
    )

    bit_2 = -128.0 * (
        3950984268.0 - 12902173599.0 * sqrt_ef
        + 275.0 * eta * (
            -1066392.0 + 57265081.0 * sqrt_ef
            + 81.0 * (-17696.0 + 16073.0 * sqrt_ef) * eta
            + 470820.0 * sqrt_ef * eta2
            - 46494.0 * (-1.0 + 45.0 * sqrt_ef) * pi2
        )
    )

    bit_3 = -32.0 * e2 * (
        -18.0 * (2603019496.0 + 19904214811.0 * sqrt_ef)
        + 55.0 * eta * (
            8147179440.0 - 5387647438.0 * sqrt_ef
            + 270.0 * (-6909392.0 + 9657701.0 * sqrt_ef) * eta
            + 901169500.0 * sqrt_ef * eta2
            - 193725.0 * (229.0 + 237.0 * sqrt_ef) * pi2
        )
    )

    bit_4 = -8.0 * e4 * (
        -6.0 * (312191560692.0 + 8654689873.0 * sqrt_ef)
        + 55.0 * eta * (
            42004763280.0 - 88628306866.0 * sqrt_ef
            - 1350.0 * (8601376.0 + 1306589.0 * sqrt_ef) * eta
            + 23638717900.0 * sqrt_ef * eta2
            + 891135.0 * (-2.0 + 627.0 * sqrt_ef) * pi2
        )
    )

    bit_5 = -2.0 * e8 * (
        4351589277552.0 - 1595548875627.0 * sqrt_ef
        + 550.0 * eta * (
            432.0 * (6368264.0 - 10627167.0 * sqrt_ef) * eta
            + 2201124800.0 * sqrt_ef * eta2
            + 9.0 * (
                8.0 * (-134041982.0 + 65136045.0 * sqrt_ef)
                + 861.0 * (14.0 + 891.0 * sqrt_ef) * pi2
            )
        )
    )

    bit_6 = -12.0 * e6 * (
        589550775792.0 - 6005081022.0 * sqrt_ef
        + 55.0 * eta * (
            90.0 * (90130656.0 - 311841025.0 * sqrt_ef) * eta
            + 17925404000.0 * sqrt_ef * eta2
            + 3.0 * (
                -2.0 * (5546517920.0 + 383583403.0 * sqrt_ef)
                + 4305.0 * (9046.0 + 19113.0 * sqrt_ef) * pi2
            )
        )
    )

    bit_7 = -40677120.0 * sqrt_ef * (
        3072.0 + 43520.0 * e2 + 82736.0 * e4 +
        28016.0 * e6 + 891.0 * e8
    ) * (
        LOG2 - np.log((1.0 / ef) + (1.0 / sqrt_ef)) - np.log(x)
    )

    num = eta * (bit_1 + bit_2 + bit_3 + bit_4 + bit_5 + bit_6 + bit_7)

    return num / den


def x_dot_3pn_SO(e: float, eta: float, m1: float, m2: float,
                  S1z: float, S2z: float) -> float:
    """3PN spin-orbit (Quentin Henry et al., arXiv:2308.13606)"""
    e2   = e * e
    e4   = e2 * e2
    e6   = e4 * e2
    e8   = e6 * e2
    ef   = 1.0 - e2
    ef_sqrt = np.sqrt(ef)
    ef_65 = ef**6 * ef_sqrt
    M2   = (m1 + m2)**2
    M4   = M2 * M2
    pre  = 64.0 * eta / 5.0

    if e:
        return (
            -9.645061728395062e-6 *
            m1 * m2 * np.pi * (
                (49766400 + 528887808*e2 + 814424832*e4 + 213166272*e6 + 3911917*e8) * m1*m1*S1z +
                (49766400 + 528887808*e2 + 814424832*e4 + 213166272*e6 + 3911917*e8) * m2*m2*S2z +
                3 * (11132928 + 133936128*e2 + 232455168*e4 + 67616624*e6 + 1479919*e8) * m1*m2*(S1z+S2z)
            ) / (ef_65 * M4)
        )
    else:
        return pre * (
            -1.0/6.0 * np.pi *
            (225*m1*m1*S1z + 225*m2*m2*S2z + 151*m1*m2*(S1z+S2z)) / M2
        )

def x_dot_3pn_SS(e: float, eta: float,
                 m1: float, m2: float,
                 S1z: float, S2z: float) -> float:
    """3PN spin-spin. Computed using inputs from energy and
        flux expressions. See Blanchet liv. rev. 
        + Bohe et al. arXiv: 1501.01529 + Marsat et al. arXiv: 1411.4118"""

    kappa1 = 1.0
    kappa2 = 1.0

    pre_factor = 64.0 * eta / 5.0

    e2 = e * e
    e4 = e2 * e2
    e6 = e4 * e2
    e8 = e6 * e2

    e_fact = 1.0 - e2
    e_fact_sqrt = np.sqrt(e_fact)

    e_fact_65 = (e_fact**6) * e_fact_sqrt

    M = m1 + m2
    M2 = M * M
    M4 = M2 * M2
    M6 = M4 * M2

    # --- Small-e fallback (better with tolerance) ---
    if abs(e) < 1e-12:
        return pre_factor * (
            (
                (36995 + 16358 * kappa1) * m1**4 * S1z * S1z +
                (36995 + 16358 * kappa2) * m2**4 * S2z * S2z +
                m1**3 * m2 * S1z *
                    ((34377 + 5864 * kappa1) * S1z + 59554 * S2z) +
                m1 * m2**3 * S2z *
                    (59554 * S1z + (34377 + 5864 * kappa2) * S2z) +
                3 * m1**2 * m2**2 *
                    ((1841 + 1318 * kappa1) * S1z * S1z +
                     35498 * S1z * S2z +
                     (1841 + 1318 * kappa2) * S2z * S2z)
            ) / (672.0 * M4)
        )

    # --- Main expression ---
    term1 = (
        27 * e8 * (15141 + 109852 * kappa1)
        + 8 * e4 * (
            22676605 + 3044160 * e_fact_sqrt
            + 32290722 * kappa1
            + 5327280 * e_fact_sqrt * kappa1
        )
        + 84 * e6 * (
            437079 + 448 * e_fact_sqrt
            + 14 * (89253 + 56 * e_fact_sqrt) * kappa1
        )
        - 576 * (
            -37891 + 896 * e_fact_sqrt
            + 2 * (-8963 + 784 * e_fact_sqrt) * kappa1
        )
        + 224 * e2 * (
            633619 + 107616 * e_fact_sqrt
            + 42 * (10029 + 4484 * e_fact_sqrt) * kappa1
        )
    )

    term2 = term1  # symmetric except kappa → kappa2
    # replace kappa1 with kappa2
    term2 = term2 + (kappa2 - kappa1) * (
        27 * e8 * 109852
        + 8 * e4 * (32290722 + 5327280 * e_fact_sqrt)
        + 84 * e6 * (14 * (89253 + 56 * e_fact_sqrt))
        - 576 * (2 * (-8963 + 784 * e_fact_sqrt))
        + 224 * e2 * (42 * (10029 + 4484 * e_fact_sqrt))
    )

    part1 = term1 * m1**4 * S1z * S1z
    part2 = term2 * m2**4 * S2z * S2z

    # Mixed terms (kept explicit for clarity)
    mix_common = (
        521613 * e8
        + 96 * (31457 - 1680 * e_fact_sqrt)
        + 294 * e6 * (72619 + 40 * e_fact_sqrt)
        + 112 * e2 * (247661 + 67260 * e_fact_sqrt)
        + e4 * (60534556 + 7610400 * e_fact_sqrt)
    )

    part3 = 6 * m1**3 * m2 * S1z * (
        (
            3 * e8 * (47103 + 202177 * kappa1)
            - 96 * (-34713 + 336 * e_fact_sqrt - 8440 * kappa1
                    + 2576 * e_fact_sqrt * kappa1)
            + 112 * e2 * (
                278677 + 13452 * e_fact_sqrt
                + 53216 * kappa1
                + 103132 * e_fact_sqrt * kappa1
            )
            + 14 * e6 * (
                772539 + 168 * e_fact_sqrt
                + 4 * (369781 + 322 * e_fact_sqrt) * kappa1
            )
            + 4 * e4 * (
                7 * (1655621 + 54360 * e_fact_sqrt)
                + 460 * (22397 + 6342 * e_fact_sqrt) * kappa1
            )
        ) * S1z
        + 2 * mix_common * S2z
    )

    part4 = 6 * m1 * m2**3 * S2z * (
        2 * mix_common * S1z +
        (
            3 * e8 * (47103 + 202177 * kappa2)
            - 96 * (-34713 + 336 * e_fact_sqrt - 8440 * kappa2
                    + 2576 * e_fact_sqrt * kappa2)
            + 112 * e2 * (
                278677 + 13452 * e_fact_sqrt
                + 53216 * kappa2
                + 103132 * e_fact_sqrt * kappa2
            )
            + 14 * e6 * (
                772539 + 168 * e_fact_sqrt
                + 4 * (369781 + 322 * e_fact_sqrt) * kappa2
            )
            + 4 * e4 * (
                7 * (1655621 + 54360 * e_fact_sqrt)
                + 460 * (22397 + 6342 * e_fact_sqrt) * kappa2
            )
        ) * S2z
    )

    part5 = m1**2 * m2**2 * (
        9 * (
            -86016 * e_fact_sqrt * kappa1
            + 192 * (1729 + 112 * e_fact_sqrt + 1766 * kappa1)
            + e8 * (58359 + 325902 * kappa1)
            + 28 * e6 * (
                121449 - 56 * e_fact_sqrt
                + (388890 + 224 * e_fact_sqrt) * kappa1
            )
            + 224 * e2 * (
                31367 - 4484 * e_fact_sqrt
                + 2 * (12189 + 8968 * e_fact_sqrt) * kappa1
            )
            + 8 * e4 * (
                1541617 - 126840 * e_fact_sqrt
                + 6 * (482187 + 84560 * e_fact_sqrt) * kappa1
            )
        ) * S1z * S1z
        + 8 * (
            8135280 + 1021401 * e8 - 467712 * e_fact_sqrt
            + 147 * e6 * (305971 + 232 * e_fact_sqrt)
            + 56 * e2 * (1084885 + 390108 * e_fact_sqrt)
            + 2 * e4 * (65163991 + 11035080 * e_fact_sqrt)
        ) * S1z * S2z
        + 9 * (
            -86016 * e_fact_sqrt * kappa2
            + 192 * (1729 + 112 * e_fact_sqrt + 1766 * kappa2)
            + e8 * (58359 + 325902 * kappa2)
            + 28 * e6 * (
                121449 - 56 * e_fact_sqrt
                + (388890 + 224 * e_fact_sqrt) * kappa2
            )
            + 224 * e2 * (
                31367 - 4484 * e_fact_sqrt
                + 2 * (12189 + 8968 * e_fact_sqrt) * kappa2
            )
            + 8 * e4 * (
                1541617 - 126840 * e_fact_sqrt
                + 6 * (482187 + 84560 * e_fact_sqrt) * kappa2
            )
        ) * S2z * S2z
    )

    numerator = m1 * m2 * (part1 + part2 + part3 + part4 + part5)
    denominator = 30240.0 * e_fact_65 * M6

    return numerator / denominator

## --------- 3.5PN terms ------------

def x_dot_3_5pnSO(e: float, eta: float, m1: float, m2: float,
                   S1z: float, S2z: float) -> float:
    """3.5PN spin-orbit"""
    pre = 64.0 * eta / 5.0
    M   = m1 + m2
    val_e0 = pre * (
        -0.00005511463844797178 *
        (3127800*m1**6*S1z + 3127800*m2**6*S2z +
         4914306*m1**3*m2**3*(S1z+S2z) +
         m1**5*m2*(6542338*S1z + 1195759*S2z) +
         m1**4*m2**2*(6694579*S1z + 3284422*S2z) +
         m1*m2**5*(1195759*S1z + 6542338*S2z) +
         m1**2*m2**4*(3284422*S1z + 6694579*S2z))
        / M**6
    )
    val_e = val_e0  #placeholder for e-dependent expression
    return val_e if e else val_e0


def x_dot_3_5pn_SS(e: float, eta: float, m1: float, m2: float,
                    S1z: float, S2z: float) -> float:
    """3.5PN spin-spin"""
    kappa1 = 1.0
    kappa2 = 1.0
    pre  = 64.0 * eta / 5.0
    M2   = (m1 + m2)**2
    val_e0 = pre * (
        np.pi * ((1 + 160*kappa1)*m1*m1*S1z*S1z +
              318*m1*m2*S1z*S2z +
              (1 + 160*kappa2)*m2*m2*S2z*S2z)
    ) / (8.0 * M2)
    val_e = val_e0  #placeholder for e-dependent expression
    return val_e if e else val_e0


def x_dot_3_5pn_cubicSpin(e: float, eta: float, m1: float, m2: float,
                            S1z: float, S2z: float) -> float:
    """3.5PN cubic-in-spin"""
    kappa1 = 1.0;  kappa2 = 1.0
    lambda1 = 1.0; lambda2 = 1.0
    pre  = 64.0 * eta / 5.0
    M    = m1 + m2
    val_e0 = pre * (
        -0.020833333333333332 *
        (2*(15+1016*kappa1+528*lambda1) * m1**4 * S1z**3 +
         2*(15+1016*kappa2+528*lambda2) * m2**4 * S2z**3 +
         m1**3*m2 * S1z*S1z * ((21+536*kappa1+1056*lambda1)*S1z + (4033+3712*kappa1)*S2z) +
         12*m1*m1*m2*m2 * S1z*S2z * ((89+434*kappa1)*S1z + (89+434*kappa2)*S2z) +
         m1*m2**3 * S2z*S2z * ((4033+3712*kappa2)*S1z + (21+536*kappa2+1056*lambda2)*S2z))
        / M**4
    )
    val_e = val_e0  #placeholder for e-dependent expression
    return val_e if e else val_e0


def x_dot_3_5_pn(e: float, eta: float) -> float:
    """3.5PN non-spinning"""
    val_e0 = (64.0 * eta * np.pi *
            (-4415.0/4032.0 + 358675.0*eta/6048.0 + 91495.0*eta*eta/1512.0)
            / 5.0)
    val_e = val_e0  #placeholder for e-dependent expression
    return val_e if e else val_e0

## --------- 4PN and 4.5PN terms ------------

def x_dot_4pn(e: float, eta: float, x: float) -> float:
    """4PN non-spinning (uses e=0 limit)"""
    euler = EULER_GAMMA
    pre   = 64.0 * eta / 5.0
    eta2  = eta * eta
    val_e0  = pre * (
        (3959271176713 - 20643291551545*eta) / 2.54270016e10 +
        eta2 * (2016887396 + 21*eta*(-1909807 + 49518*eta)) / 1.306368e6 -
        896*eta*euler/3.0 +
        (124741 + 734620*eta)*euler/4410.0 -
        361*np.pi**2/126.0 -
        eta*(1472377 + 928158*eta)*np.pi**2/16128.0 +
        127751*LOG2/1470.0 - 47385*LOG3/1568.0 +
        eta*(-850042*LOG2/2205.0 + 47385*LOG3/392.0) +
        (124741 - 582500*eta)*np.log(x)/8820.0
    )
    val_e = val_e0  #placeholder for e-dependent expression
    return val_e if e else val_e0

def x_dot_4pnSO(e: float, eta: float, m1: float, m2: float,
                 S1z: float, S2z: float) -> float:
    """4PN spin-orbit (uses e=0 limit)"""
    pre = 64.0 * eta / 5.0
    M   = m1 + m2
    val_e0 = pre * (
        -0.000496031746031746 *
        np.pi * (307708*m1**4*S1z + 307708*m2**4*S2z +
              93121*m1*m1*m2*m2*(S1z+S2z) +
              m1*m2**3*(119880*S1z + 75131*S2z) +
              m1**3*m2*(75131*S1z + 119880*S2z)) / M**4
    )
    val_e = val_e0  #placeholder for e-dependent expression
    return val_e if e else val_e0


def x_dot_4pnSS(e: float, eta: float, m1: float, m2: float,
                 S1z: float, S2z: float) -> float:
    """4PN spin-spin (uses e=0 limit)"""
    kappa1 = 1.0; kappa2 = 1.0
    pre = 64.0 * eta / 5.0
    M   = m1 + m2
    val_e0 = pre * (
        ((41400957 + 10676336*kappa1)*m1**6*S1z*S1z +
         (41400957 + 10676336*kappa2)*m2**6*S2z*S2z +
         m1**5*m2*S1z*(5*(16862889+4890568*kappa1)*S1z + 53648974*S2z) +
         m1*m2**5*S2z*(53648974*S1z + 5*(16862889+4890568*kappa2)*S2z) +
         m1**4*m2**2*((82037757+30833184*kappa1)*S1z*S1z + 168293278*S1z*S2z + (8950581+4778168*kappa2)*S2z*S2z) +
         m1**2*m2**4*((8950581+4778168*kappa1)*S1z*S1z + 168293278*S1z*S2z + 3*(27345919+10277728*kappa2)*S2z*S2z) +
         m1**3*m2**3*((43557309+18675776*kappa1)*S1z*S1z + 226571670*S1z*S2z + (43557309+18675776*kappa2)*S2z*S2z))
        / (72576.0 * M**6)
    )
    val_e = val_e0  #placeholder for e-dependent expression
    return val_e if e else val_e0

def x_dot_4_5_pn(e: float, eta: float, x: float) -> float:
    """4.5PN non-spinning (uses e=0 limit)"""
    euler = EULER_GAMMA
    pre   = 64.0 * eta / 5.0
    val_e0 = pre * (
        451*eta*np.pi**3/12.0 -
        np.pi * (700*eta*(3098001198 + eta*(525268513 + 289286988*eta)) +
              145786798080*euler +
              3*(-343801320119 + 97191198720*LOG2)) / 2.2353408e9 -
        3424*np.pi*np.log(x)/105.0
    )
    val_e = val_e0  #placeholder for e-dependent expression
    return val_e if e else val_e0


# Self-Force (SF) higher-order terms in addition to x_dot

def dxdt_4pn(x: float, eta: float) -> float:
    """4PN contribution to dx/dt"""

    x2 = x * x
    x3 = x2 * x
    x4 = x3 * x
    x5 = x4 * x

    eta2 = eta * eta
    eta3 = eta2 * eta
    eta4 = eta2 * eta2

    pi2 = np.pi * np.pi
    euler = EULER_GAMMA

    pre_fact = 64.0 * x5 * eta / 5.0

    bit_4pn = (
        x4 * (
            -(
                (-891.0 + 477.0 * eta - 11.0 * eta2)
                * (-4.928461199294532 + (9271.0 * eta) / 504.0 + (65.0 * eta2) / 18.0)
            ) / 36.0

            + (5.0 * (-3.7113095238095237 - (35.0 * eta) / 12.0)
               * (19683.0 - 66375.0 * eta + 1188.0 * eta2 + 19.0 * eta3
                  + 2214.0 * eta * pi2)
               ) / 648.0

            - (-3.0 - eta / 3.0) * (
                95.10839000836025
                - (134543.0 * eta) / 7776.0
                - (94403.0 * eta2) / 3024.0
                - (775.0 * eta3) / 324.0
                - (1712.0 * euler) / 105.0
                + (16.0 * pi2) / 3.0
                + (41.0 * eta * pi2) / 48.0
                - (3424.0 * LOG2) / 105.0
                - (856.0 * np.log(x)) / 105.0
            )

            + 2.0 * (
                -101.65745990813167
                + (232597.0 * euler) / 4410.0
                - (1369.0 * pi2) / 126.0
                + (39931.0 * LOG2) / 294.0
                - (47385.0 * LOG3) / 1568.0
                + (232597.0 * np.log(x)) / 8820.0
            )

            + 0.071630658436214 * (
                12774.514362657092
                - 39782.929590804066 * eta
                + 346.61235705608044 * eta2
                + 2.068222621184919 * eta3
                + eta4
                - 4169.536804308797 * eta * np.log(x)
            )
        )
    ) / 2.0

    return pre_fact * bit_4pn


##########################################################################
# dx_dt_4_5pn, dx_dt_5pn, dx_dt_5_5pn, dx_dt_6pn are not implemented here, 
# as they are not needed for ESIGMAHMv1. 
# (Implemented in the original C file, line number 2981 - 3468)
##########################################################################



# ================ e_dot (de/dt) terms ==================================

## ---------- 0PN -------------

def e_dot_0pn(e: float, eta: float) -> float:
    """Eq. (A31)"""
    e2  = e * e
    ef  = 1.0 - e2
    if not e:
        return 0.0
    num = -e * eta * (121.0*e2 + 304.0)
    den = 15.0 * ef**2 * np.sqrt(ef)
    return num / den

## ---------- 1PN -------------

def e_dot_1pn(e: float, eta: float) -> float:
    """Eq. (A32)"""
    if not e:
        return 0.0
    e2  = e * e
    e4  = e2 * e2
    ef  = 1.0 - e2
    t0  = 8.0 * (28588.0*eta + 8451.0)
    t2  = 12.0 * (54271.0*eta - 59834.0) * e2
    t4  = (93184.0*eta - 125361.0) * e4
    pre = e * eta / (2520.0 * ef**3 * np.sqrt(ef))
    return pre * (t0 + t2 + t4)

## ---------- 1.5PN -------------

def e_dot_1_5pn_SO(e: float, m1: float, m2: float,
                    S1z: float, S2z: float) -> float:
    """1.5PN SO eccentricity (Klein et al. arXiv:1801.08542, Eq. C1b)"""
    if not e:
        return 0.0
    e2  = e * e
    e4  = e2 * e2
    ef  = 1.0 - e2
    M   = m1 + m2
    pre = e / (ef**4) * (m1*m2) / (90.0 * M**4)
    return pre * (
        (19688 + 28256*e2 + 2367*e4)*m1*m1*S1z +
        (19688 + 28256*e2 + 2367*e4)*m2*m2*S2z +
        3*(4344 + 8090*e2 + 835*e4)*m1*m2*(S1z+S2z)
    )

def e_rad_hereditary_1_5(e: float, eta: float, x: float) -> float:
    if not e:
        return 0.0
    pre = 32.0 * eta * e * x**7 / 5.0
    return pre * (-985.0 * np.pi * phi_e_rad(e) / 48.0)

## ---------- 2PN -------------

def e_dot_2pn(e: float, eta: float) -> float:
    e_2_pn = 0.0

    eta_pow_2 = eta * eta
    e_pow_2 = e * e
    e_pow_4 = e_pow_2 * e_pow_2
    e_pow_6 = e_pow_4 * e_pow_2
    e_fact = 1.0 - e_pow_2

    zero_term = (
        -952397.0 / 1890.0
        + 5937.0 * eta / 14.0
        + 752.0 * eta_pow_2 / 5.0
    )

    e_2_term = e_pow_2 * (
        -3113989.0 / 2520.0
        - 388419.0 * eta / 280.0
        + 64433.0 * eta_pow_2 / 40.0
    )

    e_4_term = e_pow_4 * (
        4656611.0 / 3024.0
        - 13057267.0 * eta / 5040.0
        + 127411.0 * eta_pow_2 / 90.0
    )

    e_6_term = e_pow_6 * (
        420727.0 / 3360.0
        - 362071.0 * eta / 2520.0
        + 821.0 * eta_pow_2 / 9.0
    )

    zero_rt = 1336.0 / 3.0 - 2672.0 * eta / 15.0

    e_2_rt = e_pow_2 * (2321.0 / 2.0 - 2321.0 * eta / 5.0)

    e_4_rt = e_pow_4 * (565.0 / 6.0 - 113.0 * eta / 3.0)

    if e != 0.0:
        pre_factor = -e * eta / (e_fact**4 * np.sqrt(e_fact))
        e_2_pn = pre_factor * (
            zero_term
            + e_2_term
            + e_4_term
            + e_6_term
            + np.sqrt(e_fact) * (zero_rt + e_2_rt + e_4_rt)
        )
    else:
        e_2_pn = 0.0

    return e_2_pn


def e_dot_2pn_SS(e: float, m1: float, m2: float,
                  S1z: float, S2z: float) -> float:
    """2PN SS eccentricity (Quentin Henry et al., arXiv:2308.13606v1)"""
    if not e:
        return 0.0
    kappa1 = 1.0; kappa2 = 1.0
    e2  = e * e
    e4  = e2 * e2
    ef  = 1.0 - e2
    ef_45 = ef**4 * np.sqrt(ef)
    M2  = (m1+m2)**2
    M4  = M2*M2
    return (
        -0.008333333333333333 *
        e * m1 * m2 * (
            (45*(8+12*e2+e4) + 4*(3752+5950*e2+555*e4)*kappa1)*m1*m1*S1z*S1z +
            2*(14648+23260*e2+2175*e4)*m1*m2*S1z*S2z +
            (45*(8+12*e2+e4) + 4*(3752+5950*e2+555*e4)*kappa2)*m2*m2*S2z*S2z
        ) / (ef_45 * M4)
    )

## ---------- 2.5PN -------------

def e_dot_2_5pn_SO(e: float, m1: float, m2: float, S1z: float, S2z: float) -> float:
    if e == 0.0:
        return 0.0

    e_2 = e * e
    e_4 = e_2 * e_2
    e_6 = e_4 * e_2

    e_fact = 1.0 - e_2
    e_fact_2 = e_fact * e_fact
    e_fact_sqrt = np.sqrt(e_fact)
    e_fact_55 = e_fact_2 * e_fact_2 * e_fact * e_fact_sqrt

    M = m1 + m2
    M_fact_2 = M * M
    M_fact_4 = M_fact_2 * M_fact_2
    M_fact_6 = M_fact_4 * M_fact_2

    common_factor = e * m1 * m2

    term_A = (
        3 * (
            192 * (82880 + 36083 * e_fact_sqrt)
            + 168 * e_4 * (-211648 + 487675 * e_fact_sqrt)
            + 3 * e_6 * (-560896 + 1037433 * e_fact_sqrt)
            + 16 * e_2 * (1332912 + 6885449 * e_fact_sqrt)
        )
    )

    term_B = (
        3 * (
            27847680 - 9476416 * e_fact_sqrt
            + 7 * e_6 * (-420672 + 929363 * e_fact_sqrt)
            + 24 * e_4 * (-2592688 + 6508535 * e_fact_sqrt)
            + 16 * e_2 * (2332596 + 9517267 * e_fact_sqrt)
        )
    )

    term_C1 = (
        128 * (808080 - 259453 * e_fact_sqrt)
        + 3 * e_6 * (-3645824 + 6571731 * e_fact_sqrt)
        + 48 * e_2 * (2887976 + 10533923 * e_fact_sqrt)
        + 32 * e_4 * (-7222488 + 15232045 * e_fact_sqrt)
    )

    term_C2 = (
        9
        * (
            206464 * e_fact_sqrt
            + 21830032 * e_2 * e_fact_sqrt
            + 22413824 * e_4 * e_fact_sqrt
            + 1071519 * e_6 * e_fact_sqrt
            - 896 * (1 - e) * (1 + e) * (2960 + 6927 * e_2 + 313 * e_4)
        )
    )

    term_D1 = (
        9
        * (
            128 * (20720 + 1613 * e_fact_sqrt)
            + 256 * e_4 * (-23149 + 87554 * e_fact_sqrt)
            + 112 * e_2 * (31736 + 194911 * e_fact_sqrt)
            + e_6 * (-280448 + 1071519 * e_fact_sqrt)
        )
    )

    term_D2 = term_C1

    numerator = (
        common_factor
        * (
            term_A * (m1**4) * S1z
            + term_A * (m2**4) * S2z
            + term_B * (m1**2) * (m2**2) * (S1z + S2z)
            + (m1**3) * m2 * (term_C1 * S1z + term_C2 * S2z)
            + m1 * (m2**3) * (term_D1 * S1z + term_D2 * S2z)
        )
    )

    result = numerator / (60480.0 * e_fact_55 * M_fact_6)

    return result

def e_rad_hereditary_2_5(e: float, eta: float, x: float) -> float:
    if not e:
        return 0.0
    pre = 32.0 * eta * e * x**4 * x**4 * np.sqrt(x) / 5.0
    a2  = (55691.0 * psi_e_rad(e) / 1344.0 + 19067.0 * eta * zed_e_rad(e) / 126.0) * np.pi
    return pre * a2

## --------- 3PN -------------

def e_rad_hereditary_3(e: float, eta: float, x: float) -> float:
    if not e:
        return 0.0
    pre  = 32.0 * eta * e * x**7 / 5.0
    x32  = x * np.sqrt(x)
    a3   = (89789209.0/352800.0 - 87419.0*LOG2/630.0 + 78003.0*LOG3/560.0) * kappa_e_rad(e)
    a4   = (-769.0/96.0) * (16.0*np.pi**2/3.0 - 1712.0*EULER_GAMMA/105.0 - 1712.0*np.log(4.0*x32)/105.0) * capital_f_e(e)
    return pre * (a3 + a4)

def e_dot_3pn(e: float, eta: float, x: float) -> float:
    if e == 0.0:
        return 0.0

    eta_pow_2 = eta * eta
    eta_pow_3 = eta_pow_2 * eta

    e_pow_2 = e * e
    e_pow_4 = e_pow_2 * e_pow_2
    e_pow_6 = e_pow_4 * e_pow_2
    e_pow_8 = e_pow_6 * e_pow_2

    e_fact = 1.0 - e_pow_2
    sqrt_e_fact = np.sqrt(e_fact)

    pi_pow_2 =  np.np.pi *  np.np.pi

    zero_term = (
        7742634967.0 / 891000.0
        + (43386337.0 / 113400.0 + 1017.0 * pi_pow_2 / 10.0) * eta
        - 4148897.0 * eta_pow_2 / 2520.0
        - 61001.0 * eta_pow_3 / 486.0
    )

    e_2_term = e_pow_2 * (
        6556829759.0 / 891000.0
        + (770214901.0 / 25200.0 - 15727.0 * pi_pow_2 / 192.0) * eta
        - 80915371.0 * eta_pow_2 / 15120.0
        - 86910509.0 * eta_pow_3 / 19440.0
    )

    e_4_term = e_pow_4 * (
        -17072216761.0 / 2376000.0
        + (8799500893.0 / 907200.0 - 295559.0 * pi_pow_2 / 1920.0) * eta
        + 351962207.0 * eta_pow_2 / 20160.0
        - 2223241.0 * eta_pow_3 / 180.0
    )

    e_6_term = e_pow_6 * (
        17657772379.0 / 3696000.0
        + (-91818931.0 / 10080.0 - 6519.0 * pi_pow_2 / 640.0) * eta
        + 2495471.0 * eta_pow_2 / 252.0
        - 11792069.0 * eta_pow_3 / 2430.0
    )

    e_8_term = e_pow_8 * (
        302322169.0 / 1774080.0
        - 1921387.0 * eta / 10080.0
        + 41179.0 * eta_pow_2 / 216.0
        - 193396.0 * eta_pow_3 / 1215.0
    )

    zero_rt = (
        -22713049.0 / 15750.0
        + (-5526991.0 / 945.0 + 8323.0 * pi_pow_2 / 180.0) * eta
        + 54332.0 * eta_pow_2 / 45.0
    )

    e_2_rt = e_pow_2 * (
        89395687.0 / 7875.0
        + (-38295557.0 / 1260.0 + 94177.0 * pi_pow_2 / 960.0) * eta
        + 681989.0 * eta_pow_2 / 90.0
    )

    e_4_rt = e_pow_4 * (
        5321445613.0 / 378000.0
        + (-26478311.0 / 1512.0 + 2501.0 * pi_pow_2 / 2880.0) * eta
        + 225106.0 * eta_pow_2 / 45.0
    )

    e_6_rt = e_pow_6 * (
        186961.0 / 336.0
        - 289691.0 * eta / 504.0
        + 3197.0 * eta_pow_2 / 18.0
    )

    free_term = 730168.0 / (23625.0 * (1.0 + sqrt_e_fact))

    log_term = (
        (304.0 / 15.0)
        * (
            82283.0 / 1995.0
            + 297674.0 * e_pow_2 / 1995.0
            + 1147147.0 * e_pow_4 / 15960.0
            + 61311.0 * e_pow_6 / 21280.0
        )
        *  np.np.log(x * (1.0 + sqrt_e_fact) / (2.0 * e_fact))
    )

    pre_factor = -e * eta / (e_fact**5 * sqrt_e_fact)

    return pre_factor * (
        zero_term
        + e_2_term
        + e_4_term
        + e_6_term
        + e_8_term
        + sqrt_e_fact * (zero_rt + e_2_rt + e_4_rt + e_6_rt)
        + free_term
        + log_term
    )

def e_dot_3pn_SO(e: float, m1: float, m2: float,
                  S1z: float, S2z: float) -> float:
    """3PN SO eccentricity (Quentin Henry et al., arXiv:2308.13606v1)"""
    if not e:
        return 0.0
    e2  = e * e
    e4  = e2 * e2
    e6  = e4 * e2
    ef  = 1.0 - e2
    ef_sqrt = np.sqrt(ef)
    ef_55 = ef**4 * ef * ef_sqrt
    M2  = (m1+m2)**2
    M4  = M2*M2
    return (
        e * m1 * m2 * np.pi * (
            (64622592 + 238783104*e2 + 96887280*e4 + 2313613*e6)*m1*m1*S1z +
            (64622592 + 238783104*e2 + 96887280*e4 + 2313613*e6)*m2*m2*S2z +
            24*(1744704 + 8150400*e2 + 3941409*e4 + 122714*e6)*m1*m2*(S1z+S2z)
        )
    ) / (51840.0 * ef_55 * M4)


import math

def e_dot_3pn_SS(e: float, m1: float, m2: float, S1z: float, S2z: float) -> float:
    if e == 0.0:
        return 0.0

    kappa1 = 1.0
    kappa2 = 1.0

    e_2 = e * e
    e_4 = e_2 * e_2
    e_6 = e_4 * e_2

    e_fact = 1.0 - e_2
    e_fact_2 = e_fact * e_fact
    e_fact_sqrt = np.sqrt(e_fact)
    e_fact_55 = e_fact_2 * e_fact_2 * e_fact * e_fact_sqrt

    M = m1 + m2
    M_fact_2 = M * M
    M_fact_4 = M_fact_2 * M_fact_2
    M_fact_6 = M_fact_4 * M_fact_2

    prefactor_const = -0.000016534391534391536

    term = (
        e * m1 * m2
        * (
            # S1z^2 sector
            (
                27 * e_6 * (53837 + 368940 * kappa1)
                + 12 * e_4 * (
                    6350925 + 27328 * e_fact_sqrt + 16634634 * kappa1 + 47824 * e_fact_sqrt * kappa1
                )
                + 32 * (
                    2226847 + 545664 * e_fact_sqrt + 538758 * kappa1 + 954912 * e_fact_sqrt * kappa1
                )
                + 32 * e_2 * (
                    7292761 + 1157688 * e_fact_sqrt
                    + 9 * (847619 + 225106 * e_fact_sqrt) * kappa1
                )
            ) * m1**4 * S1z * S1z

            # S2z^2 sector
            + (
                27 * e_6 * (53837 + 368940 * kappa2)
                + 12 * e_4 * (
                    6350925 + 27328 * e_fact_sqrt + 16634634 * kappa2 + 47824 * e_fact_sqrt * kappa2
                )
                + 32 * (
                    2226847 + 545664 * e_fact_sqrt + 538758 * kappa2 + 954912 * e_fact_sqrt * kappa2
                )
                + 32 * e_2 * (
                    7292761 + 1157688 * e_fact_sqrt
                    + 9 * (847619 + 225106 * e_fact_sqrt) * kappa2
                )
            ) * m2**4 * S2z * S2z

            # m1^3 m2 S1z sector
            + 6 * m1**3 * m2 * S1z * (
                (
                    9 * e_6 * (55293 + 221219 * kappa1)
                    + 2 * e_4 * (
                        11036963 + 10248 * e_fact_sqrt + 19572200 * kappa1
                        + 78568 * e_fact_sqrt * kappa1
                    )
                    + 16 * (
                        798819 + 68208 * e_fact_sqrt - 336460 * kappa1
                        + 522928 * e_fact_sqrt * kappa1
                    )
                    + 8 * e_2 * (
                        7063231 + 289422 * e_fact_sqrt
                        + 6 * (698816 + 369817 * e_fact_sqrt) * kappa1
                    )
                ) * S1z

                + 2 * (
                    1818549 * e_6
                    + 240 * (31249 + 22736 * e_fact_sqrt)
                    + 2 * e_4 * (20863999 + 51240 * e_fact_sqrt)
                    + 8 * e_2 * (7764719 + 1447110 * e_fact_sqrt)
                ) * S2z
            )

            # m1 m2^3 S2z sector
            + 6 * m1 * m2**3 * S2z * (
                2 * (
                    1818549 * e_6
                    + 240 * (31249 + 22736 * e_fact_sqrt)
                    + 2 * e_4 * (20863999 + 51240 * e_fact_sqrt)
                    + 8 * e_2 * (7764719 + 1447110 * e_fact_sqrt)
                ) * S1z

                + (
                    9 * e_6 * (55293 + 221219 * kappa2)
                    + 2 * e_4 * (
                        11036963 + 10248 * e_fact_sqrt + 19572200 * kappa2
                        + 78568 * e_fact_sqrt * kappa2
                    )
                    + 16 * (
                        798819 + 68208 * e_fact_sqrt - 336460 * kappa2
                        + 522928 * e_fact_sqrt * kappa2
                    )
                    + 8 * e_2 * (
                        7063231 + 289422 * e_fact_sqrt
                        + 6 * (698816 + 369817 * e_fact_sqrt) * kappa2
                    )
                ) * S2z
            )

            # m1^2 m2^2 sector
            + m1**2 * m2**2 * (
                9 * (
                    3 * e_6 * (63301 + 361786 * kappa1)
                    + 4 * e_4 * (
                        1667883 - 3416 * e_fact_sqrt + 5133138 * kappa1
                        + 13664 * e_fact_sqrt * kappa1
                    )
                    + 32 * (
                        69895 - 22736 * e_fact_sqrt - 37046 * kappa1
                        + 90944 * e_fact_sqrt * kappa1
                    )
                    + 32 * e_2 * (
                        439124 - 48237 * e_fact_sqrt + 616472 * kappa1
                        + 192948 * e_fact_sqrt * kappa1
                    )
                ) * S1z * S1z

                + 8 * (
                    3553929 * e_6
                    + 3 * e_4 * (29583761 + 99064 * e_fact_sqrt)
                    + 116 * e_2 * (1176325 + 289422 * e_fact_sqrt)
                    + 8 * (2057131 + 1978032 * e_fact_sqrt)
                ) * S1z * S2z

                + 9 * (
                    3 * e_6 * (63301 + 361786 * kappa2)
                    + 4 * e_4 * (
                        1667883 - 3416 * e_fact_sqrt + 5133138 * kappa2
                        + 13664 * e_fact_sqrt * kappa2
                    )
                    + 32 * (
                        69895 - 22736 * e_fact_sqrt - 37046 * kappa2
                        + 90944 * e_fact_sqrt * kappa2
                    )
                    + 32 * e_2 * (
                        439124 - 48237 * e_fact_sqrt + 616472 * kappa2
                        + 192948 * e_fact_sqrt * kappa2
                    )
                ) * S2z * S2z
            )
        )
    )

    return prefactor_const * term / (e_fact_55 * M_fact_6)

## ---------- 3.5PN -------------

def e_dot_3_5pn(e: float, eta: float) -> float:
    """3.5PN eccentricity — zero."""
    return 0.0


# ================= l_dot (dl/dt) terms ===================================

## ---------- 1PN -------------

def l_dot_1pn(e: float, eta: float) -> float:
    """Eq. (A2)"""
    return 3.0 / (e*e - 1.0)

## ---------- 1.5PN -------------

def l_dot_1_5pn_SO(e: float, m1: float, m2: float,
                    S1z: float, S2z: float) -> float:
    """SO correction (Klein et al. arXiv:1801.08542, Eq. B1e/B2e)"""
    e2  = e * e
    ef  = 1.0 - e2
    pf  = 1.0 / (ef * np.sqrt(ef) * (m1+m2)**2)
    return pf * (4*m1*m1*S1z + 4*m2*m2*S2z + 3*m1*m2*(S1z+S2z))

## ---------- 2PN -------------

def l_dot_2pn_SS(e: float, m1: float, m2: float,
                  S1z: float, S2z: float) -> float:
    """SS correction (Klein et al. arXiv:1801.08542, Eq. B1e/B2e)"""
    kappa1 = 1.0; kappa2 = 1.0
    return (
        -3.0 * (m1*S1z*(2*m2*S2z + m1*S1z*kappa1) + m2*m2*S2z*S2z*kappa2)
    ) / (2.0 * (-1 + e**2)**2 * (m1+m2)**2)

def l_dot_2pn(e: float, eta: float) -> float:
    """Eq. (A3)"""
    ef  = 1.0 - e*e
    ef2 = ef*ef
    return ((26.0*eta - 51.0)*e*e + 28.0*eta - 18.0) / (4.0 * ef2)

## ---------- 2.5PN -------------

def l_dot_2_5pn_SO(e: float, m1: float, m2: float,
                    S1z: float, S2z: float) -> float:
    """2.5PN SO (Quentin Henry et al., arXiv:2308.13606v1)"""
    e2  = e * e
    ef  = 1.0 - e2
    ef2 = ef * ef
    ef_sqrt = np.sqrt(ef)
    ef_25 = ef2 * ef_sqrt
    M2  = (m1+m2)**2
    M4  = M2*M2
    return (
        (20*m1**4*(S1z + 3*e2*S1z) +
         20*(1+3*e2)*m2**4*S2z +
         (5+129*e2)*m1*m1*m2*m2*(S1z+S2z) +
         m1**3*m2*((23+137*e2)*S1z + 45*e2*S2z) +
         m1*m2**3*(23*S2z + e2*(45*S1z+137*S2z)))
    ) / (2.0 * ef_25 * M4)

## ---------- 3PN -------------

def l_dot_3pn(e: float, eta: float) -> float:
    """Eq. (A4)"""
    ef    = 1.0 - e*e
    ef3   = ef*ef*ef
    pre   = -1.0 / (128.0 * np.sqrt(ef) * ef3)
    eta2  = eta*eta
    pi2   = np.pi*np.pi
    e2    = e*e
    e4    = e2*e2
    t0    = 1920.0 - 768.0*eta
    t2    = (1920.0 - 768.0*eta)*e2
    t4    = (1536.0*eta - 3840.0)*e4
    r0    = 896.0*eta2 - 14624.0*eta + 492.0*pi2*eta - 192.0
    r2    = (5120.0*eta2 + 123.0*pi2*eta - 17856.0*eta + 8544.0)*e2
    r4    = (1040.0*eta2 - 1760.0*eta + 2496.0)*e4
    return pre * (t0 + t2 + t4 + np.sqrt(ef)*(r0 + r2 + r4))


def l_dot_3pn_SS(e: float, m1: float, m2: float,
                  S1z: float, S2z: float) -> float:
    """3PN SS (Quentin Henry et al., arXiv:2308.13606v1)"""
    kappa1 = 1.0; kappa2 = 1.0
    e2  = e * e
    M2  = (m1+m2)**2
    M4  = M2*M2
    return (
        (m1*m1*((-8+42*kappa1+6*e2*(4+13*kappa1))*m1*m1 +
                (-60+50*kappa1+11*e2*(3+11*kappa1))*m1*m2 +
                3*(-14+6*kappa1+3*e2*(1+8*kappa1))*m2*m2) * S1z*S1z +
         2*m1*m2*((6+93*e2)*m1*m1 + (12+157*e2)*m1*m2 + 3*(2+31*e2)*m2*m2)*S1z*S2z +
         m2*m2*(3*(-14+6*kappa2+3*e2*(1+8*kappa2))*m1*m1 +
                (-60+50*kappa2+11*e2*(3+11*kappa2))*m1*m2 +
                2*(-4+21*kappa2+3*e2*(4+13*kappa2))*m2*m2)*S2z*S2z)
    ) / (4.0 * (-1+e2)**3 * M4)


# =========== phi_dot (dphi/dt) terms =====================================================

def cosu_factor(e: float, u: float) -> float:
    return e * np.cos(u) - 1.0

## --------- 0PN -------------

def phi_dot_0pn(e: float, eta: float, u: float) -> float:
    """Eq. (A11)"""
    cf = cosu_factor(e, u)
    return np.sqrt(1.0 - e*e) / (cf * cf)

## --------- 1PN -------------

def phi_dot_1pn(e: float, eta: float, u: float) -> float:
    """Eq. (A12)"""
    cf = cosu_factor(e, u)
    return -(e*(eta - 4.0)*(e - np.cos(u))) / (np.sqrt(1.0 - e*e) * cf**3)

## --------- 1.5PN -------------

def phi_dot_1_5_pnSO_ecc(e: float, m1: float, m2: float,
                           S1z: float, S2z: float, u: float) -> float:
    """1.5PN SO phi_dot (Klein et al. arXiv:1801.08542, Eq. B1/B2)"""
    if not e:
        return 0.0
    cf = -1.0 + e * np.cos(u)
    return (2*e*(m1*S1z + m2*S2z)*(e - np.cos(u))) / ((-1+e*e)*(m1+m2)*cf**3)

## --------- 2PN -------------

def phi_dot_2pn(e: float, eta: float, u: float) -> float:
    """Eq. (A13)"""
    cf   = cosu_factor(e, u)
    cf5  = cf**5
    e2   = e*e; e3=e2*e; e4=e2*e2; e5=e4*e; e6=e4*e2
    ef   = 1.0 - e2
    eta2 = eta*eta
    cu   = np.cos(u); cu2=cu*cu; cu3=cu2*cu

    t0  = 90.0 - 36.0*eta
    t2  = (-2.0*eta2 + 50.0*eta + 75.0)*e2
    t4  = (20.0*eta2 - 26.0*eta - 60.0)*e4
    t6  = (-12.0*eta - 18.0)*eta*e6
    c1  = ((-eta2+97.0*eta+12.0)*e5 + (-16.0*eta2-74.0*eta-81.0)*e3 + (-eta2+67.0*eta-246.0)*e)*cu
    c2  = ((17.0*eta2-17.0*eta+48.0)*e6 + (-4.0*eta2-38.0*eta+153.0)*e4 + (5.0*eta2-35.0*eta+114.0)*e2)*cu2
    c3  = ((-14.0*eta2+8.0*eta-147.0)*e5 + (8.0*eta2+22.0*eta+42.0)*e3)*cu3

    r0  = (180.0-72.0*eta)*e2 + 36.0*eta - 90.0
    rc1 = ((144.0*eta-360.0)*e3 + (90.0-36.0*eta)*e)*cu
    rc2 = ((180.0-72.0*eta)*e4 + (90.0-36.0*eta)*e2)*cu2
    rc3 = e3*(36.0*eta-90.0)*cu3

    pre = 1.0 / (12.0 * np.sqrt(ef) * ef * cf5)
    return pre * (t0+t2+t4+t6+c1+c2+c3 + np.sqrt(ef)*(r0+rc1+rc2+rc3))


def phi_dot_2_pnSS_ecc(e: float, m1: float, m2: float,
                        S1z: float, S2z: float, u: float) -> float:
    """2PN SS phi_dot (Klein et al. arXiv:1801.08542, Eq. B1/B2)"""
    if not e: 
        return 0.0
    kappa1 = 1.0; kappa2 = 1.0
    return (
        e * (kappa2*m2*m2*S2z*S2z + m1*S1z*(kappa1*m1*S1z + 2*m2*S2z)) *
        (e - np.cos(u))
    ) / ((1-e**2)**1.5 * (m1+m2)**2 * (-1+e*np.cos(u))**3)

## --------- 2.5PN -------------

def phi_dot_2_5pn_SO(e: float, m1: float, m2: float, S1z: float, S2z: float, u: float) -> float:
    e_2 = e * e
    e_3 = e_2 * e
    e_4 = e_2 * e_2
    e_6 = e_4 * e_2

    e_fact = 1.0 - e_2
    e_fact_sqrt = np.sqrt(e_fact)

    M = m1 + m2
    M_fact_2 = M * M
    M_fact_4 = M_fact_2 * M_fact_2

    cos_u = np.cos(u)

    numerator = (
        # (m1 - m2)(m1 + m2)(...)
        (m1 - m2) * M * (
            24 * (m1 * m2 - 3 * M_fact_2)
            - 6 * e_6 * (m1 * m2 - 2 * M_fact_2)
            + 18 * e_4 * (3 * m1 * m2 - 2 * M_fact_2)
            - e_2 * (31 * m1 * m2 + 56 * M_fact_2)
        ) * (S1z - S2z)

        # symmetric spin part
        + (
            e_2 * (
                50 * m1*m1*m2*m2 - 57 * m1*m2*M_fact_2 - 56 * M_fact_4
            )
            + 6 * e_6 * (
                2 * m1*m1*m2*m2 - 3 * m1*m2*M_fact_2 + 2 * M_fact_4
            )
            - 6 * e_4 * (
                8 * m1*m1*m2*m2 - 27 * m1*m2*M_fact_2 + 6 * M_fact_4
            )
            - 12 * (
                m1*m1*m2*m2 - 8 * m1*m2*M_fact_2 + 6 * M_fact_4
            )
        ) * (S1z + S2z)

        # cos(u) term
        - e * (
            -8 * (56 + 49 * e_2 + 9 * e_4) * m1**4 * S1z
            -8 * (56 + 49 * e_2 + 9 * e_4) * m2**4 * S2z
            + 4 * (-217 - 200 * e_2 + 9 * e_4) * m1*m1*m2*m2 * (S1z + S2z)
            - 2 * m1**3 * m2 * (
                (530 + 496 * e_2 + 6 * e_4) * S1z
                + 9 * (14 + 13 * e_2) * S2z
            )
            - 2 * m1 * m2**3 * (
                9 * (14 + 13 * e_2) * S1z
                + 2 * (265 + 248 * e_2 + 3 * e_4) * S2z
            )
        ) * cos_u

        # cos^2(u) term
        + e_2 * (
            -8 * (2 + e_2) * (29 + 9 * e_2) * m1**4 * S1z
            -8 * (2 + e_2) * (29 + 9 * e_2) * m2**4 * S2z
            -8 * (2 + e_2) * (47 + 21 * e_2) * m1*m1*m2*m2 * (S1z + S2z)
            -2 * m1**3 * m2 * (
                (514 + 440 * e_2 + 78 * e_4) * S1z
                + 9 * (2 + e_2) * (5 + 4 * e_2) * S2z
            )
            -2 * m1 * m2**3 * (
                9 * (2 + e_2) * (5 + 4 * e_2) * S1z
                + 2 * (257 + 220 * e_2 + 39 * e_4) * S2z
            )
        ) * (cos_u ** 2)

        # cos^3(u) term
        - e_3 * (
            -8 * (17 + 21 * e_2) * m1**4 * S1z
            -8 * (17 + 21 * e_2) * m2**4 * S2z
            -8 * (11 + 57 * e_2) * m1*m1*m2*m2 * (S1z + S2z)
            -2 * m1**3 * m2 * (
                4 * (29 + 57 * e_2) * S1z - 6 * S2z + 87 * e_2 * S2z
            )
            -2 * m1 * m2**3 * (
                (-6 + 87 * e_2) * S1z + 4 * (29 + 57 * e_2) * S2z
            )
        ) * (cos_u ** 3)

        # sqrt(e_fact) term
        - 12 * e_fact_sqrt * (
            -12 * m1**4 * S1z
            -12 * m2**4 * S2z
            -21 * m1*m1*m2*m2 * (S1z + S2z)
            -2 * m1**3 * m2 * (13 * S1z + 3 * S2z)
            -2 * m1 * m2**3 * (3 * S1z + 13 * S2z)
        ) * ((-1 + e * cos_u) ** 2) * (1 - 2 * e_2 + e * cos_u)
    )

    denominator = (
        12.0
        * ((-1 + e_2) ** 2)
        * M_fact_4
        * ((-1 + e * cos_u) ** 5)
    )

    return numerator / denominator


## --------- 3PN -------------

def phi_dot_3pn(e: float, eta: float, u: float) -> float:
    u_factor = cosu_factor(e, u)
    u_factor_pow_7 = u_factor ** 7
    pi_pow_2 = np.pi ** 2
    eta_pow_2 = eta ** 2
    eta_pow_3 = eta ** 3
    e_pow_2 = e ** 2
    e_fact = 1.0 - e_pow_2
    e_pow_3 = e ** 3
    e_pow_4 = e ** 4
    e_pow_5 = e ** 5
    e_pow_6 = e ** 6
    e_pow_7 = e ** 7
    e_pow_8 = e ** 8
    e_pow_9 = e ** 9
    e_pow_10 = e ** 10
    e_factor = 1.0 - e_pow_2
    cos_u = np.cos(u)
    cos_u_pow_2 = cos_u ** 2
    cos_u_pow_3 = cos_u ** 3
    cos_u_pow_4 = cos_u ** 4
    cos_u_pow_5 = cos_u ** 5

    pre_factor = 1.0 / (13440.0 * np.sqrt(e_factor) * e_factor ** 2 * u_factor_pow_7)

    e_0_term = 67200.0 * eta_pow_2 - 761600.0 * eta + 8610.0 * eta * pi_pow_2 + 201600.0
    e_2_term = (4480.0 * eta_pow_3 - 412160.0 * eta_pow_2
                - 30135.0 * pi_pow_2 * eta + 553008.0 * eta + 342720.0) * e_pow_2
    e_4_term = (-52640.0 * eta_pow_3 + 516880.0 * eta_pow_2
                + 68880.0 * pi_pow_2 * eta - 1916048.0 * eta + 262080.0) * e_pow_4
    e_6_term = (84000.0 * eta_pow_3 - 190400.0 * eta_pow_2
                - 17220.0 * pi_pow_2 * eta - 50048.0 * eta - 241920.0) * e_pow_6
    e_8_term = (-52640.0 * eta_pow_2 - 13440.0 * eta + 483280.0) * eta * e_pow_8
    e_10_term = (10080.0 * eta_pow_2 + 40320.0 * eta - 15120.0) * eta * e_pow_10

    cosu_1 = (
        (-2240.0 * eta_pow_3 - 168000.0 * eta_pow_2 - 424480.0 * eta) * e_pow_9
        + (28560.0 * eta_pow_3 + 242480.0 * eta_pow_2 + 34440.0 * pi_pow_2 * eta
           - 1340224.0 * eta + 725760.0) * e_pow_7
        + (-33040.0 * eta_pow_3 - 754880.0 * eta_pow_2 - 172200.0 * pi_pow_2 * eta
           + 5458480.0 * eta - 221760.0) * e_pow_5
        + (40880.0 * eta_pow_3 + 738640.0 * eta_pow_2 + 30135.0 * pi_pow_2 * eta
           + 1554048.0 * eta - 2936640.0) * e_pow_3
        + (-560.0 * eta_pow_3 - 100240.0 * eta_pow_2 - 43050.0 * pi_pow_2 * eta
           + 3284816.0 * eta - 389760.0) * e
    ) * cos_u

    cosu_2 = (
        (4480.0 * eta_pow_3 - 20160.0 * eta_pow_2 + 16800.0 * eta) * e_pow_10
        + (3920.0 * eta_pow_3 + 475440.0 * eta_pow_2 - 17220.0 * pi_pow_2 * eta
           + 831952.0 * eta - 7257600.0) * e_pow_8
        + (-75600.0 * eta_pow_3 + 96880.0 * eta_pow_2 + 154980.0 * pi_pow_2 * eta
           - 3249488.0 * eta - 685440.0) * e_pow_6
        + (5040.0 * eta_pow_3 - 659120.0 * eta_pow_2 + 25830.0 * pi_pow_2 * eta
           - 7356624.0 * eta + 6948480.0) * e_pow_4
        + (-5040.0 * eta_pow_3 + 190960.0 * eta_pow_2 + 137760.0 * pi_pow_2 * eta
           - 7307920.0 * eta + 107520.0) * e_pow_2
    ) * cos_u_pow_2

    cosu_3 = (
        (560.0 * eta_pow_3 - 137200.0 * eta_pow_2 + 388640.0 * eta + 241920.0) * e_pow_9
        + (30800.0 * eta_pow_3 - 264880.0 * eta_pow_2 - 68880.0 * pi_pow_2 * eta
           + 624128.0 * eta + 766080.0) * e_pow_7
        + (66640.0 * eta_pow_3 + 612080.0 * eta_pow_2 - 8610.0 * pi_pow_2 * eta
           + 6666080.0 * eta - 6652800.0) * e_pow_5
        + (-30800.0 * eta_pow_3 - 294000.0 * eta_pow_2 - 223860.0 * pi_pow_2 * eta
           + 9386432.0 * eta) * e_pow_3
    ) * cos_u_pow_3

    cosu_4 = (
        (-16240.0 * eta_pow_3 + 12880.0 * eta_pow_2 + 18480.0 * eta) * e_pow_10
        + (16240.0 * eta_pow_3 - 91840.0 * eta_pow_2 + 17220.0 * pi_pow_2 * eta
           - 652192.0 * eta + 100800.0) * e_pow_8
        + (-55440.0 * eta_pow_3 + 34160.0 * eta_pow_2 - 30135.0 * pi_pow_2 * eta
           - 2185040.0 * eta + 2493120.0) * e_pow_6
        + (21480.0 * eta_pow_3 + 86800.0 * eta_pow_2 + 163590.0 * pi_pow_2 * eta
           - 5713888.0 * eta + 228480.0) * e_pow_4
    ) * cos_u_pow_4

    cosu_5 = (
        (13440.0 * eta_pow_3 + 94640.0 * eta_pow_2 - 113680.0 * eta - 221760.0) * e_pow_9
        + (-11200.0 * eta_pow_3 - 112000.0 * eta_pow_2 + 12915.0 * pi_pow_2 * eta
           + 692928.0 * eta - 194880.0) * e_pow_7
        + (4480.0 * eta_pow_3 + 8960.0 * eta_pow_2 - 43050.0 * pi_pow_2 * eta
           + 1127280.0 * eta - 147840.0) * e_pow_5
    ) * cos_u_pow_5

    rt_zero = (
        -67200.0 * eta_pow_2 + 761600.0 * eta
        + e_pow_4 * (40320.0 * eta_pow_2 + 309120.0 * eta - 672000.0)
        + e_pow_2 * (208320.0 * eta_pow_2 + 17220.0 * pi_pow_2 * eta
                     - 2289280.0 * eta + 1680000.0)
        - 8610.0 * pi_pow_2 * eta - 201600.0
    )

    rt_cosu_1 = (
        (-282240.0 * eta_pow_2 - 450240.0 * eta + 1478400.0) * e_pow_5
        + (-719040.0 * eta_pow_2 - 68880.0 * pi_pow_2 * eta
           + 8128960.0 * eta - 5040000.0) * e_pow_3
        + (94080.0 * eta_pow_2 + 25830.0 * pi_pow_2 * eta
           - 1585920.0 * eta - 470400.0) * e
    ) * cos_u

    rt_cosu_2 = (
        (604800.0 * eta_pow_2 - 504000.0 * eta - 403200.0) * e_pow_6
        + (1034880.0 * eta_pow_2 + 103320.0 * pi_pow_2 * eta
           - 11195520.0 * eta + 5779200.0) * e_pow_4
        + (174720.0 * eta_pow_2 - 17220.0 * pi_pow_2 * eta
           - 486080.0 * eta + 2688000.0) * e_pow_2
    ) * cos_u_pow_2

    rt_cosu_3 = (
        (-524160.0 * eta_pow_2 + 1122240.0 * eta - 940800.0) * e_pow_7
        + (-873600.0 * eta_pow_2 - 68880.0 * pi_pow_2 * eta
           + 7705600.0 * eta - 3897600.0) * e_pow_5
        + (-416640.0 * eta_pow_2 - 17220.0 * pi_pow_2 * eta
           + 3357760.0 * eta - 3225600.0) * e_pow_3
    ) * cos_u_pow_3

    rt_cosu_4 = (
        (161280.0 * eta_pow_2 - 477120.0 * eta + 537600.0) * e_pow_8
        + (477120.0 * eta_pow_2 + 17220.0 * pi_pow_2 * eta
           - 2894080.0 * eta + 2217600.0) * e_pow_6
        + (268800.0 * eta_pow_2 + 25830.0 * pi_pow_2 * eta
           - 2721600.0 * eta + 1276800.0) * e_pow_4
    ) * cos_u_pow_4

    rt_cosu_5 = (
        (-127680.0 * eta_pow_2 + 544320.0 * eta - 739200.0) * e_pow_7
        + (-53760.0 * eta_pow_2 - 8610.0 * pi_pow_2 * eta
           + 674240.0 * eta - 67200.0) * e_pow_5
    ) * cos_u_pow_5

    return pre_factor * (
        e_0_term + e_2_term + e_4_term + e_6_term + e_8_term + e_10_term
        + cosu_1 + cosu_2 + cosu_3 + cosu_4 + cosu_5
        + np.sqrt(e_fact) * (
            rt_zero + rt_cosu_1 + rt_cosu_2 + rt_cosu_3 + rt_cosu_4 + rt_cosu_5
        )
    )


def phi_dot_3pn_SS(e: float, m1: float, m2: float, S1z: float, S2z: float, u: float) -> float:
    """Quentin Henry et al terms, arXiv:2308.13606v1"""
    kappa1 = 1.0
    kappa2 = 1.0
    e_2 = e * e
    e_3 = e_2 * e
    e_4 = e_2 * e_2
    e_6 = e_4 * e_2
    e_fact = 1.0 - e_2
    e_fact_sqrt = np.sqrt(e_fact)
    M_fact_2 = (m1 + m2) ** 2
    M_fact_4 = M_fact_2 * M_fact_2
    cos_u = np.cos(u)
    cos_u_2 = cos_u * cos_u
    cos_u_3 = cos_u_2 * cos_u

    # --- constant terms (no cos(u) power) ---
    s1z2_m1_4 = 6 * (
        4 * e_6 * e_fact_sqrt * kappa1
        - 2 * (-1 + e_fact_sqrt) * (4 + 7 * kappa1)
        - e_2 * (24 + 20 * e_fact_sqrt + 42 * kappa1 + 19 * e_fact_sqrt * kappa1)
        - 4 * e_4 * (-4 + (-7 + 3 * e_fact_sqrt) * kappa1)
    ) * m1**4 * S1z * S1z

    s2z2_m2_4 = 6 * (
        4 * e_6 * e_fact_sqrt * kappa2
        - 2 * (-1 + e_fact_sqrt) * (4 + 7 * kappa2)
        - e_2 * (24 + 20 * e_fact_sqrt + 42 * kappa2 + 19 * e_fact_sqrt * kappa2)
        - 4 * e_4 * (-4 + (-7 + 3 * e_fact_sqrt) * kappa2)
    ) * m2**4 * S2z * S2z

    cross_m1_3_m2 = m1**3 * m2 * S1z * (
        (
            48 * e_6 * e_fact_sqrt * (-1 + kappa1)
            - 6 * (-1 + e_fact_sqrt) * (3 + 23 * kappa1)
            - 4 * e_4 * (-9 - 60 * e_fact_sqrt - 69 * kappa1 + 32 * e_fact_sqrt * kappa1)
            - e_2 * (54 + 297 * e_fact_sqrt + 414 * kappa1 + 145 * e_fact_sqrt * kappa1)
        ) * S1z
        + 6 * (
            60 * e_4 + 4 * e_6 * e_fact_sqrt - 30 * (-1 + e_fact_sqrt)
            - e_2 * (90 + 67 * e_fact_sqrt)
        ) * S2z
    )

    cross_m1_m2_3 = m1 * m2**3 * S2z * (
        6 * (
            60 * e_4 + 4 * e_6 * e_fact_sqrt - 30 * (-1 + e_fact_sqrt)
            - e_2 * (90 + 67 * e_fact_sqrt)
        ) * S1z
        + (
            48 * e_6 * e_fact_sqrt * (-1 + kappa2)
            - 6 * (-1 + e_fact_sqrt) * (3 + 23 * kappa2)
            - 4 * e_4 * (-9 - 60 * e_fact_sqrt - 69 * kappa2 + 32 * e_fact_sqrt * kappa2)
            - e_2 * (54 + 297 * e_fact_sqrt + 414 * kappa2 + 145 * e_fact_sqrt * kappa2)
        ) * S2z
    )

    cross_m1_2_m2_2 = m1**2 * m2**2 * (
        3 * (
            4 * e_6 * e_fact_sqrt * (-4 + 3 * kappa1)
            - 6 * (-1 + e_fact_sqrt) * (-1 + 4 * kappa1)
            - e_2 * (-18 + 55 * e_fact_sqrt + 4 * (18 + e_fact_sqrt) * kappa1)
            + e_4 * (-12 + 68 * e_fact_sqrt - 8 * (-6 + 5 * e_fact_sqrt) * kappa1)
        ) * S1z * S1z
        + 2 * (
            12 * e_6 * e_fact_sqrt - 174 * (-1 + e_fact_sqrt)
            + 4 * e_4 * (87 + 7 * e_fact_sqrt)
            - e_2 * (522 + 409 * e_fact_sqrt)
        ) * S1z * S2z
        + 3 * (
            4 * e_6 * e_fact_sqrt * (-4 + 3 * kappa2)
            - 6 * (-1 + e_fact_sqrt) * (-1 + 4 * kappa2)
            - e_2 * (-18 + 55 * e_fact_sqrt + 4 * (18 + e_fact_sqrt) * kappa2)
            + e_4 * (-12 + 68 * e_fact_sqrt - 8 * (-6 + 5 * e_fact_sqrt) * kappa2)
        ) * S2z * S2z
    )

    term_const = s1z2_m1_4 + s2z2_m2_4 + cross_m1_3_m2 + cross_m1_m2_3 + cross_m1_2_m2_2

    # --- cos(u) terms ---
    cosu_m1_4_S1z2 = 6 * (
        -8 + 40 * e_fact_sqrt + 2 * (-7 + 25 * e_fact_sqrt) * kappa1
        + 4 * e_4 * (-8 + (-14 + 3 * e_fact_sqrt) * kappa1)
        + e_2 * (40 + 44 * e_fact_sqrt + (70 + 61 * e_fact_sqrt) * kappa1)
    ) * m1**4 * S1z * S1z

    cosu_m2_4_S2z2 = 6 * (
        -8 + 40 * e_fact_sqrt + 2 * (-7 + 25 * e_fact_sqrt) * kappa2
        + 4 * e_4 * (-8 + (-14 + 3 * e_fact_sqrt) * kappa2)
        + e_2 * (40 + 44 * e_fact_sqrt + (70 + 61 * e_fact_sqrt) * kappa2)
    ) * m2**4 * S2z * S2z

    cosu_m1_3_m2 = m1**3 * m2 * S1z * (
        (
            -18 + 246 * e_fact_sqrt + 46 * (-3 + 10 * e_fact_sqrt) * kappa1
            + 2 * e_4 * (-36 - 60 * e_fact_sqrt - 276 * kappa1 + 47 * e_fact_sqrt * kappa1)
            + e_2 * (90 + 243 * e_fact_sqrt + (690 + 535 * e_fact_sqrt) * kappa1)
        ) * S1z
        + 18 * (
            -10 + 42 * e_fact_sqrt + 4 * e_4 * (-10 + e_fact_sqrt)
            + e_2 * (50 + 47 * e_fact_sqrt)
        ) * S2z
    )

    cosu_m1_m2_3 = m1 * m2**3 * S2z * (
        18 * (
            -10 + 42 * e_fact_sqrt + 4 * e_4 * (-10 + e_fact_sqrt)
            + e_2 * (50 + 47 * e_fact_sqrt)
        ) * S1z
        + (
            -18 + 246 * e_fact_sqrt + 46 * (-3 + 10 * e_fact_sqrt) * kappa2
            + 2 * e_4 * (-36 - 60 * e_fact_sqrt - 276 * kappa2 + 47 * e_fact_sqrt * kappa2)
            + e_2 * (90 + 243 * e_fact_sqrt + (690 + 535 * e_fact_sqrt) * kappa2)
        ) * S2z
    )

    cosu_m1_2_m2_2 = m1**2 * m2**2 * (
        3 * (
            6 + 14 * e_fact_sqrt + 8 * (-3 + 8 * e_fact_sqrt) * kappa1
            + 4 * e_4 * (6 - 7 * e_fact_sqrt + 8 * (-3 + e_fact_sqrt) * kappa1)
            + e_2 * (5 * (-6 + e_fact_sqrt) + 24 * (5 + 3 * e_fact_sqrt) * kappa1)
        ) * S1z * S1z
        + 2 * (
            -174 + 760 * e_fact_sqrt
            + e_4 * (-696 + 34 * e_fact_sqrt)
            + e_2 * (870 + 835 * e_fact_sqrt)
        ) * S1z * S2z
        + 3 * (
            6 + 14 * e_fact_sqrt + 8 * (-3 + 8 * e_fact_sqrt) * kappa2
            + 4 * e_4 * (6 - 7 * e_fact_sqrt + 8 * (-3 + e_fact_sqrt) * kappa2)
            + e_2 * (5 * (-6 + e_fact_sqrt) + 24 * (5 + 3 * e_fact_sqrt) * kappa2)
        ) * S2z * S2z
    )

    term_cosu = e * (
        cosu_m1_4_S1z2 + cosu_m2_4_S2z2 + cosu_m1_3_m2 + cosu_m1_m2_3 + cosu_m1_2_m2_2
    ) * cos_u

    # --- cos(u)^2 terms ---
    cosu2_m1_4_S1z2 = 6 * (
        8 + 56 * e_fact_sqrt + 14 * kappa1 + 58 * e_fact_sqrt * kappa1
        + 4 * e_4 * (-4 - 7 * kappa1 + 3 * e_fact_sqrt * kappa1)
        + e_2 * (8 + 28 * e_fact_sqrt + 14 * kappa1 + 53 * e_fact_sqrt * kappa1)
    ) * m1**4 * S1z * S1z

    cosu2_m2_4_S2z2 = 6 * (
        8 + 56 * e_fact_sqrt + 14 * kappa2 + 58 * e_fact_sqrt * kappa2
        + 4 * e_4 * (-4 - 7 * kappa2 + 3 * e_fact_sqrt * kappa2)
        + e_2 * (8 + 28 * e_fact_sqrt + 14 * kappa2 + 53 * e_fact_sqrt * kappa2)
    ) * m2**4 * S2z * S2z

    cosu2_m1_3_m2 = m1**3 * m2 * S1z * (
        (
            2 * e_4 * (-18 - 12 * e_fact_sqrt - 138 * kappa1 + 55 * e_fact_sqrt * kappa1)
            + 2 * (9 + 111 * e_fact_sqrt + 69 * kappa1 + 262 * e_fact_sqrt * kappa1)
            + e_2 * (18 + 171 * e_fact_sqrt + 138 * kappa1 + 455 * e_fact_sqrt * kappa1)
        ) * S1z
        + 18 * (
            10 + 46 * e_fact_sqrt
            + 4 * e_4 * (-5 + 2 * e_fact_sqrt)
            + e_2 * (10 + 39 * e_fact_sqrt)
        ) * S2z
    )

    cosu2_m1_m2_3 = m1 * m2**3 * S2z * (
        18 * (
            10 + 46 * e_fact_sqrt
            + 4 * e_4 * (-5 + 2 * e_fact_sqrt)
            + e_2 * (10 + 39 * e_fact_sqrt)
        ) * S1z
        + (
            2 * e_4 * (-18 - 12 * e_fact_sqrt - 138 * kappa2 + 55 * e_fact_sqrt * kappa2)
            + 2 * (9 + 111 * e_fact_sqrt + 69 * kappa2 + 262 * e_fact_sqrt * kappa2)
            + e_2 * (18 + 171 * e_fact_sqrt + 138 * kappa2 + 455 * e_fact_sqrt * kappa2)
        ) * S2z
    )

    cosu2_m1_2_m2_2 = m1**2 * m2**2 * (
        3 * (
            -6 - 14 * e_fact_sqrt + 24 * kappa1 + 68 * e_fact_sqrt * kappa1
            + 4 * e_4 * (3 - 2 * e_fact_sqrt - 12 * kappa1 + 7 * e_fact_sqrt * kappa1)
            + e_2 * (-6 + 13 * e_fact_sqrt + 24 * kappa1 + 72 * e_fact_sqrt * kappa1)
        ) * S1z * S1z
        + 2 * (
            174 + 872 * e_fact_sqrt
            + e_4 * (-348 + 98 * e_fact_sqrt)
            + e_2 * (174 + 659 * e_fact_sqrt)
        ) * S1z * S2z
        + 3 * (
            -6 - 14 * e_fact_sqrt + 24 * kappa2 + 68 * e_fact_sqrt * kappa2
            + 4 * e_4 * (3 - 2 * e_fact_sqrt - 12 * kappa2 + 7 * e_fact_sqrt * kappa2)
            + e_2 * (-6 + 13 * e_fact_sqrt + 24 * kappa2 + 72 * e_fact_sqrt * kappa2)
        ) * S2z * S2z
    )

    term_cosu2 = -e_2 * (
        cosu2_m1_4_S1z2 + cosu2_m2_4_S2z2 + cosu2_m1_3_m2 + cosu2_m1_m2_3 + cosu2_m1_2_m2_2
    ) * cos_u_2

    # --- cos(u)^3 terms ---
    cosu3_m1_4_S1z2 = 6 * (
        8 + 24 * e_fact_sqrt + 2 * (7 + 9 * e_fact_sqrt) * kappa1
        + e_2 * (-8 + 4 * e_fact_sqrt - 14 * kappa1 + 23 * e_fact_sqrt * kappa1)
    ) * m1**4 * S1z * S1z

    cosu3_m2_4_S2z2 = 6 * (
        8 + 24 * e_fact_sqrt + 2 * (7 + 9 * e_fact_sqrt) * kappa2
        + e_2 * (-8 + 4 * e_fact_sqrt - 14 * kappa2 + 23 * e_fact_sqrt * kappa2)
    ) * m2**4 * S2z * S2z

    cosu3_m1_3_m2 = m1**3 * m2 * S1z * (
        (
            18 + 42 * e_fact_sqrt + 2 * (69 + 77 * e_fact_sqrt) * kappa1
            + e_2 * (-18 + 81 * e_fact_sqrt - 138 * kappa1 + 209 * e_fact_sqrt * kappa1)
        ) * S1z
        + 6 * (
            30 + 38 * e_fact_sqrt + 5 * e_2 * (-6 + 11 * e_fact_sqrt)
        ) * S2z
    )

    cosu3_m1_m2_3 = m1 * m2**3 * S2z * (
        6 * (
            30 + 38 * e_fact_sqrt + 5 * e_2 * (-6 + 11 * e_fact_sqrt)
        ) * S1z
        + (
            18 + 42 * e_fact_sqrt + 2 * (69 + 77 * e_fact_sqrt) * kappa2
            + e_2 * (-18 + 81 * e_fact_sqrt - 138 * kappa2 + 209 * e_fact_sqrt * kappa2)
        ) * S2z
    )

    cosu3_m1_2_m2_2 = m1**2 * m2**2 * (
        3 * (
            -6 - 18 * e_fact_sqrt + 8 * (3 + 2 * e_fact_sqrt) * kappa1
            + e_2 * (6 + 15 * e_fact_sqrt + 8 * (-3 + 5 * e_fact_sqrt) * kappa1)
        ) * S1z * S1z
        + 2 * (
            174 + 274 * e_fact_sqrt + e_2 * (-174 + 269 * e_fact_sqrt)
        ) * S1z * S2z
        + 3 * (
            -6 - 18 * e_fact_sqrt + 8 * (3 + 2 * e_fact_sqrt) * kappa2
            + e_2 * (6 + 15 * e_fact_sqrt + 8 * (-3 + 5 * e_fact_sqrt) * kappa2)
        ) * S2z * S2z
    )

    term_cosu3 = e_3 * (
        cosu3_m1_4_S1z2 + cosu3_m2_4_S2z2 + cosu3_m1_3_m2 + cosu3_m1_m2_3 + cosu3_m1_2_m2_2
    ) * cos_u_3

    denominator = 12.0 * (e_2 - 1.0) ** 3 * M_fact_4 * (1.0 - e * cos_u) ** 5

    phi_3pn_SS = (term_const + term_cosu + term_cosu2 + term_cosu3) / denominator

    return phi_3pn_SS

## ------- 4PN & 4.5PN---------------
def phi_dot_4pn_SS(e: float, m1: float, m2: float,
                    S1z: float, S2z: float) -> float:
    """4PN SS phi_dot — returns 0 (same as active code in original)."""
    return 0.0


def phi_dot_4_5_pn(e: float, eta: float, x: float) -> float:
    """4.5PN phi_dot — returns 0."""
    return 0.0


# ================ Relative separation ======================================

## ---------- 0PN --------------------

def rel_sep_0pn(e: float, u: float) -> float:
    return 1.0 - e * np.cos(u)

## ---------- 1PN --------------------

def rel_sep_1pn(e: float, u: float, eta: float) -> float:
    ef   = 1.0 - e*e
    b1   = 2.0*(1.0 - e*np.cos(u)) / ef
    b2   = (-18.0 + 2.0*eta - e*(6.0 - 7.0*eta)*np.cos(u)) / 6.0
    return b1 + b2

## ---------- 1.5PN --------------------

def rel_sep_1_5pn(e: float, u: float, m1: float, m2: float,
                   S1z: float, S2z: float) -> float:
    """1.5PN SO (Klein et al. arXiv:1801.08542, Eq. B1a/B2a/B2b)"""
    e2  = e*e
    ef  = 1.0 - e2
    M   = m1 + m2
    return (
        -1.0/3.0 *
        (2*m1*m1*(S1z + 3*e2*S1z) +
         2*(1+3*e2)*m2*m2*S2z +
         3*(1+e2)*m1*m2*(S1z+S2z) -
         2*e*(4*m1*m1*S1z + 4*m2*m2*S2z + 3*m1*m2*(S1z+S2z))*np.cos(u))
        / (ef**1.5 * M**2)
    )

## ---------- 2PN --------------------

def rel_sep_2pn(e: float, u: float, eta: float) -> float:
    eta2  = eta*eta
    e2    = e*e
    ef    = 1.0 - e2
    n1    = (-48.0 + 28.0*eta + e2*(-51.0+26.0*eta))*(-1.0+e*np.cos(u))
    d1    = 6.0*ef**2
    n2    = (72.0*(-4.0+7.0*eta) +
             36.0*np.sqrt(ef)*(-5.0+2.0*eta)*(2.0+e*np.cos(u)) +
             ef*(72.0+30.0*eta+8.0*eta2 + e*(-72.0+7*(33.0-5.0*eta)*eta)*np.cos(u)))
    d2    = 72.0*ef
    return n1/d1 + n2/d2


def rel_sep_2pnSS(e: float, u: float, m1: float, m2: float,
                   S1z: float, S2z: float) -> float:
    """2PN SS relative separation"""
    kappa1 = 1.0; kappa2 = 1.0
    return (
        (m1*S1z*(2*m2*S2z + m1*S1z*kappa1) + m2*m2*S2z*S2z*kappa2) *
        (1 + e*e - 2*e*np.cos(u))
    ) / (2.0 * (-1+e**2)**2 * (m1+m2)**2)

## ---------- 2.5PN --------------------

def rel_sep_2_5pn_SO(e: float, u: float, m1: float, m2: float,
                      S1z: float, S2z: float) -> float:
    """2.5PN SO relative separation."""
    e2   = e*e; e4=e2*e2
    ef   = 1.0-e2; ef_sqrt=np.sqrt(ef)
    M2   = (m1+m2)**2; M4=M2*M2
    return (
        (2*(-1+e2)**2 *
         (-12*m1**4*S1z - 12*m2**4*S2z -
          21*m1*m1*m2*m2*(S1z+S2z) -
          2*m1**3*m2*(13*S1z+3*S2z) -
          2*m1*m2**3*(3*S1z+13*S2z)) * (2+e*np.cos(u)) +
         2*ef_sqrt *
         (12*(2+7*e2+e4)*m1**4*S1z +
          12*(2+7*e2+e4)*m2**4*S2z +
          2*(22+85*e2+10*e4)*m1*m1*m2*m2*(S1z+S2z) +
          m1**3*m2*(2*(28+96*e2+11*e4)*S1z + 3*(4+19*e2+2*e4)*S2z) +
          m1*m2**3*(3*(4+19*e2+2*e4)*S1z + 2*(28+96*e2+11*e4)*S2z) -
          e*(60*(1+e2)*m1**4*S1z + 60*(1+e2)*m2**4*S2z +
             3*(35+43*e2)*m1*m1*m2*m2*(S1z+S2z) +
             m1**3*m2*(133*S1z+137*e2*S1z+30*S2z+45*e2*S2z) +
             m1*m2**3*(30*S1z+45*e2*S1z+133*S2z+137*e2*S2z))*np.cos(u)))
    ) / (6.0 * (-1+e2)**3 * M4)

## ---------- 3PN --------------------

def rel_sep_3pn(e: float, u: float, eta: float) -> float:
    pi_pow_2 = np.pi * np.pi

    e_factor = 1.0 - e * e
    eta_pow_2 = eta * eta
    e_pow_2 = e * e
    e_pow_4 = e_pow_2 * e_pow_2

    def pow7_2(x):
        return x ** 3.5

    term1 = (
        (-665280.0 * eta_pow_2 + 1753920.0 * eta - 1814400.0) * e_pow_4
        + (
            (725760.0 * eta_pow_2 - 77490.0 * pi_pow_2 + 5523840.0) * eta
            - 3628800.0
        ) * e_pow_2
        + (
            (544320.0 * eta_pow_2 + 154980.0 * pi_pow_2 - 14132160.0) * eta
            + 7257600.0
        )
    ) * e_pow_2

    term2 = -604800.0 * eta_pow_2 + 6854400.0 * eta

    term3_cos = (
        (
            (302400.0 * eta_pow_2 - 1254960.0 * eta + 453600.0) * e_pow_4
            + (
                (-1542240.0 * eta_pow_2 - 38745.0 * pi_pow_2 + 6980400.0) * eta
                - 453600.0
            ) * e_pow_2
            + (
                (2177280.0 * eta_pow_2 + 77490.0 * pi_pow_2 - 12373200.0) * eta
                + 4989600.0
            )
        ) * e * e_pow_2
        + (
            (-937440.0 * eta_pow_2 - 37845.0 * pi_pow_2 + 6647760.0) * eta
            - 4989600.0
        ) * e
    ) * np.cos(u)

    term4_sqrt = np.sqrt(e_factor) * (
        (
            (
                (-4480.0 * eta_pow_2 - 25200.0 * eta + 22680.0) * eta
                - 120960.0
            ) * e_pow_4
            + (
                13440.0 * eta_pow_2 * eta
                + 4404960.0 * eta * eta
                + 116235.0 * pi_pow_2
                - 12718296.0 * eta
                + 5261760.0
            ) * e_pow_2
            + (
                (-13440.0 * eta_pow_2 + 2242800.0 * eta + 348705.0 * pi_pow_2
                 - 19225080.0) * eta
                + 1614160.0
            )
        ) * e_pow_2
        + (
            4480.0 * eta_pow_2 + 45360.0 * eta - 8600904.0
        ) * eta
        + (
            (
                (
                    (-6860.0 * eta_pow_2 - 550620.0 * eta - 986580.0) * eta
                    + 120960.0
                ) * e_pow_4
                + (
                    (20580.0 * eta_pow_2 - 2458260.0 * eta + 3458700.0) * eta
                    - 2358720.0
                ) * e_pow_2
                + (
                    (-20580.0 * eta_pow_2 - 3539340.0 * eta
                     - 116235.0 * pi_pow_2 + 20173860.0) * eta
                    - 16148160.0
                )
            ) * e * e_pow_2
            + (
                (6860.0 * eta_pow_2 - 1220940.0 * eta
                 + 464940.0 * pi_pow_2 + 17875620.0) * eta
                - 417440.0
            ) * e
        ) * np.cos(u)
        + 116235.0 * pi_pow_2 * eta
        + 1814400.0
    )

    term5 = -77490.0 * pi_pow_2 * eta - 1814400.0

    numerator = term1 + term2 + term3_cos + term4_sqrt + term5

    denominator = 181440.0 * pow7_2(e_factor)

    return numerator / denominator

def rel_sep_3pn_SS(e: float, u: float, m1: float, m2: float, S1z: float, S2z: float) -> float:
    kappa1 = 1.0
    kappa2 = 1.0
    e_2 = e * e
    e_4 = e_2 * e_2
    e_fact = 1.0 - e_2
    e_fact_sqrt = np.sqrt(e_fact)
    e_fact_35 = e_fact_sqrt * e_fact * e_fact * e_fact
    M_fact_2 = (m1 + m2) * (m1 + m2)
    M_fact_4 = M_fact_2 * M_fact_2
    cos_u = np.cos(u)
    e_fact_m1 = (-1 + e_2) ** 2  # i.e. (e^2 - 1)^2, used for pow(-1+e_2, 2)

    # --- constant terms ---
    m1_4_S1z2 = (
        -48 + 40 * e_fact_sqrt - 84 * kappa1 + 99 * e_fact_sqrt * kappa1
        + 6 * e_2 * (16 + 28 * e_fact_sqrt + 28 * kappa1 + 51 * e_fact_sqrt * kappa1)
        + e_4 * (-48 + (-84 + 45 * e_fact_sqrt) * kappa1)
    ) * m1**4 * S1z * S1z

    m2_4_S2z2 = (
        -48 + 40 * e_fact_sqrt - 84 * kappa2 + 99 * e_fact_sqrt * kappa2
        + 6 * e_2 * (16 + 28 * e_fact_sqrt + 28 * kappa2 + 51 * e_fact_sqrt * kappa2)
        + e_4 * (-48 + (-84 + 45 * e_fact_sqrt) * kappa2)
    ) * m2**4 * S2z * S2z

    m1_3_m2 = 3 * m1**3 * m2 * S1z * (
        (
            -6 + 4 * e_fact_sqrt - 46 * kappa1 + 48 * e_fact_sqrt * kappa1
            + e_4 * (-6 - 6 * e_fact_sqrt - 46 * kappa1 + 25 * e_fact_sqrt * kappa1)
            + e_2 * (12 + 55 * e_fact_sqrt + 92 * kappa1 + 158 * e_fact_sqrt * kappa1)
        ) * S1z
        + 2 * (
            -30 + 29 * e_fact_sqrt
            + 5 * e_2 * (12 + 25 * e_fact_sqrt + 3 * e_2 * (-2 + e_fact_sqrt))
        ) * S2z
    )

    m1_m2_3 = 3 * m1 * m2**3 * S2z * (
        2 * (
            -30 + 29 * e_fact_sqrt
            + 5 * e_2 * (12 + 25 * e_fact_sqrt + 3 * e_2 * (-2 + e_fact_sqrt))
        ) * S1z
        + (
            -6 + 4 * e_fact_sqrt - 46 * kappa2 + 48 * e_fact_sqrt * kappa2
            + e_4 * (-6 - 6 * e_fact_sqrt - 46 * kappa2 + 25 * e_fact_sqrt * kappa2)
            + e_2 * (12 + 55 * e_fact_sqrt + 92 * kappa2 + 158 * e_fact_sqrt * kappa2)
        ) * S2z
    )

    m1_2_m2_2 = m1**2 * m2**2 * (
        9 * (
            2 - 2 * e_fact_sqrt - 8 * kappa1 + 6 * e_fact_sqrt * kappa1
            + e_4 * (2 - 2 * e_fact_sqrt - 8 * kappa1 + 6 * e_fact_sqrt * kappa1)
            + e_2 * (-4 + 3 * e_fact_sqrt + 4 * (4 + 7 * e_fact_sqrt) * kappa1)
        ) * S1z * S1z
        + 2 * (
            -174 + 175 * e_fact_sqrt
            + 348 * e_2 * (1 + 2 * e_fact_sqrt)
            + 6 * e_4 * (-29 + 11 * e_fact_sqrt)
        ) * S1z * S2z
        + 9 * (
            2 - 2 * e_fact_sqrt - 8 * kappa2 + 6 * e_fact_sqrt * kappa2
            + e_4 * (2 - 2 * e_fact_sqrt - 8 * kappa2 + 6 * e_fact_sqrt * kappa2)
            + e_2 * (-4 + 3 * e_fact_sqrt + 4 * (4 + 7 * e_fact_sqrt) * kappa2)
        ) * S2z * S2z
    )

    term_const = -(m1_4_S1z2 + m2_4_S2z2 + m1_3_m2 + m1_m2_3 + m1_2_m2_2)

    # --- cos(u) terms ---
    cosu_m1_4_S1z2 = 2 * (
        12 + 68 * e_fact_sqrt + 3 * (7 + 36 * e_fact_sqrt) * kappa1
        + 3 * e_4 * (4 + 7 * kappa1)
        + 3 * e_2 * (-8 + 12 * e_fact_sqrt - 14 * kappa1 + 39 * e_fact_sqrt * kappa1)
    ) * m1**4 * S1z * S1z

    cosu_m2_4_S2z2 = 2 * (
        12 + 68 * e_fact_sqrt + 3 * (7 + 36 * e_fact_sqrt) * kappa2
        + 3 * e_4 * (4 + 7 * kappa2)
        + 3 * e_2 * (-8 + 12 * e_fact_sqrt - 14 * kappa2 + 39 * e_fact_sqrt * kappa2)
    ) * m2**4 * S2z * S2z

    cosu_m1_3_m2 = 3 * m1**3 * m2 * S1z * (
        (
            20 * e_fact_sqrt + 33 * e_2 * e_fact_sqrt
            + 110 * e_fact_sqrt * kappa1 + 121 * e_2 * e_fact_sqrt * kappa1
            + e_fact_m1 * (3 + 23 * kappa1)
        ) * S1z
        + (152 * e_fact_sqrt + 186 * e_2 * e_fact_sqrt + 30 * e_fact_m1) * S2z
    )

    cosu_m1_m2_3 = 3 * m1 * m2**3 * S2z * (
        (152 * e_fact_sqrt + 186 * e_2 * e_fact_sqrt + 30 * e_fact_m1) * S1z
        + (
            20 * e_fact_sqrt + 33 * e_2 * e_fact_sqrt
            + 110 * e_fact_sqrt * kappa2 + 121 * e_2 * e_fact_sqrt * kappa2
            + e_fact_m1 * (3 + 23 * kappa2)
        ) * S2z
    )

    cosu_m1_2_m2_2 = m1**2 * m2**2 * (
        9 * (
            e_4 * (-1 + 4 * kappa1)
            + (1 + 4 * e_fact_sqrt) * (-1 + 4 * kappa1)
            + e_2 * (2 + 3 * e_fact_sqrt - 8 * kappa1 + 24 * e_fact_sqrt * kappa1)
        ) * S1z * S1z
        + 2 * (
            87 + 87 * e_4 + 466 * e_fact_sqrt
            + 3 * e_2 * (-58 + 157 * e_fact_sqrt)
        ) * S1z * S2z
        + 9 * (
            e_4 * (-1 + 4 * kappa2)
            + (1 + 4 * e_fact_sqrt) * (-1 + 4 * kappa2)
            + e_2 * (2 + 3 * e_fact_sqrt - 8 * kappa2 + 24 * e_fact_sqrt * kappa2)
        ) * S2z * S2z
    )

    term_cosu = e * (
        cosu_m1_4_S1z2 + cosu_m2_4_S2z2 + cosu_m1_3_m2 + cosu_m1_m2_3 + cosu_m1_2_m2_2
    ) * cos_u

    r_3pn_SS = -0.05555555555555555 * (term_const + term_cosu) / (e_fact_35 * M_fact_4)

    return r_3pn_SS


# ============ Utility functions =========================================

def SymMassRatio(q: float) -> float:
    """Symmetric mass ratio from mass ratio q (can be >1 or <1)."""
    return q / (1.0 + q)**2


def SmallMassRatio(eta: float) -> float:
    """q<1 mass ratio from symmetric mass ratio eta."""
    return (1.0 - 2.0*eta - np.sqrt(1.0 - 4.0*eta)) / (2.0*eta)


# ================== ODE right-hand-side dispatchers ==========================

def dx_dt(radiation_pn_order: int,
          eta: float, m1: float, m2: float, S1z: float, S2z: float,
          x: float, e: float) -> float:

    x2 = x * x
    x3 = x2 * x
    xsqx = x * np.sqrt(x)
    x5 = x**5

    # -------------------------
    # Instantaneous PN part
    # -------------------------
    inst = x_dot_0pn(e, eta)

    if radiation_pn_order >= 2:
        inst += x_dot_1pn(e, eta) * x

    if radiation_pn_order >= 3:
        inst += x_dot_1_5_pn(e, eta, m1, m2, S1z, S2z) * xsqx

    if radiation_pn_order >= 4:
        inst += x_dot_2pn(e, eta, x) * x2
        inst += x_dot_2pn_SS(e, eta, m1, m2, S1z, S2z) * x2

    if radiation_pn_order >= 5:
        inst += x_dot_2_5pn_SO(e, eta, m1, m2, S1z, S2z) * x2 * np.sqrt(x)

    if radiation_pn_order >= 6:
        inst += x_dot_3pn(e, eta, x) * x3
        inst += x_dot_3pn_SO(e, eta, m1, m2, S1z, S2z) * x3
        inst += x_dot_3pn_SS(e, eta, m1, m2, S1z, S2z) * x3

    if radiation_pn_order >= 7:
        inst += x_dot_3_5pnSO(e, eta, m1, m2, S1z, S2z) * x3 * np.sqrt(x)
        inst += x_dot_3_5_pn(e, eta) * x3 * np.sqrt(x)
        inst += x_dot_3_5pn_SS(e, eta, m1, m2, S1z, S2z) * x3 * np.sqrt(x)
        inst += x_dot_3_5pn_cubicSpin(e, eta, m1, m2, S1z, S2z) * x3 * np.sqrt(x)
    
    if radiation_pn_order >= 8:
        inst += (
            x_dot_4pn(e, eta, x)
            + x_dot_4pnSO(e, eta, m1, m2, S1z, S2z)
            + x_dot_4pnSS(e, eta, m1, m2, S1z, S2z)
        ) * (x2 * x2)

    if radiation_pn_order >= 9:
        inst += x_dot_4_5_pn(e, eta, x) * (x2 * x2) * np.sqrt(x)

    if radiation_pn_order >= 10:
        inst += 0.0 #dxdt_5pn(x, eta), not implemented

    if radiation_pn_order >= 11:
        inst += 0.0 #dxdt_5_5pn(x, eta), not implemented

    if radiation_pn_order >= 12:
        inst += 0.0 #dxdt_6pn(x, eta), not implemented 

    # multiply ONLY instantaneous part
    result = inst * x5

    # -------------------------
    # Hereditary part (NO x^5)
    # -------------------------
    if radiation_pn_order >= 3:
        result += x_dot_hereditary_1_5(e, eta, x)

    if radiation_pn_order >= 5:
        result += x_dot_hereditary_2_5(e, eta, x)

    if radiation_pn_order >= 6:
        result += x_dot_hereditary_3(e, eta, x)

    return result


def de_dt(radiation_pn_order: int,
          eta: float, m1: float, m2: float, S1z: float, S2z: float,
          x: float, e: float) -> float:

    x2 = x * x
    x3 = x2 * x
    x4 = x2 * x2
    xsqx = x * np.sqrt(x)  

    # -------------------------
    # Instantaneous PN part
    # -------------------------
    inst = e_dot_0pn(e, eta)

    if radiation_pn_order >= 1:
        inst += e_dot_1pn(e, eta) * x

    if radiation_pn_order >= 3:
        inst += e_dot_1_5pn_SO(e, m1, m2, S1z, S2z) * xsqx

    if radiation_pn_order >= 4:
        inst += e_dot_2pn(e, eta) * x2
        inst += e_dot_2pn_SS(e, m1, m2, S1z, S2z) * x2

    if radiation_pn_order >= 5:
        inst += e_dot_2_5pn_SO(e, m1, m2, S1z, S2z) * x2 * np.sqrt(x)

    if radiation_pn_order >= 6:
        inst += e_dot_3pn(e, eta, x) * x3
        inst += e_dot_3_5pn(e, eta) * x3 * np.sqrt(x)
        inst += e_dot_3pn_SO(e, m1, m2, S1z, S2z) * x3
        inst += e_dot_3pn_SS(e, m1, m2, S1z, S2z) * x3

    if radiation_pn_order >= 7:
        inst += 0.0 # placeholder for higher order terms, not implemented

    if radiation_pn_order >= 8:
        inst += 0.0 # placeholder for higher order terms, not implemented

    if radiation_pn_order >= 9:
        inst += 0.0 # placeholder for higher order terms, not implemented

    if radiation_pn_order >= 10:
        inst += 0.0 # placeholder for higher order terms, not implemented

    if radiation_pn_order >= 11:
        inst += 0.0 # placeholder for higher order terms, not implemented

    if radiation_pn_order >= 12:
        inst += 0.0 # placeholder for higher order terms, not implemented

    # multiply only instantaneous part
    result = inst * x4

    # -------------------------
    # Hereditary part (NOT multiplied by x^4)
    # -------------------------
    if radiation_pn_order >= 3:
        result += e_rad_hereditary_1_5(e, eta, x)

    if radiation_pn_order >= 5:
        result += e_rad_hereditary_2_5(e, eta, x)

    if radiation_pn_order >= 6:
        result += e_rad_hereditary_3(e, eta, x)

    return result


def dl_dt(eta: float, m1: float, m2: float, S1z: float, S2z: float,
          x: float, e: float) -> float:
    """dl/dt — 3PN accurate with spin corrections."""
    x32 = np.sqrt(x) * x
    return (
        (1.0 +
         x * l_dot_1pn(e, eta) +
         x32 * l_dot_1_5pn_SO(e, m1, m2, S1z, S2z) +
         x*x * l_dot_2pn(e, eta) +
         x*x * l_dot_2pn_SS(e, m1, m2, S1z, S2z) +
         l_dot_2_5pn_SO(e, m1, m2, S1z, S2z)*x32*x +
         x**3 * l_dot_3pn(e, eta) +
         x**3 * l_dot_3pn_SS(e, m1, m2, S1z, S2z)) * x32
    )


def dphi_dt(u: float, eta: float, m1: float, m2: float, S1z: float, S2z: float,
            x: float, e: float) -> float:
    """dphi/dt — 4PN accurate."""
    x32 = np.sqrt(x) * x
    return (
        (phi_dot_0pn(e, eta, u) +
         x * phi_dot_1pn(e, eta, u) +
         x32 * phi_dot_1_5_pnSO_ecc(e, m1, m2, S1z, S2z, u) +
         x*x * phi_dot_2_pnSS_ecc(e, m1, m2, S1z, S2z, u) +
         x*x * phi_dot_2pn(e, eta, u) +
         x32*x * phi_dot_2_5pn_SO(e, m1, m2, S1z, S2z, u) +
         x**3 * phi_dot_3pn(e, eta, u) +
         x32*x32 * phi_dot_3pn_SS(e, m1, m2, S1z, S2z, u) +
         phi_dot_4pn_SS(e, m1, m2, S1z, S2z)*x32*x*x32 +
         phi_dot_4_5_pn(e, eta, x)*x32**3) * x32
    )


# ================== Kepler equation solvers ==========================

def pow1_3(x): return x ** (1.0/3.0)
def pow3(x): return x * x * x
def pow5(x): return x * x * x * x * x

def mikkola_finder(eccentricity, mean_anomaly):
    """
    Solves Kepler's equation using Mikkola's method for an initial guess.
    """
    # Range reduction of mean_anomaly to [-pi, pi]
    while mean_anomaly >  np.np.pi:
        mean_anomaly -= 2 *  np.np.pi
    while mean_anomaly < - np.np.pi:
        mean_anomaly += 2 *  np.np.pi

    # Compute the sign of l
    sgn_mean_anomaly = 1.0 if mean_anomaly >= 0.0 else -1.0
    mean_anomaly *= sgn_mean_anomaly

    # compute alpha and beta of Mikkola Eq. (9a)
    a = (1.0 - eccentricity) / (4.0 * eccentricity + 0.5)
    b = (0.5 * mean_anomaly) / (4.0 * eccentricity + 0.5)

    # compute the sign of beta needed in Eq. (9b)
    sgn_b = 1.0 if b >= 0.0 else -1.0

    # Mikkola Eq. (9b)
    z = pow1_3(b + sgn_b * np.sqrt(b * b + a * a * a))
    
    # Mikkola Eq. (9c)
    s = z - a / z
    
    # add the correction given in Mikkola Eq. (7)
    s = s - 0.078 * pow5(s) / (1.0 + eccentricity)
    
    # finally Mikkola Eq. (8) gives u
    ecc_anomaly = mean_anomaly + eccentricity * (3.0 * s - 4.0 * pow3(s))
    
    # correct the sign of u
    return ecc_anomaly * sgn_mean_anomaly


def pn_kepler_equation(eta, x, e, l):
    """
    3PN accurate Kepler equation solver using Newton's method.
    """
    mean_anom_negative = False
    u = 0.0
    tol = 1.0e-12

    if l == 0:
        return u

    # range reduction of the l
    while l >  np.np.pi:
        l -= 2 *  np.np.pi
    while l < - np.np.pi:
        l += 2 *  np.np.pi

    # solve in the positive part of the orbit
    if l < 0.0:
        l = -l
        mean_anom_negative = True

    newt_thresh = tol * abs(1.0 - e)
    
    # use mikkola for a guess
    u = mikkola_finder(e, l)

    # high eccentricity case
    if (e > 0.8) and (l <  np.np.pi / 3.0):
        trial = l / abs(1.0 - e)
        if trial * trial > 6.0 * abs(1.0 - e):
            # cubic term is dominant
            if l <  np.np.pi:
                trial = pow1_3(6.0 * l)
        u = trial

    # iterate using Newton's method to get solution
    if e < 1.0:
        newt_err = u - e * np.sin(u) - l
        while abs(newt_err) > newt_thresh:
            u -= newt_err / (1.0 - e * np.cos(u))
            newt_err = u - e * np.sin(u) - l
            
    return -u if mean_anom_negative else u

# ================== ODE system for eccentric models ==========================

def eccentric_x_model_odes(t, y, params):
    """
    ODE system for eccentric gravitational wave models.
    y = [x, e, l, phi]
    """
    # Assuming params is an object or dictionary
    eta = params.eta
    m1 = params.m1
    m2 = params.m2
    S1z = params.S1z
    S2z = params.S2z
    radiation_pn_order = params.radiation_pn_order

    # Input variables
    x, e, l = y[0], y[1], y[2]

    # Calculate eccentric anomaly
    u = pn_kepler_equation(eta, x, e, l)

    dydt = [0.0, 0.0, 0.0, 0.0]

    if e != 0:
        dydt[0] = dx_dt(radiation_pn_order, eta, m1, m2, S1z, S2z, x, e)
        dydt[1] = de_dt(radiation_pn_order, eta, m1, m2, S1z, S2z, x, e)
        dydt[2] = dl_dt(eta, m1, m2, S1z, S2z, x, e)
        dydt[3] = dphi_dt(u, eta, m1, m2, S1z, S2z, x, e)
    else:
        # zero eccentricity limit (arXiv:0909.0066)
        dydt[0] = dx_dt(radiation_pn_order, eta, m1, m2, S1z, S2z, x, e)
        dydt[1] = de_dt(radiation_pn_order, eta, m1, m2, S1z, S2z, x, e)
        dydt[2] = dl_dt(eta, m1, m2, S1z, S2z, x, e)
        dydt[3] = x * np.sqrt(x)

    # Check for NaN (equivalent to XLAL_REAL8_FAIL_NAN)
    if  np.isnan(dydt[0]) or  np.isnan(dydt[1]):
        # Raise an exception or return a specific error code
        raise ValueError("ODE derivative calculated as NaN")

    return dydt

