import numpy as np
import lal


def freq(v, MTOT):
    """
    convert velocity to frequency
    """
    return v**3 / (np.pi * MTOT * lal.MTSUN_SI)


def f0_from_a0(a0, m1, m2):
    """
    Get 22-mode frequency from initial separation
    """
    return (
        2
        / 2
        / np.pi
        * np.sqrt(lal.G_SI * (m1 + m2) * lal.MSUN_SI / a0**3 / lal.AU_SI**3)
    )


def time_to_coalescence_circular(f0, m1, m2):
    """
    Get 22-mode frequency from initial separation
    """
    M = (m1 + m2) * lal.MTSUN_SI
    eta = m1 * m2 / (m1 + m2) ** 2
    v = (np.pi * M * f0) ** (1 / 3)

    return 5 / 256 * M / eta * v ** (-8)


def nmax_fit_Wen(e):
    """
    Get peak harmonic from eccentricity using Wen 2003
    """
    return 2 * (1 + e) ** 1.1954 / (1 - e**2) ** 1.5


def nmax_fit_Hamers(e):
    """
    Get peak harmonic from eccentricity using Hamers 2021
    """
    e = np.atleast_1d(e)
    dict_of_ck = {1: -1.01678, 2: 5.57372, 3: -4.9271, 4: 1.68506}
    sum_arr = np.asarray(e)
    for k in dict_of_ck.keys():
        sum_arr = dict_of_ck[k] * e**k

    return 2 * (1 + sum_arr) / (1 - e**2) ** 1.5


def evolution_equations_0PN(t, y, MTOT, eta):
    """
    Evolution equations without spin at leading (0 PN) order
    """
    MT = MTOT * lal.MTSUN_SI
    et, l, v = y
    et2 = et * et
    et4 = et2 * et2

    OTS = np.sqrt(1.0 - et2)
    OTS5 = OTS * OTS * OTS * OTS * OTS
    OTS7 = OTS5 * OTS * OTS

    v2 = v * v
    v3 = v2 * v
    v8 = v3 * v3 * v2
    v9 = v8 * v

    # Eccentricity Evolution

    CFdedt = -et * eta * v8 / (MT)
    POLdedt0 = (304.0 + 121.0 * et2) / (15.0 * OTS5)

    ans_e = CFdedt * POLdedt0

    # Mean Anomaly Evolution

    CFdldt = v3 / (MT)
    POLdldt0 = 1.0

    ans_l = CFdldt * POLdldt0

    # Velocity Evolution

    CFdvdt = eta * v9 / (2.0 * MT)
    POLdvdt0 = (192.0 + 584.0 * et2 + 74 * et4) / (15.0 * OTS7)

    ans_v = CFdvdt * POLdvdt0

    return ans_e, ans_l, ans_v
