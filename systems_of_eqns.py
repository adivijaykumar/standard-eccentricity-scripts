import numpy as np

MTSUN_SI = 4.925491025543576e-06
LAL_PI = 3.141592653589793
LAL_G_SI = 6.67430e-11  # Gravitational constant, N m^2 kg^-2
LAL_GMSUN_SI = 1.3271244 * 10 ** (20)  # Nominal solar mass parameter, m^3 s^-2
LAL_MSUN_SI = 1.9884099021470415 * 10 ** (30)
LAL_AU_SI = 149597870700e0  # Astronomical unit, m


def evolution_equations_0PN(t, y, MTOT, eta):
    """
    Evolution equations without spin at leading (0 PN) order
    """
    MT = MTOT * MTSUN_SI
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


def evolution_equations_2PN(t, y, MTOT, eta):
    """
    Evolution equations without spin at 2 PN order
    """

    MT = MTOT * MTSUN_SI
    et, l, v = y
    et2 = et * et
    et4 = et2 * et2
    et6 = et4 * et2
    et8 = et4 * et4
    et10 = et8 * et2
    et12 = et10 * et2
    et14 = et12 * et2

    eta2 = eta * eta

    OTS = np.sqrt(1.0 - et2)
    OTS2 = OTS * OTS
    OTS4 = OTS2 * OTS2
    OTS5 = OTS * OTS4
    OTS7 = OTS5 * OTS2
    OTS9 = OTS7 * OTS2
    OTS11 = OTS9 * OTS2

    v2 = v * v
    v3 = v2 * v
    v4 = v3 * v
    v8 = v4 * v4
    v9 = v8 * v

    # Eccentricity Evolution

    CFdedt = -et * eta * v8 / (MT)
    POLdedt0 = (304.0 + 121.0 * et2) / (15.0 * OTS5)
    POLdedt2 = -(
        (
            67608.0
            + 228704.0 * eta
            + et4 * (-125361.0 + 93184.0 * eta)
            + et2 * (-718008.0 + 651252.0 * eta)
        )
        * v2
    ) / (2520.0 * OTS7)
    FACdedt3 = (394.0 / 3.0) * LAL_PI * v3
    RFTdedt3 = (
        192.0
        / 985.0
        * OTS
        * (
            OTS
            * (
                1.0
                + 7.260831042 * et2
                + 5.844370473 * et4
                + 0.8452020270 * et6
                + 0.07580633432 * et8
                + 0.002034045037 * et10
            )
            / (
                1.0
                - 4.900627291 * et2
                + 9.512155497 * et4
                - 9.051368575 * et6
                + 4.096465525 * et8
                - 0.5933309609 * et10
                - 0.05427399445 * et12
                - 0.009020225634 * et14
            )
            - (
                1.0
                + 1.893242666 * et2
                - 2.708117333 * et4
                + 0.6192474531 * et6
                + 0.0500847462 * et8
                - 0.01059040781 * et10
            )
            / (
                1.0
                - 4.638007334 * et2
                + 8.716680569 * et4
                - 8.451197591 * et6
                + 4.435922348 * et8
                - 1.199023304 * et10
                + 0.1398678608 * et12
                - 0.004254544193 * et14
            )
        )
        / et2
    )
    POLdedt4 = (
        (
            -15238352.0
            + 12823920.0 * eta
            + 4548096.0 * eta2
            + et6 * (3786543.0 - 4344852.0 * eta + 2758560.0 * eta2)
            + et4 * (46566110.0 - 78343602.0 * eta + 42810096 * eta2)
            + et2 * (-37367868 - 41949252 * eta + 48711348 * eta2)
            - 1008.0 * (2672.0 + 6963.0 * et2 + 565.0 * et4) * (-5.0 + 2.0 * eta) * OTS
        )
        * v4
    ) / (30240.0 * OTS9)

    ans_e = CFdedt * (POLdedt0 + POLdedt2 + (FACdedt3 * RFTdedt3) + POLdedt4)

    # Mean Anomaly Evolution

    CFdldt = v3 / MT
    POLdldt0 = 1.0
    POLdldt2 = -(3.0 * v2) / OTS2
    POLdldt4 = ((-18.0 + 28.0 * eta + et2 * (-51.0 + 26.0 * eta)) * v4) / (4.0 * OTS4)

    ans_l = CFdldt * (POLdldt0 + POLdldt2 + POLdldt4)

    # Velocity Evolution

    CFdvdt = eta * v9 / (2.0 * MT)
    POLdvdt0 = (192.0 + 584.0 * et2 + 74 * et4) / (15.0 * OTS7)
    POLdvdt2 = (
        (
            -11888.0
            + et2 * (87720.0 - 159600.0 * eta)
            + et4 * (171038.0 - 141708.0 * eta)
            + et6 * (11717.0 - 8288.0 * eta)
            - 14784.0 * eta
        )
        * v2
    ) / (420.0 * OTS9)
    FACdvdt3 = (256.0 / 5.0) * LAL_PI * v3
    RFTdvdt3 = (
        1.0
        + 7.260831042 * et2
        + 5.844370473 * et4
        + 0.8452020270 * et6
        + 0.07580633432 * et8
        + 0.002034045037 * et10
    ) / (
        1.0
        - 4.900627291 * et2
        + 9.512155497 * et4
        - 9.051368575 * et6
        + 4.096465525 * et8
        - 0.5933309609 * et10
        - 0.05427399445 * et12
        - 0.009020225634 * et14
    )
    POLdvdt4 = (
        (
            -360224.0
            + 4514976.0 * eta
            + 1903104.0 * eta2
            + et8 * (3523113.0 - 3259980.0 * eta + 1964256.0 * eta2)
            + et2 * (-92846560.0 + 15464736.0 * eta + 61282032.0 * eta2)
            + et6 * (83424402.0 - 123108426.0 * eta + 64828848.0 * eta2)
            + et4 * (783768.0 - 207204264.0 * eta + 166506060.0 * eta2)
            - 3024.0
            * (96.0 + 4268.0 * et2 + 4386.0 * et4 + 175.0 * et6)
            * (-5.0 + 2.0 * eta)
            * OTS
        )
        * v4
    ) / (45360.0 * OTS11)

    ans_v = CFdvdt * (POLdvdt0 + POLdvdt2 + (FACdvdt3 * RFTdvdt3) + POLdvdt4)

    return ans_e, ans_l, ans_v
