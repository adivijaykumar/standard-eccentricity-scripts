from systems_of_eqns import evolution_equations_2PN, MTSUN_SI
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import numpy as np
from utils import freq, f0_from_a0, time_to_coalescence_circular

e0 = 0.9
a0_in_AU = 1e-4
m1_src = 35  # IN SOLAR MASSES
m2_src = 20  # IN SOLAR MASSES
redshift = 0.5
extraction_f22 = [20, 10, 5]
extraction_M_times_f22 = [1000, 500, 300]


m1, m2 = m1_src * (1 + redshift), m2_src * (1 + redshift)
MTOT = m1 + m2
eta = m1 * m2 / (m1 + m2) ** 2
f0 = f0_from_a0(a0_in_AU, m1, m2)
v0 = (MTOT * MTSUN_SI * np.pi * f0) ** (1 / 3)
tmax = time_to_coalescence_circular(f0, m1, m2)  # APPROXIMATE USING CIRCULAR CASE
if tmax < 1e6:
    tmax = 1e6  # INTEGRATE AT LEAST FOR 1e6 seconds

kwargs = dict(
    t_span=[0, tmax], y0=[e0, 0, v0], args=[MTOT, eta], t_eval=np.arange(0, tmax, 1 / 8)
)
sol_2PN = solve_ivp(evolution_equations_2PN, **kwargs)
e_vs_f = interp1d(
    freq(sol_2PN.y[2], MTOT),
    sol_2PN.y[0],
)

print("\n#############\n")
for f22 in extraction_f22:
    e_at_f22 = e_vs_f(f22)
    print(f"Eccentricity at f22 = {f22} Hz is {e_at_f22:.3f}")

print("\n#############\n")

for M_times_f22 in extraction_M_times_f22:
    e_at_M_times_f22 = e_vs_f(M_times_f22 / MTOT)
    print(f"Eccentricity at M*f22 = {M_times_f22} Hz MSun is {e_at_M_times_f22:.3f}")

print("\n#############\n")
