import astropy.units as u
import numpy as np

# CONSTANTS:

G = (6.6726e-11 * (u.m**3 / (u.kg * u.s**2))).value
m_E = (5.9742e24 * u.kg).to_value(u.kg)
m_M = (7.35e22 * u.kg).to_value(u.kg)
d_EM = (3.844e5 * u.km).to_value(u.m)

a_EM = d_EM # semi-major axis
eccentricity_EM = 0.0549

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def calc_r_COM(m1, m2, d): # where r1 is the distance from m1 to COM etc.
    r1 = d * (m2 / (m1 + m2))
    r2 = d * (m1 / (m1 + m2))
    return r1, r2

def calc_periapsis(a, eccentricity):
    return a * (1.0 - eccentricity)

r_perigee = calc_periapsis(a_EM, eccentricity_EM) # perigee = closest E-M approach
r_E, r_M = calc_r_COM(m_E, m_M, r_perigee) # COM distances with elliptical orbit
r_E_circular, r_M_circular = calc_r_COM(m_E, m_M, d_EM) # COM "" circular orbit

pos_E = np.array([-r_E, 0.0, 0.0])
pos_M = np.array([r_M, 0.0, 0.0])

print("Earth Position COM frame:", pos_E)
print("Moon Position COM frame: ", pos_M)

print("Old Earth COM Position:", np.array([-r_E_circular, 0.0, 0.0]))
print("Old Moon COM Position: ", np.array([-r_M_circular, 0.0, 0.0]))