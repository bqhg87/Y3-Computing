import astropy.units as u
import numpy as np

G           =    (6.6726e-11 * (u.m**3 / (u.kg * u.s**2))).value
AU          =    (1.49597871e11 * u.m).to_value(u.m)
m_E         =    (5.9742e24 * u.kg).to_value(u.kg)
m_M         =    (7.35e22 * u.kg).to_value(u.kg)
m_Sun       =    (1.9884e30 * u.kg).to_value(u.kg)
m_Mercury   =    (3.3011e23 * u.kg).to_value(u.kg)
m_Venus     =    (4.8675e24 * u.kg).to_value(u.kg)
m_Mars      =    (6.4171e23 * u.kg).to_value(u.kg)
m_Jupiter   =    (1.8982e27 * u.kg).to_value(u.kg)
m_Saturn    =    (5.6834e26 * u.kg).to_value(u.kg)
a_EM, e_EM              =    (3.844e5 * u.km).to_value(u.m), 0.0549
a_Mercury, e_Mercury    =    (5.7909e7 * u.km).to_value(u.m), 0.2056
a_Venus, e_Venus        =    (1.0821e8 * u.km).to_value(u.m), 0.0068
a_Earth, e_Earth        =    AU, 0.0167
a_Mars, e_Mars          =    (2.2794e8 * u.km).to_value(u.m), 0.0934
a_Jupiter, e_Jupiter    =    (7.7833e8 * u.km).to_value(u.m), 0.0489
a_Saturn, e_Saturn      =    (1.4291e9 * u.km).to_value(u.m), 0.0547

VIS = {
    "Sun":           {"color": "gold",    "size": 250},
    "Mercury":       {"color": "silver",  "size": 40},
    "Venus":         {"color": "orange",  "size": 80},
    "Earth":         {"color": "blue",    "size": 100},
    "Moon":          {"color": "grey",    "size": 30},
    "Mars":          {"color": "red",     "size": 70},
    "Jupiter":       {"color": "sienna",  "size": 180},
    "Saturn":        {"color": "khaki",   "size": 150},
    "COM":           {"color": "white",   "size": 25,  "marker": "x"},
    "Test Particle": {"color": "cyan",  "size": 40,  "marker": "h"}
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# SIMULATION:

NAME_1, MASS_1 = "Earth", m_E
NAME_2, MASS_2 = "Moon", m_M
A_12, E_12 = a_EM, e_EM

T = 0.0
DT = 100.0
STEPS_PER_FRAME = 50
PAUSE_TIME = 0.001

SHOW_TEST_PARTICLE = True
USE_L2 = True
R_TEST = np.array([100.0 * a_EM, 200.0 * a_EM, 50 * a_EM])
V_TEST = np.array([500.0, -1000.0, -200.0])
M_TEST = 1000.0

SHOW_COM = False
SHOW_TELEMETRY = True
SHOW_LEGEND = False
SHOW_GRID = False