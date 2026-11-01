import astropy.units as u
import numpy as np
from visualiser import create_orbital_plot
import matplotlib.pyplot as plt
from collections import deque

G = (6.6726e-11 * (u.m**3 / (u.kg * u.s**2))).value
m_E = (5.9742e24 * u.kg).to_value(u.kg)
m_M = (7.35e22 * u.kg).to_value(u.kg)
d_EM = (3.844e5 * u.km).to_value(u.m)
a_EM = d_EM # semi-major axis
eccentricity_EM = 0.0549

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# ORBITS + VELOCITIES:

def calc_r_COM(m1, m2, d): # where r1 is the distance from m1 to COM etc.
    r1 = d * (m2 / (m1 + m2))
    r2 = d * (m1 / (m1 + m2))
    return r1, r2

def calc_orbital_state(a, eccentricity, m1, m2, true_anomaly_deg, G): # pos + vel for 2-body system at a given true anomaly
    nu = np.radians(true_anomaly_deg)
    
    p = a * (1.0 - eccentricity**2) # semi-latus rectum
    r_rel = p / (1.0 + eccentricity * np.cos(nu))

    r1_COM, r2_COM = calc_r_COM(m1, m2, r_rel)

    # unit vectors in xy plane
    unit_r = np.array([np.cos(nu), np.sin(nu), 0.0]) # aligned with r
    unit_t = np.array([-np.sin(nu), np.cos(nu), 0.0]) # tangential

    pos_1 = -r1_COM * unit_r
    pos_2 = r2_COM * unit_r

    mu = G * (m1 + m2)
    v_radial = np.sqrt(mu / p) * eccentricity * np.sin(nu)
    v_transverse = np.sqrt(mu / p) * (1.0 + eccentricity * np.cos(nu))

    vel_rel = v_radial * unit_r + v_transverse * unit_t
    vel_1 = -vel_rel * (m2 / (m1 + m2))
    vel_2 = vel_rel * (m1 / (m1 + m2))

    return pos_1, pos_2, vel_1, vel_2

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# GRAVITY:

def calc_test_potential(r_test, r_body, m_body, G):
    r_vec = r_test - r_body
    r_dist = np.linalg.norm(r_vec)
    if r_dist == 0.0:
        return 0.0 # prevent /0
    return -G * m_body / r_dist

def calc_test_accel(r_test, r_body, m_body, G):
    r_vec = r_test - r_body
    r_dist = np.linalg.norm(r_vec)
    if r_dist == 0.0:
        return np.zeros(3) # prevent /0
    return -G * m_body * r_vec / r_dist**3

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# SIMULATION:

# time system
dt = 100.0
steps_per_frame = 5
pause_time = 0.001
t = 0.0

mu_EM = G * (m_E + m_M) # gravitational parameter
n_EM = np.sqrt(mu_EM / a_EM**3) # mean orbital angular velocity (K3L)

def mean_to_true(M, e): # convert mean anomaly to true anomaly
    E = M + e * np.sin(M) # eccentric anomaly
    for _ in range(10): # newton-raphson convergence (small e)
        f = E - e * np.sin(E) - M # kepler eq
        df = 1.0 - e * np.cos(E)
        E -= f / df
    return 2.0 * np.arctan2(np.sqrt(1.0 + e) * np.sin(E / 2.0), np.sqrt(1.0 - e) * np.cos(E / 2.0))

def update_state(t):
    M = (n_EM * t) % (2.0 * np.pi) # mean anomaly (rad)
    nu_deg = np.degrees(mean_to_true(M, eccentricity_EM)) # true anomaly at t
    return calc_orbital_state(a_EM, eccentricity_EM, m_E, m_M, nu_deg, G)

def calc_accel(r, pos_E, pos_M): # earth + moon combined acceleration
    return calc_test_accel(r, pos_E, m_E, G) + calc_test_accel(r, pos_M, m_M, G)

pos_E, pos_M, vel_E, vel_M = update_state(t)
 
# scalar L2 distance in terms of EM dist from COM (circular approximation)
m_ratio = m_M / (m_E + m_M) # needs to be generalised
x_L2 = (1.0 - m_ratio) + (m_ratio / 3.0)**(1/3)
for _ in range(10): # newton-raphson convergence
   f = (1.0 - m_ratio) / (x_L2 + m_ratio)**2 + m_ratio / (x_L2 - 1.0 + m_ratio)**2 - x_L2
   df = -2.0 * (1.0 - m_ratio) / (x_L2 + m_ratio)**3 - 2.0 * m_ratio / (x_L2 - 1.0 + m_ratio)**3 - 1.0
   x_L2 -= f / df
 
# test particle at L2:
r_test = x_L2 * (pos_M - pos_E)
v_test = x_L2 * (vel_M - vel_E)
m_test = 1000.0 # (irrelevant)
 
acc = calc_accel(r_test, pos_E, pos_M)

system = {
    "Earth": {"pos": pos_E, "color": "blue", "size": 100},
    "Moon": {"pos": pos_M, "color": "grey", "size": 50},
    "COM": {"pos": np.array([0.0, 0.0, 0.0]), "color": "white", "size": 25, "marker": "x"},
    "Test Particle": {"pos": r_test, "color": "yellow", "size": 40, "marker": "*"},
}

fig, ax, scatters = create_orbital_plot(system)

limit = a_EM * 1.5
ax.set_xlim([-limit, limit]), ax.set_ylim([-limit, limit]), ax.set_zlim([-limit, limit])

max_points = 10000
earth_pts = deque(maxlen=max_points)
moon_pts = deque(maxlen=max_points)
test_pts = deque(maxlen=max_points)

earth_trail, = ax.plot([], [], [], color="blue", alpha=0.4, linewidth=1.5)
moon_trail, = ax.plot([], [], [], color="grey", alpha=0.4, linewidth=1.5)
test_trail, = ax.plot([], [], [], color="yellow", alpha=0.6, linewidth=1.5)

info_text = ax.text2D(0.05, 0.95, "", transform=ax.transAxes, color="white", fontsize=9, family="monospace", va="top")

try:
    while plt.fignum_exists(fig.number):
        for _ in range(steps_per_frame): # repeat 4th-order runge-kutta for each step
            h = dt
            t_n = t
            
            # starting point (t_n)
            pos_En, pos_Mn, _, _ = update_state(t_n) # don't care abt velocities (use _)
            ddot_n = calc_accel(r_test, pos_En, pos_Mn) # acceleration at t_n
            
            # substep 1 (t_n + h/2)
            t_1 = t_n + 0.5 * h
            pos_E1, pos_M1, _, _ = update_state(t_1)
            z1 = r_test + 0.5 * h * v_test
            zdot_1 = v_test + 0.5 * h * ddot_n
            ddot_1 = calc_accel(z1, pos_E1, pos_M1)
            
            # substep 2 (t_n + h/2)
            t_2 = t_n + 0.5 * h
            pos_E2, pos_M2, _, _ = update_state(t_2)
            z2 = r_test + 0.5 * h * zdot_1
            zdot_2 = v_test + 0.5 * h * ddot_1
            ddot_2 = calc_accel(z2, pos_E2, pos_M2)
            
            # substep 3 (t_n + h)
            t_3 = t_n + h
            pos_E3, pos_M3, _, _ = update_state(t_3)
            z3 = r_test + h * zdot_2
            zdot_3 = v_test + h * ddot_2
            ddot_3 = calc_accel(z3, pos_E3, pos_M3)
            
            # updates
            r_test += (h / 6.0) * (v_test + 2.0 * zdot_1 + 2.0 * zdot_2 + zdot_3)
            v_test += (h / 6.0) * (ddot_n + 2.0 * ddot_1 + 2.0 * ddot_2 + ddot_3)
            t += h
            pos_E, pos_M, vel_E, vel_M = update_state(t)
            acc = calc_accel(r_test, pos_E, pos_M)

        # info display
        info_text.set_text(
            f"Time: {t / 86400.0:.2f} days\n"
            f"Test Position: [{r_test[0]:.2e}, {r_test[1]:.2e}, {r_test[2]:.2e}] m\n"
            f"Gravitational Force: {(np.linalg.norm(acc) * m_test):.4e} N"
        )

        # add new trail point
        earth_pts.append(pos_E), moon_pts.append(pos_M), test_pts.append(r_test.copy())

        # update particle positions
        scatters["Earth"]._offsets3d = (np.array([pos_E[0]]), np.array([pos_E[1]]), np.array([pos_E[2]]))
        scatters["Moon"]._offsets3d = (np.array([pos_M[0]]), np.array([pos_M[1]]), np.array([pos_M[2]]))
        scatters["Test Particle"]._offsets3d = (np.array([r_test[0]]), np.array([r_test[1]]), np.array([r_test[2]]))

        # update trails
        earth_trail.set_data_3d(*np.array(earth_pts).T)
        moon_trail.set_data_3d(*np.array(moon_pts).T)
        test_trail.set_data_3d(*np.array(test_pts).T)

        plt.pause(pause_time)

except KeyboardInterrupt:
    pass