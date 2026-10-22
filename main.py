import astropy.units as u
import numpy as np
from visualiser import create_orbital_plot
import matplotlib.pyplot as plt

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
    
    p = a * (1.0 - eccentricity**2)  # semi-latus rectum
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

true_anomaly_deg = 0.0 # 0 perigee, 180 apogee

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# VISUALISER:

# UPDATE TO REMOVE EARTH / MOON DEPENDENCE
# CUSTOM PARTICLES?
# CHANGE REFERENCE POINT / Origin

a_vis = a_EM
eccentricity_vis = eccentricity_EM
deg_vis = 0.0
m1_vis = m_E
m2_vis = m_M
G_vis = G
trail_vis = 20.0 # length of trail in seconds (0 = no trail)
increment_time_vis = 0.01 # how many seconds to increment the orbit
increment_amt_vis = 0.5 # how many degrees per increment time

def update_state():
    return calc_orbital_state(a_vis, eccentricity_vis, m1_vis, m2_vis, deg_vis, G_vis)

pos_E, pos_M, vel_E, vel_M = update_state()

system = {
    "Earth": {"pos": pos_E, "color": "blue", "size": 100},
    "Moon": {"pos": pos_M, "color": "grey", "size": 50},
    "COM": {"pos": np.array([0.0, 0.0, 0.0]), "color": "white", "size": 25, "marker": "x"},
}

fig, ax, scatters = create_orbital_plot(system)

# debug: locking axis limits prevents visualisation from collapsing
limit = a_vis * 1.2
ax.set_xlim([-limit, limit]), ax.set_ylim([-limit, limit]), ax.set_zlim([-limit, limit])

# trail:
earth_x, earth_y, earth_z = [], [], []
moon_x, moon_y, moon_z = [], [], []

earth_trail, = ax.plot([], [], [], color="blue", alpha=0.4, linewidth=1.5)
moon_trail, = ax.plot([], [], [], color="grey", alpha=0.4, linewidth=1.5)

max_points = int(trail_vis / increment_time_vis)

try:
    while True:
        pos_E, pos_M, vel_E, vel_M = update_state()

        # add current positions
        earth_x.append(pos_E[0]), earth_y.append(pos_E[1]), earth_z.append(pos_E[2])
        moon_x.append(pos_M[0]), moon_y.append(pos_M[1]), moon_z.append(pos_M[2])

        # update positions
        scatters["Earth"]._offsets3d = (np.array([pos_E[0]]), np.array([pos_E[1]]), np.array([pos_E[2]]))
        scatters["Moon"]._offsets3d = (np.array([pos_M[0]]), np.array([pos_M[1]]), np.array([pos_M[2]]))

        # update trails
        earth_trail.set_data_3d(earth_x[-max_points:], earth_y[-max_points:], earth_z[-max_points:])
        moon_trail.set_data_3d(moon_x[-max_points:], moon_y[-max_points:], moon_z[-max_points:])

        fig.canvas.draw_idle()
        fig.canvas.flush_events()
        plt.pause(increment_time_vis) # animation update period

        deg_vis = (deg_vis + increment_amt_vis) % 360.0 # loops

except KeyboardInterrupt:
    print("Simulation ended")