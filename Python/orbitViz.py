import numpy as np 
from vpython import *

def compute_R_ECI_2_ECEF(alpha):
    R = np.arr([[np.cos(alpha), np.sin(alpha), 0],
                [-np.sin(alpha), np.cos(alpha), 0],
                [0, 0, 1]])
    return R

def LLA_2_ECEF(lat, lon, alt):
    x_ECEF = (6378*1000 + alt)*np.cos(lat)*np.cos(lon)
    y_ECEF = (6378*1000 + alt)*np.cos(lat)*np.sin(lon)
    z_ECEF = (6378*1000 + alt)*np.sin(lat)

    r_ECEF = np.arr([x_ECEF, y_ECEF, z_ECEF])
    return r_ECEF

if __name__ == 'main':
    data = np.loadtxt("../C++/orbs.csv")

    # Time array
    ts = data[:, 0]

    # State of Earth w.r.t. Sun
    S2E_px = data[:, 1]
    S2E_py = data[:, 2]
    S2E_pz = data[:, 3]
    S2E_vx = data[:, 4]
    S2E_vy = data[:, 5]
    S2E_vz = data[:, 6]

    # State of actual satellite in ECI frame
    Sat_ECI_px = data[:, 7]
    Sat_ECI_py = data[:, 8]
    Sat_ECI_pz = data[:, 9]
    Sat_ECI_vx = data[:, 10]
    Sat_ECI_vy = data[:, 11]
    Sat_ECI_vz = data[:, 12]

    # State of reference satellite in ECI frame
    refSat_ECI_px = data[:, 13]
    refSat_ECI_py = data[:, 14]
    refSat_ECI_pz = data[:, 15]
    refSat_ECI_vx = data[:, 16]
    refSat_ECI_vy = data[:, 17]
    refSat_ECI_vz = data[:, 18]

    # Actual satellite control inputs in ECI frame
    u_ECI_x = data[:, 19]
    u_ECI_y = data[:, 20]
    u_ECI_z = data[:, 21]

    # Actual satellite control inputs in LVLH frame
    u_LVLH_x = data[:, 22]
    u_LVLH_y = data[:, 23]
    u_LVLH_z = data[:, 24]

    # Actual satellite LVLH state
    Sat_LVLH_px = data[:, 25]
    Sat_LVLH_py = data[:, 26]
    Sat_LVLH_pz = data[:, 27]
    Sat_LVLH_vx = data[:, 28]
    Sat_LVLH_vy = data[:, 29]
    Sat_LVLH_vz = data[:, 30]

    # Actual satellite LLA
    Sat_lat = data[:, 31]
    Sat_lon = data[:, 32]
    Sat_alt = data[:, 33]

    # Alpha (angle between ECI and ECEF frames)
    alphas = data[:, 34]

    #helioCanvas = canvas(title="Heliocentric View", width=400, height=300, background=color.black)
    geoCanvas = canvas(title="Geocentric View", width=400, height=300, background=color.black)
    #lvlhCanvas = canvas(title="Local Vertical Line Heading View", width=400, height=300, background=color.gray)

    # Set up geoCanvas

    R = 6378*1000 # m
    earth = sphere(canvas=geoCanvas, pos=vector(0,0,0), radius=R, texture=textures.earth, shininess=0)
    earth.rotate(origin=vector(0,0,0), axis=vector(0,0,1))
    geoCanvas.lights = []
    geoCanvas.lights = geoCanvas.lights + [distant_light(direction=vector(-cos(23.5*pi/180),-sin(23.5*pi/180),0), color=vector(0.9,0.9,0.9))]

    lat = 38.836587*(pi/180)
    lon = -77.196691*(pi/180)
    alt = R
    r_target_ECEF = LLA_2_ECEF(lat, lon, alt)
    R_ECI_2_ECEF = compute_R_ECI_2_ECEF(alphas[0])
    R_ECEF_2_ECI = np.transpose(R_ECI_2_ECEF)
    r_target_ECI = np.dot(R_ECEF_2_ECI, r_target_ECEF)
    ball = sphere(pos=vector(r_target_ECI[0], r_target_ECI[1], r_target_ECI[2]), radius = 6378*30, color=color.red, make_trail=False)
    ball.rotate(origin=vector(0,0,0), axis=vector(0,0,1))

    w = 2*pi/86164 # rotation rate of Earth, rad/s
    count = 0
    t = ts[count]
    dt = ts[2] - ts[1]

    while t<ts[-1]:
        rate(100)
        count = count + 1
        earth.rotate(origin=vector(0,0,0), axis=vector(0,0,1), angle = w*dt)
        t = ts[count]
        alpha = alphas[count]
        R_ECI_2_ECEF = compute_R_ECI_2_ECEF(alpha)
        R_ECEF_2_ECI = np.transpose(R_ECI_2_ECEF)
        r_target_ECI = np.dot(R_ECEF_2_ECI, r_target_ECEF)
        ball.pos = vector(r_target_ECI[0], r_target_ECI[1], r_target_ECI[2])









