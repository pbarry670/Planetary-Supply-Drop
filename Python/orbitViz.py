import numpy as np
import pandas as pd
from vpython import canvas, sphere, color, vector, textures, distant_light, local_light, arrow, rate

def compute_R_ECI_2_ECEF(alpha):
    R = np.array([[np.cos(alpha), np.sin(alpha), 0],
                [-np.sin(alpha), np.cos(alpha), 0],
                [0, 0, 1]])
    return R

def LLA_2_ECEF(lat, lon, alt):
    x_ECEF = (6378*1000 + alt)*np.cos(lat)*np.cos(lon)
    y_ECEF = (6378*1000 + alt)*np.cos(lat)*np.sin(lon)
    z_ECEF = (6378*1000 + alt)*np.sin(lat)

    r_ECEF = np.array([x_ECEF, y_ECEF, z_ECEF])
    return r_ECEF

def ECI_2_VPython(x, y, z):
    newx = y
    newy = z
    newz = x

    return newx, newy, newz

def LVLH_2_VPython(x, y, z):
    newx = -y
    newy = z
    newz = x

    return newx, newy, newz

# Read orbital data .csv file
df = pd.read_csv("../C++/results/orbs.csv", header=None, names=['t',
                                                        'S2Epx','S2Epy','S2Epz','S2Evx','S2Evy','S2Evz',
                                                        'satECIpx','satECIpy','satECIpz','satECIvx','satECIvy','satECIvz',
                                                        'refECIpx', 'refECIpy','refECIpz','refECIvx','refECIvy','refECIvz',
                                                        'uECIx','uECIy','uECIz',
                                                        'uLVLHx','uLVLHy','uLVLHz',
                                                        'satLVLHpx','satLVLHpy','satLVLHpz','satLVLHvx','satLVLHvy','satLVLHvz',
                                                        'satlat', 'satlon', 'satalt',
                                                        'alpha'])
numrows, numcols = df.shape

# Time array
ts = df['t']
ts = np.array(ts)

# State of Earth w.r.t. Sun
S2E_px = df['S2Epx']
S2E_py = df['S2Epy']
S2E_pz = df['S2Epz']
S2E_vx = df['S2Evx']
S2E_vy = df['S2Evy']
S2E_vz = df['S2Evz']

# State of actual satellite in ECI frame
Sat_ECI_px = df['satECIpx']
Sat_ECI_py = df['satECIpy']
Sat_ECI_pz = df['satECIpz']
Sat_ECI_vx = df['satECIvx']
Sat_ECI_vy = df['satECIvy']
Sat_ECI_vz = df['satECIvz']

# State of reference satellite in ECI frame
refSat_ECI_px = df['refECIpx']
refSat_ECI_py = df['refECIpy']
refSat_ECI_pz = df['refECIpz']
refSat_ECI_vx = df['refECIvx']
refSat_ECI_vy = df['refECIvy']
refSat_ECI_vz = df['refECIvz']

# Actual satellite control inputs in ECI frame
u_ECI_x = df['uECIx']
u_ECI_y = df['uECIy']
u_ECI_z = df['uECIz']

# Actual satellite control inputs in LVLH frame
u_LVLH_x = df['uLVLHx']
u_LVLH_y = df['uLVLHy']
u_LVLH_z = df['uLVLHz']

# Actual satellite LVLH state
Sat_LVLH_px = df['satLVLHpx']
Sat_LVLH_py = df['satLVLHpy']
Sat_LVLH_pz = df['satLVLHpz']
Sat_LVLH_vx = df['satLVLHvx']
Sat_LVLH_vy = df['satLVLHvy']
Sat_LVLH_vz = df['satLVLHvz']

# Actual satellite LLA
Sat_lat = df['satlat']
Sat_lon = df['satlon']
Sat_alt = df['satalt']

# Alpha (angle between ECI and ECEF frames)
alphas = df['alpha']

geoCanvas = canvas(title="Geocentric View", width=400, height=300, background=color.black)
helioCanvas = canvas(title="Heliocentric View", width=400, height=300, background=color.black)
lvlhCanvas = canvas(title="Local Vertical Line Heading View", width=400, height=300, background=vector(0.5,0.5,0.5))

# Set up geoCanvas
earth = sphere(canvas=geoCanvas, pos=vector(0,0,0), radius=6378*1000, texture=textures.earth, shininess=0)
geoCanvas.lights = []
distant_light(canvas=geoCanvas, direction=vector(-np.cos(23.5*np.pi/180),-np.sin(23.5*np.pi/180),0), color=vector(0.9,0.9,0.9))

lat = 38.836587*(np.pi/180)
lon = -77.196691*(np.pi/180)
alt = 0
r_target_ECEF = LLA_2_ECEF(lat, lon, alt)
R_ECI_2_ECEF = compute_R_ECI_2_ECEF(alphas[0])
R_ECEF_2_ECI = np.transpose(R_ECI_2_ECEF)
r_target_ECI = np.dot(R_ECEF_2_ECI, r_target_ECEF)
targetx, targety, targetz = ECI_2_VPython(r_target_ECI[0], r_target_ECI[1], r_target_ECI[2])
target = sphere(canvas=geoCanvas, pos=vector(targetx, targety, targetz), radius = 6378*30, color=color.red, make_trail=False)

satECIx, satECIy, satECIz = ECI_2_VPython(Sat_ECI_px[0], Sat_ECI_py[0], Sat_ECI_pz[0])
sat = sphere(canvas=geoCanvas, pos=vector(satECIx, satECIy, satECIz), radius = 6378*50, color=color.magenta, make_trail = False)

refSatECIx, refSatECIy, refSatECIz = ECI_2_VPython(refSat_ECI_px[0], refSat_ECI_py[0], refSat_ECI_pz[0])
refSat = sphere(canvas=geoCanvas, pos=vector(refSatECIx, refSatECIy, refSatECIz), radius = 6378*10, color=color.green, make_trail = True)

# Set up helioCanvas
distance_scalefactor = 100
sun = sphere(canvas=helioCanvas, pos=vector(0,0,0), radius = 696340000, shininess=0, color=color.yellow)
helioEx, helioEy, helioEz = ECI_2_VPython(S2E_px[0], S2E_py[0], S2E_pz[0])
helioEarth = sphere(canvas=helioCanvas, pos=vector(helioEx/distance_scalefactor, helioEy/distance_scalefactor, helioEz/distance_scalefactor), radius = 6378*1000, shininess=0, texture=textures.earth, make_trail=True, trail_radius=6378*100)
helioCanvas.lights = []
local_light(canvas=helioCanvas, pos=vector(696340000,696340000,0), color=vector(0.9,0.9,0.9))
local_light(canvas=helioCanvas, pos=vector(696340000,0,696340000), color=vector(0.9,0.9,0.9))
local_light(canvas=helioCanvas, pos=vector(0,696340000,696340000), color=vector(0.9,0.9,0.9))
local_light(canvas=helioCanvas, pos=vector(-696340000,696340000,0), color=vector(0.9,0.9,0.9))
local_light(canvas=helioCanvas, pos=vector(-696340000,0,696340000), color=vector(0.9,0.9,0.9))
local_light(canvas=helioCanvas, pos=vector(0,-696340000,696340000), color=vector(0.9,0.9,0.9))
local_light(canvas=helioCanvas, pos=vector(696340000,-696340000,0), color=vector(0.9,0.9,0.9))
local_light(canvas=helioCanvas, pos=vector(696340000,0,-696340000), color=vector(0.9,0.9,0.9))
local_light(canvas=helioCanvas, pos=vector(0,696340000,-696340000), color=vector(0.9,0.9,0.9))
local_light(canvas=helioCanvas, pos=vector(-696340000,-696340000,0), color=vector(0.9,0.9,0.9))
local_light(canvas=helioCanvas, pos=vector(-696340000,0,-696340000), color=vector(0.9,0.9,0.9))
local_light(canvas=helioCanvas, pos=vector(0,-696340000,-696340000), color=vector(0.9,0.9,0.9))
local_light(canvas=helioCanvas, pos=vector(696340000,696340000,696340000), color=vector(0.9,0.9,0.9))
local_light(canvas=helioCanvas, pos=vector(696340000,696340000,-696340000), color=vector(0.9,0.9,0.9))
local_light(canvas=helioCanvas, pos=vector(696340000,-696340000,-696340000), color=vector(0.9,0.9,0.9))
local_light(canvas=helioCanvas, pos=vector(-696340000,-696340000,-696340000), color=vector(0.9,0.9,0.9))
local_light(canvas=helioCanvas, pos=vector(-696340000,696340000,696340000), color=vector(0.9,0.9,0.9))
local_light(canvas=helioCanvas, pos=vector(-696340000,-696340000,696340000), color=vector(0.9,0.9,0.9))
helioSat = sphere(canvas=helioCanvas, pos=vector(helioEx/distance_scalefactor + satECIx, helioEy/distance_scalefactor + satECIy, helioEz/distance_scalefactor + satECIz), radius = 6378*100, shininess=0, color=color.magenta, make_trail=True, trail_radius=6378*50, retain=200)

# Set up lvlhCanvas
origin = sphere(canvas=lvlhCanvas, pos=vector(0,0,0), radius=250, shininess=0, color=color.green)
relposx, relposy, relposz = LVLH_2_VPython(Sat_LVLH_px[0], Sat_LVLH_py[0], Sat_LVLH_pz[0])
relctrlx, relctrly, relctrlz = LVLH_2_VPython(u_LVLH_x[0], u_LVLH_y[0], u_LVLH_z[0])
relativeSat = sphere(canvas=lvlhCanvas, pos=vector(relposx, relposy, relposz), radius=250, shininess=0, color=color.magenta)
lvlhControl = arrow(canvas=lvlhCanvas, pos=vector(relposx, relposy, relposz), axis=vector(relctrlx, relctrly, relctrlz), color = color.cyan)

w = 2*np.pi/86164 # rotation rate of Earth, rad/s
count = 0
t = ts[count]
dt = ts[1] - ts[0]


hasReachedFinalTime = False
while not hasReachedFinalTime:
    rate(500)
    count = count + 1

    # Update geoCanvas
    earth.rotate(origin=vector(0,0,0), axis=vector(0,1,0), angle = w*dt)
    t = ts[count]
    alpha = alphas[count]

    R_ECI_2_ECEF = compute_R_ECI_2_ECEF(alpha)
    R_ECEF_2_ECI = np.transpose(R_ECI_2_ECEF)
    r_target_ECI = np.dot(R_ECEF_2_ECI, r_target_ECEF)
    targetx, targety, targetz = ECI_2_VPython(r_target_ECI[0], r_target_ECI[1], r_target_ECI[2])
    target.pos = vector(targetx, targety, targetz)

    satECIx, satECIy, satECIz = ECI_2_VPython(Sat_ECI_px[count], Sat_ECI_py[count], Sat_ECI_pz[count])
    sat.pos = vector(satECIx, satECIy, satECIz)

    refSatECIx, refSatECIy, refSatECIz = ECI_2_VPython(refSat_ECI_px[count], refSat_ECI_py[count], refSat_ECI_pz[count])
    refSat.pos = vector(refSatECIx, refSatECIy, refSatECIz)

    # Update helioCanvas
    helioEx, helioEy, helioEz = ECI_2_VPython(S2E_px[count], S2E_py[count], S2E_pz[count])
    helioEarth.pos = vector(helioEx/distance_scalefactor, helioEy/distance_scalefactor, helioEz/distance_scalefactor)
    helioSat.pos = vector(helioEx/distance_scalefactor +satECIx, helioEy/distance_scalefactor +satECIy, helioEz/distance_scalefactor +satECIz)
    helioEarth.rotate(origin=vector(helioEx/distance_scalefactor + satECIx, helioEy/distance_scalefactor + satECIy, helioEz/distance_scalefactor + satECIz), axis=vector(0,1,0), angle = w*dt)

    # Update lvlhCanvas
    relposx, relposy, relposz = LVLH_2_VPython(Sat_LVLH_px[count], Sat_LVLH_py[count], Sat_LVLH_pz[count])
    relctrlx, relctrly, relctrlz = LVLH_2_VPython(u_LVLH_x[count], u_LVLH_y[count], u_LVLH_z[count])
    relativeSat.pos = vector(relposx, relposy, relposz)
    lvlhControl.pos = vector(relposx, relposy, relposz)
    lvlhControl.axis = vector(relctrlx, relctrly, relctrlz)
    lvlhControl.headwidth = lvlhControl.shaftwidth

    if t == ts[-1]:
        hasReachedFinalTime = True

print("Orbital visualization complete! \n")










