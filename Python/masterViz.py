import numpy as np
import pandas as pd
from vpython import canvas, sphere, arrow, color, vector, distant_light, textures, local_light, graph, gcurve, box, rate
import matplotlib.pyplot as plt

def visualizeOrbits(filepath, animationRate, latTarget, lonTarget):

    # Read orbital data .csv file
    df = pd.read_csv(filepath, header=None, names=['t',
                                                    'S2Epx','S2Epy','S2Epz','S2Evx','S2Evy','S2Evz',
                                                    'satECIpx','satECIpy','satECIpz','satECIvx','satECIvy','satECIvz',
                                                    'refECIpx', 'refECIpy','refECIpz','refECIvx','refECIvy','refECIvz',
                                                    'uECIx','uECIy','uECIz',
                                                    'uLVLHx','uLVLHy','uLVLHz',
                                                    'satLVLHpx','satLVLHpy','satLVLHpz','satLVLHvx','satLVLHvy','satLVLHvz',
                                                    'satlat', 'satlon', 'satalt',
                                                    'alpha'])

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

    lat = latTarget
    lon = lonTarget
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
        rate(animationRate)
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



def visualizeBallisticDescent(filepath, animationRate):

    # Read ballistic descent data from a .csv file
    df_bd = pd.read_csv(filepath, header=None, names = ['t', 'V', 'gamma', 'h', 'totalDTheta'])

    ts_b = df_bd['t']
    Vs = df_bd['V']
    gammas = df_bd['gamma']
    hs = df_bd['h']
    cumulativeDTheta = df_bd['totalDTheta']

    ts_b = np.array(ts_b)
    Vs = np.array(Vs)
    gammas = np.array(gammas)
    hs = np.array(hs)
    cumulativeDTheta = np.array(cumulativeDTheta)

    dps = np.empty((ts_b.shape[0], 2)) # allocate space for holding positions, x y
    dps[0, 0] = 0
    dps[0, 1] = 0
    dt = 0.05

    for i in range(ts_b.shape[0] - 1):
        Vy = -1*Vs[i]*np.cos(gammas[i])
        Vx = -1*Vs[i]*np.sin(gammas[i])
        dpx = Vx*dt
        dpy = Vy*dt

        dps[i+1, 0] = dpx
        dps[i+1, 1] = dpy

    positions = np.cumsum(dps, axis=0)
    positions[:, 0] = positions[:, 0] - positions[-1, 0] + hs[-1]
    positions[:, 1] = positions[:, 1] - positions[-1, 1]

    horizontal_offset = -1*np.mean(positions[:,1])
    vertical_offset = np.mean(positions[:,0])

    # Set up ballistic descent canvas
    scene = canvas(title="Ballistic Descent, Local Level View", width=1000, height=400, background=color.cyan)
    scene.lights = []
    distant_light(canvas=scene, direction=vector(0,0,1), color=vector(0.9,0.9,0.9))
    distant_light(canvas=scene, direction=vector(0,0,-1), color=vector(0.9,0.9,0.9))
    capx, capy = LL_2_VPython(positions[0,0], positions[0,1])
    capsule_b = sphere(canvas=scene, pos=vector(capx-horizontal_offset, capy-vertical_offset, 0), radius=20000, shininess=0, color=color.blue, make_trail=True, trail_color=color.blue)
    target_b = sphere(canvas=scene, pos=vector(-horizontal_offset, -vertical_offset, 0), radius=50000, shininess=0, color=color.red)
    ground = box(canvas=scene, pos=vector(-horizontal_offset, -vertical_offset, 0), length=3000000, height=50, width=3000000, color=color.green)

    # Set up ballistic descent velocity graph
    velGraph = graph(title = 'Capsule Speed', xtitle='Global Time (s)', ytitle='Speed (m/s)')
    Vcap = gcurve(graph=velGraph, color=color.blue)
    Vcap.plot(ts_b[0], Vs[0])

    # Set up ballistic descent flight path angle graph
    fpaGraph = graph(title = 'Capsule Flight Path Angle', xtitle='Global Time (s)', ytitle='Flight Path Angle (deg)')
    fpa = gcurve(graph=fpaGraph, color=color.red)
    fpa.plot(ts_b[0], gammas[0]*(180/np.pi))

    # Set up ballistic descent altitude graph
    altGraph = graph(title = 'Capsule Altitude', xtitle='Global Time (s)', ytitle='Altitude (m)')
    alt = gcurve(graph=altGraph, color=color.purple)
    alt.plot(ts_b[0], hs[0])

    hasReachedFinalBallisticTime = False
    count = 0
    while not hasReachedFinalBallisticTime:
        rate(animationRate)
        count = count + 1
        t = ts_b[count]
        capx, capy = LL_2_VPython(positions[count, 0], positions[count,1])
        capsule_b.pos = vector(capx-horizontal_offset, capy-vertical_offset, 0)

        Vcap.plot(ts_b[count], Vs[count])
        fpa.plot(ts_b[count], gammas[count]*(180/np.pi))
        alt.plot(ts_b[count], hs[count])

        if t == ts_b[-1]:
            hasReachedFinalBallisticTime = True


def visualizeChuteDescent(filepath, animationRate):

    df_c = pd.read_csv(filepath, header=None)
    ts_c = np.array(df_c[0])
    vxs = np.array(df_c[1])
    pxs = np.array(df_c[2])
    vys = np.array(df_c[3])
    pys = np.array(df_c[4])
    vzs = np.array(df_c[5])
    pzs = np.array(df_c[6])

    horizontal_offset = -1*pys.mean()
    vertical_offset = pxs.mean()

    # Set up chute descent local-level canvas
    cdCanvas = canvas(title="Chute Descent, Local Level View", width=450, height=450, background=color.cyan)
    cdTopCanvas = canvas(title="Chute Descent, Birds' Eye View", width=450, height=450, background=color.green)

    # Set up chute descent local-level canvas
    cdCanvas.lights = []
    distant_light(canvas=cdCanvas, direction=vector(0,0,1), color=vector(0.9,0.9,0.9))
    distant_light(canvas=cdCanvas, direction=vector(0,0,-1), color=vector(0.9,0.9,0.9))
    capx, capy = LL_2_VPython(pxs[0], pys[0])
    capsule_c = sphere(canvas=cdCanvas, pos=vector(capx-horizontal_offset, capy-vertical_offset, pzs[0]), radius=100, shininess=0, color=color.blue, make_trail=True, trail_radius=50)
    target_c = sphere(canvas=cdCanvas, pos=vector(-horizontal_offset,-vertical_offset,0), radius=50, color=color.red)
    ground = box(canvas=cdCanvas, pos=vector(-horizontal_offset, -vertical_offset, 0), length=4000, height=10, width=4000, color=color.green)

    # Set up chute descent birds-eye canvas
    cdTopCanvas.lights = []
    distant_light(canvas=cdTopCanvas, direction=vector(0,0,1), color=vector(0.9,0.9,0.9))
    capsule_ctop = sphere(canvas=cdTopCanvas, pos=vector(pys[0], pzs[0], pxs[0]), radius=100, shininess=0, color=color.blue, make_trail=True, trail_radius=10)
    target_ctop = sphere(canvas=cdTopCanvas, pos=vector(0,0,0), radius=50, color=color.red)

    # Set up chute descent position graph
    posGraph = graph(title = 'Lander Position', xtitle='Global Time (s)', ytitle='Position (m)')
    xpos = gcurve(graph=posGraph, color=color.blue, label='X')
    ypos = gcurve(graph=posGraph, color=color.green, label='Y')
    zpos = gcurve(graph=posGraph, color=color.magenta, label='Z')
    xpos.plot(ts_c[0], vxs[0])
    ypos.plot(ts_c[0], vys[0])
    zpos.plot(ts_c[0], vzs[0])

    # Set up chute descent velocity graph
    velGraph = graph(title = 'Lander Velocity', xtitle='Global Time (s)', ytitle='Velocity (m)')
    xvel = gcurve(graph=velGraph, color=color.blue, label='X')
    yvel = gcurve(graph=velGraph, color=color.green, label='Y')
    zvel = gcurve(graph=velGraph, color=color.magenta, label='Z')
    xvel.plot(ts_c[0], vxs[0])
    yvel.plot(ts_c[0], vys[0])
    zvel.plot(ts_c[0], vzs[0])

    hasReachedFinalChuteTime = False
    count = 0
    while not hasReachedFinalChuteTime:
        rate(animationRate)
        count = count+1
        t = ts_c[count]

        capx, capy = LL_2_VPython(pxs[count], pys[count])
        capsule_c.pos = vector(capx-horizontal_offset, capy-vertical_offset, pzs[count])

        capsule_ctop.pos = vector(pys[count], pzs[count], pxs[count])

        xpos.plot(ts_c[count], pxs[count])
        ypos.plot(ts_c[count], pys[count])
        zpos.plot(ts_c[count], pzs[count])

        xvel.plot(ts_c[count], vxs[count])
        yvel.plot(ts_c[count], vys[count])
        zvel.plot(ts_c[count], vzs[count])        

        if t == ts_c[-1]:
            hasReachedFinalChuteTime = True


def visualizePoweredDescent(filepath, animationRate):

    df_p = pd.read_csv(filepath, header=None)
    ts_p = np.array(df_p[0])
    pxs = np.array(df_p[1])
    pys = np.array(df_p[2])
    pzs = np.array(df_p[3])
    vxs = np.array(df_p[4])
    vys = np.array(df_p[5])
    vzs = np.array(df_p[6])
    ms = np.array(df_p[7])
    Txs = np.array(df_p[8])
    Tys = np.array(df_p[9])
    Tzs = np.array(df_p[10])

    horizontal_offset = -1*pys.mean()
    vertical_offset = pxs.mean()

    # Set up powered descent canvas
    pdCanvas = canvas(title="Powered Descent, Local Level View", width=600, height=600, background=color.cyan)
    pdCanvas.lights = []
    distant_light(canvas=pdCanvas, direction=vector(0,0,1), color=vector(0.9,0.9,0.9))
    distant_light(canvas=pdCanvas, direction=vector(0,0,-1), color=vector(0.9,0.9,0.9))
    landerx, landery = LL_2_VPython(pxs[0], pys[0])
    lander = sphere(canvas=pdCanvas, pos=vector(landerx-horizontal_offset, landery-vertical_offset, pzs[0]), radius=10, shininess=0, color=color.blue, make_trail=True, trail_radius=2)
    target = sphere(canvas=pdCanvas, pos=vector(-horizontal_offset,-vertical_offset,0), radius=10, color=color.red)
    ground = box(canvas=pdCanvas, pos=vector(-horizontal_offset, -vertical_offset, 0), length=1000, height=5, width=1000, color=color.green)
    thrustx, thrusty = LL_2_VPython(Txs[0], Tys[0])
    thrust_sf = 200
    ctrlThrusters = arrow(canvas=pdCanvas, pos=vector(landerx-horizontal_offset, landery-vertical_offset, pzs[0]), axis=vector(thrustx/thrust_sf, thrusty/thrust_sf, Tzs[0]/thrust_sf), color = color.orange)

    # Set up position graph
    posGraph = graph(title = 'Lander Position', xtitle='Global Time (s)', ytitle='Position (m)')
    xpos = gcurve(graph=posGraph, color=color.blue, label='X')
    ypos = gcurve(graph=posGraph, color=color.green, label='Y')
    zpos = gcurve(graph=posGraph, color=color.magenta, label='Z')
    xpos.plot(ts_p[0], pxs[0])
    ypos.plot(ts_p[0], pys[0])
    zpos.plot(ts_p[0], pzs[0])

    # Set up velocity graph
    velGraph = graph(title = 'Lander Velocity', xtitle='Global Time (s)', ytitle='Velocity (m)', ymax=50, ymin=-50)
    xvel = gcurve(graph=velGraph, color=color.blue, label='X')
    yvel = gcurve(graph=velGraph, color=color.green, label='Y')
    zvel = gcurve(graph=velGraph, color=color.magenta, label='Z')
    xvel.plot(ts_p[0], vxs[0])
    yvel.plot(ts_p[0], vys[0])
    zvel.plot(ts_p[0], vzs[0])

    # Set up thrust graph
    thrustGraph = graph(title = 'Lander Thrust Profile', xtitle='Global Time (s)', ytitle='Thrust (N)')
    xthrust = gcurve(graph=thrustGraph, color=color.blue, label='X')
    ythrust = gcurve(graph=thrustGraph, color=color.green, label='Y')
    zthrust = gcurve(graph=thrustGraph, color=color.magenta, label='Z')
    xvel.plot(ts_p[0], Txs[0])
    yvel.plot(ts_p[0], Tys[0])
    zvel.plot(ts_p[0], Tzs[0])

    # Set up mass graph
    massGraph = graph(title = 'Lander Mass History', xtitle='Global Time (s)', ytitle='Mass (kg)')
    mlander = gcurve(graph=massGraph, color=color.blue, label='Mass')
    mlander.plot(ts_p[0], ms[0])

    hasReachedFinalLandingTime = False
    count = 0
    while not hasReachedFinalLandingTime:
        rate(animationRate)
        count = count+1
        t = ts_p[count]

        landerx, landery = LL_2_VPython(pxs[count], pys[count])
        lander.pos = vector(landerx-horizontal_offset, landery-vertical_offset, pzs[count])

        thrustx, thrusty = LL_2_VPython(Txs[count], Tys[count])
        ctrlThrusters.pos = vector(landerx-horizontal_offset, landery-vertical_offset, pzs[count])
        ctrlThrusters.axis = vector(thrustx/thrust_sf, thrusty/thrust_sf, Tzs[count]/thrust_sf)
        ctrlThrusters.headwidth = ctrlThrusters.shaftwidth

        xpos.plot(ts_p[count], pxs[count])
        ypos.plot(ts_p[count], pys[count])
        zpos.plot(ts_p[count], pzs[count])

        xvel.plot(ts_p[count], vxs[count])
        yvel.plot(ts_p[count], vys[count])
        zvel.plot(ts_p[count], vzs[count])

        xthrust.plot(ts_p[count], Txs[count])
        ythrust.plot(ts_p[count], Tys[count])
        zthrust.plot(ts_p[count], Tzs[count])

        mlander.plot(ts_p[count], ms[count])

        if t == ts_p[-1]:
            hasReachedFinalLandingTime = True

    return pxs, pys, pzs

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

def LL_2_VPython(x, y):
    
    ynew = x
    xnew = -y

    return xnew, ynew

if __name__ == "__main__":

    # Target Drop Location
    latTarget = 38.836587*(np.pi/180)
    lonTarget = -77.196691*(np.pi/180)

    # Establish Filepaths
    pathToOrbitData = "../C++/results/orbs.csv"
    pathToBallisticData = "../C++/results/bd.csv"
    pathToChuteData = "../C++/results/cd.csv"
    pathToPoweredData = "../C++/results/pd.csv"
    
    # Set animation rates
    orbitAnimationRate = 500
    bdAnimationRate = 500
    cdAnimationRate = 150
    pdAnimationRate = 40

    print("Beginning animations for planetary supply drop...")

    visualizeOrbits(pathToOrbitData, orbitAnimationRate, latTarget, lonTarget)
    print("Orbit visualizations complete.")
    visualizeBallisticDescent(pathToBallisticData, bdAnimationRate)
    print("Ballistic descent visualization complete.")
    visualizeChuteDescent(pathToChuteData, cdAnimationRate)
    print("Chute descent visualization complete.")
    pxs, pys, pzs = visualizePoweredDescent(pathToPoweredData, pdAnimationRate)
    print("Powered descent visualization complete.")

    fig1 = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(pzs,pys,pxs)
    ax.set_title('Powered Descent 3D Trajectory')
    plt.show()

    print("Animations complete!")



# Functions for each visualization that take in: filepath to .csv, rate of animation
# Main script