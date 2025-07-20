import numpy as np
import pandas as pd
from vpython import canvas, sphere, arrow, color, vector, distant_light, graph, gcurve, rate
import matplotlib.pyplot as plt

def LL_2_VPython(x, y):
    
    ynew = x
    xnew = -y

    return xnew, ynew


df_p = pd.read_csv("../C++/results/pd.csv", header=None)
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

fig1 = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(pzs,pys,pxs)
ax.set_title('3D Trajectory')


# Set up powered descent canvas
pdCanvas = canvas(title="Powered Descent, Local Level View", width=300, height=300, background=color.white)
pdCanvas.lights = []
distant_light(canvas=pdCanvas, direction=vector(0,0,1), color=vector(0.9,0.9,0.9))
distant_light(canvas=pdCanvas, direction=vector(0,0,-1), color=vector(0.9,0.9,0.9))
landerx, landery = LL_2_VPython(pxs[0], pys[0])
lander = sphere(canvas=pdCanvas, pos=vector(landerx, landery, pzs[0]), radius=10, shininess=0, color=color.blue, make_trail=True, trail_radius=2)
target = sphere(canvas=pdCanvas, pos=vector(0,0,0), radius=10, color=color.red)
thrustx, thrusty = LL_2_VPython(Txs[0], Tys[0])
thrust_sf = 200
ctrlThrusters = arrow(canvas=pdCanvas, pos=vector(landerx, landery, pzs[0]), axis=vector(thrustx/thrust_sf, thrusty/thrust_sf, Tzs[0]/thrust_sf), color = color.orange)

# Set up position graph
posGraph = graph(title = 'Lander Position', xtitle='Global Time (s)', ytitle='Position (m)')
xpos = gcurve(graph=posGraph, color=color.blue)
ypos = gcurve(graph=posGraph, color=color.green)
zpos = gcurve(graph=posGraph, color=color.magenta)
xpos.plot(ts_p[0], pxs[0])
ypos.plot(ts_p[0], pys[0])
zpos.plot(ts_p[0], pzs[0])

print(vxs[0])
print(vys[0])
print(vzs[0])
# Set up velocity graph
velGraph = graph(title = 'Lander Velocity', xtitle='Global Time (s)', ytitle='Velocity (m)', ymax=50, ymin=-50)
xvel = gcurve(graph=velGraph, color=color.blue)
yvel = gcurve(graph=velGraph, color=color.green)
zvel = gcurve(graph=velGraph, color=color.magenta)
xvel.plot(ts_p[0], vxs[0])
yvel.plot(ts_p[0], vys[0])
zvel.plot(ts_p[0], vzs[0])

# Set up thrust graph
thrustGraph = graph(title = 'Lander Thrust Profile', xtitle='Global Time (s)', ytitle='Thrust (N)')
xthrust = gcurve(graph=thrustGraph, color=color.blue)
ythrust = gcurve(graph=thrustGraph, color=color.green)
zthrust = gcurve(graph=thrustGraph, color=color.magenta)
xvel.plot(ts_p[0], Txs[0])
yvel.plot(ts_p[0], Tys[0])
zvel.plot(ts_p[0], Tzs[0])

# Set up mass graph
massGraph = graph(title = 'Lander Mass History', xtitle='Global Time (s)', ytitle='Mass (kg)')
mlander = gcurve(graph=massGraph, color=color.blue)
mlander.plot(ts_p[0], ms[0])

hasReachedFinalLandingTime = False
count = 0
while not hasReachedFinalLandingTime:
    rate(60)
    count = count+1
    t = ts_p[count]

    landerx, landery = LL_2_VPython(pxs[count], pys[count])
    lander.pos = vector(landerx, landery, pzs[count])

    thrustx, thrusty = LL_2_VPython(Txs[count], Tys[count])
    ctrlThrusters.pos = vector(landerx, landery, pzs[count])
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

print("Powered descent visualization completed! \n")
plt.show()
