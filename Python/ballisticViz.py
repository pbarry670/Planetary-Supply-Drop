import numpy as np
import pandas as pd
from vpython import canvas, sphere, color, vector, distant_light, rate

def LL_2_VPython(x, y):
    
    ynew = x
    xnew = -y

    return xnew, ynew

# Read ballistic descent data from a .csv file
df_bd = pd.read_csv("../C++/bd.csv", header=None, names = ['t', 'V', 'gamma', 'h', 'totalDTheta'])
numrows_bd, numcols_bd = df_bd.shape

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

# Set up ballistic descent canvas
scene = canvas(title="Ballistic Descent, Local Level View", width=1000, height=400, background=color.cyan)
scene.lights = []
distant_light(canvas=scene, direction=vector(0,0,1), color=vector(0.9,0.9,0.9))
capx, capy = LL_2_VPython(positions[0, 0], positions[0,1])
capsule_b = sphere(canvas=scene, pos=vector(capx, capy, 0), radius=20000, shininess=0, color=color.blue, make_trail=True, trail_color=color.blue)
target_b = sphere(canvas=scene, pos=vector(0, 0, 0), radius=50000, shininess=0, color=color.red)

hasReachedFinalBallisticTime = False
count = 0
while not hasReachedFinalBallisticTime:
    rate(1500)
    count = count + 1
    t = ts_b[count]
    capx, capy = LL_2_VPython(positions[count, 0], positions[count,1])
    capsule_b.pos = vector(capx, capy, 0)

    if t == ts_b[-1]:
        hasReachedFinalBallisticTime = True

print("Ballistic descent visualization completed! \n")












