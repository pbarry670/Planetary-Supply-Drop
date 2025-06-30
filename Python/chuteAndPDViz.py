import numpy as np
import pandas as pd
from vpython import canvas, sphere, color, vector, distant_light, rate
import time
#TODO: Once powered descent data is available, concatenate it to the chute descent data and merge it into the visualization.

def LL_2_VPython(x, y):
    
    ynew = x
    xnew = -y

    return xnew, ynew

df_c = pd.read_csv("../C++/cd.csv", header=None)
ts_c = np.array(df_c[0])
vxs = np.array(df_c[1])
pxs = np.array(df_c[2])
vys = np.array(df_c[3])
pys = np.array(df_c[4])
vzs = np.array(df_c[5])
pzs = np.array(df_c[6])

# Set up chute descent local-level canvas
cdCanvas = canvas(title="Chute Descent, Local Level View", width=300, height=300, background=color.cyan)
cdTopCanvas = canvas(title="Chute Descent, Birds' Eye View", width=300, height=300, background=color.cyan)

# Set up chute descent local-level canvas
cdCanvas.lights = []
distant_light(canvas=cdCanvas, direction=vector(0,0,1), color=vector(0.9,0.9,0.9))
capx, capy = LL_2_VPython(pxs[0], pys[0])
capsule_c = sphere(canvas=cdCanvas, pos=vector(capx, capy, pzs[0]), radius=100, shininess=0, color=color.blue, make_trail=True, trail_radius=50)
target_c = sphere(canvas=cdCanvas, pos=vector(0,0,0), radius=50, color=color.red)

# Set up chute descent birds-eye canvas
cdTopCanvas.lights = []
distant_light(canvas=cdTopCanvas, direction=vector(0,0,1), color=vector(0.9,0.9,0.9))
capsule_ctop = sphere(canvas=cdTopCanvas, pos=vector(pys[0], pzs[0], pxs[0]), radius=100, shininess=0, color=color.blue, make_trail=True, trail_radius=10)
target_ctop = sphere(canvas=cdTopCanvas, pos=vector(0,0,0), radius=50, color=color.red)

hasReachedFinalChuteTime = False
count = 0
while not hasReachedFinalChuteTime:
    rate(250)
    count = count+1
    t = ts_c[count]

    capx, capy = LL_2_VPython(pxs[count], pys[count])
    capsule_c.pos = vector(capx, capy, pzs[count])

    capsule_ctop.pos = vector(pys[count], pzs[count], pxs[count])

    if t == ts_c[-1]:
        hasReachedFinalChuteTime = True

print("Chute descent visualization completed! \n")