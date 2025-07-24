import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt

def solve_p3(state0, n_T, m_wet, m_fuel, T_min, T_max, alpha, gamma_gs, theta, Vmax, g_mag, tf, dt):

    # Problem Data
    xpos0 = state0[0]
    ypos0 = state0[1]
    zpos0 = state0[2]
    xdot0 = state0[3]
    ydot0 = state0[4]
    zdot0 = state0[5]

    n_Tx = n_T[0]
    n_Ty = n_T[1]
    n_Tz = n_T[2]

    #g = np.array([-1*g_mag, 0, 0])
    g = np.zeros(3)
    g[0] = -1*g_mag

    N = int(tf/dt)

    # Construct the problem
    x = cp.Variable((N+1, 6))
    z = cp.Variable((N+1))
    u = cp.Variable((N+1, 3))
    sigma = cp.Variable(N+1)

    A = np.array([[0, 0, 0, 1, 0, 0],
                  [0, 0, 0, 0, 1, 0],
                  [0, 0, 0, 0, 0, 1],
                  [0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0]])


    B = np.array([[0, 0, 0],
                  [0, 0, 0],
                  [0, 0, 0],
                  [1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 1]])
    

    objective = cp.norm(x[N, 1:3]) # Minimize the yz-distance from the origin, which is the target (indexing is exclusive for the end)
    
    # Construct the constraints
    constraints = []
    for i in range(N):
        dyn_constr = x[i+1, :] == x[i] + dt*(A@x[i, :] + B@g + B@(u[i, :]))
        m_constr = z[i+1] == z[i] -alpha*sigma[i]*dt
        vmax_constr = cp.norm(x[i, 3:6]) <= Vmax
        glideslope_constr = x[i, 0]*np.cos(gamma_gs) >= cp.norm(x[i, 1:3])
        Tnorm_constr = cp.norm(u[i, :]) <= sigma[i]

        z0 = np.log(m_wet - alpha*T_max*i*dt)

        Tmin_constr = sigma[i] >= T_min*np.exp(-z0)*(1 - (z[i] - z0) + 0.5*(z[i] - z0)**2)
        Tmax_constr = sigma[i] <= T_max*np.exp(-z0)*(1 - (z[i] - z0))
        Tpoint_constr = n_Tx*u[i,0] + n_Ty*u[i,1] + n_Tz*u[i,2] <= np.cos(theta)*sigma[i]

        constraints.extend([dyn_constr, m_constr, Tmin_constr, Tmax_constr, Tnorm_constr, vmax_constr, glideslope_constr, Tpoint_constr])

    initial_m_constr = z[0] == np.log(m_wet)
    final_m_constr = z[N] >= np.log(m_wet - m_fuel)
    initial_xpos_constr = x[0,0] == xpos0
    initial_ypos_constr = x[0,1] == ypos0
    initial_zpos_constr = x[0,2] == zpos0
    initial_xvel_constr = x[0,3] == xdot0
    initial_yvel_constr = x[0,4] == ydot0
    initial_zvel_constr = x[0,5] == zdot0
    final_xpos_constr = x[-1, 0] == 0
    final_xvel_constr = x[-1, 3] == 0
    final_yvel_constr = x[-1, 4] == 0
    final_zvel_constr = x[-1, 5] == 0

    constraints.extend([initial_m_constr, final_m_constr])
    constraints.extend([initial_xpos_constr, initial_ypos_constr, initial_zpos_constr, initial_xvel_constr, initial_yvel_constr, initial_zvel_constr])
    constraints.extend([final_xpos_constr, final_xvel_constr, final_yvel_constr, final_zvel_constr])

    # Solve the problem
    p3 = cp.Problem(cp.Minimize(objective), constraints)
    print("Solving minimum landing error problem (P3)...")
    p3.solve(verbose=False)

    status = p3.status
    miss_dist = p3.value
    traj_hist = np.array(x.value)
    z_hist = np.array(z.value)
    u_hist = np.array(u.value)

    return status, miss_dist, traj_hist, z_hist, u_hist



def solve_p4(state0, n_T, m_wet, m_fuel, T_min, T_max, alpha, gamma_gs, theta, Vmax, g_mag, tf, dt, d_p3star):

     # Problem Data
    xpos0 = state0[0]
    ypos0 = state0[1]
    zpos0 = state0[2]
    xdot0 = state0[3]
    ydot0 = state0[4]
    zdot0 = state0[5]

    n_Tx = n_T[0]
    n_Ty = n_T[1]
    n_Tz = n_T[2]

    #g = np.array([-1*g_mag, 0, 0])
    g = np.zeros(3)
    g[0] = -1*g_mag

    N = int(tf/dt)

    # Construct the problem
    x = cp.Variable((N+1, 6))
    z = cp.Variable((N+1))
    u = cp.Variable((N+1, 3))
    sigma = cp.Variable(N+1)

    A = np.array([[0, 0, 0, 1, 0, 0],
                  [0, 0, 0, 0, 1, 0],
                  [0, 0, 0, 0, 0, 1],
                  [0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0]])


    B = np.array([[0, 0, 0],
                  [0, 0, 0],
                  [0, 0, 0],
                  [1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 1]])
    

    objective = cp.sum(z[N])
    
    # Construct the constraints
    constraints = []
    for i in range(N):
        dyn_constr = x[i+1, :] == x[i] + dt*(A@x[i, :] + B@g + B@(u[i, :]))
        m_constr = z[i+1] == z[i] -alpha*sigma[i]*dt
        vmax_constr = cp.norm(x[i, 3:6]) <= Vmax
        glideslope_constr = x[i, 0]*np.cos(gamma_gs) >= cp.norm(x[i, 1:3])
        Tnorm_constr = cp.norm(u[i, :]) <= sigma[i]

        z0 = np.log(m_wet - alpha*T_max*i*dt)

        Tmin_constr = sigma[i] >= T_min*np.exp(-z0)*(1 - (z[i] - z0) + 0.5*(z[i] - z0)**2)
        Tmax_constr = sigma[i] <= T_max*np.exp(-z0)*(1 - (z[i] - z0))
        Tpoint_constr = n_Tx*u[i,0] + n_Ty*u[i,1] + n_Tz*u[i,2] <= np.cos(theta)*sigma[i]

        constraints.extend([dyn_constr, m_constr, Tmin_constr, Tmax_constr, Tnorm_constr, vmax_constr, glideslope_constr, Tpoint_constr])

    initial_m_constr = z[0] == np.log(m_wet)
    final_m_constr = z[N] >= np.log(m_wet - m_fuel)
    initial_xpos_constr = x[0,0] == xpos0
    initial_ypos_constr = x[0,1] == ypos0
    initial_zpos_constr = x[0,2] == zpos0
    initial_xvel_constr = x[0,3] == xdot0
    initial_yvel_constr = x[0,4] == ydot0
    initial_zvel_constr = x[0,5] == zdot0
    final_xpos_constr = x[-1, 0] == 0
    final_xvel_constr = x[-1, 3] == 0
    final_yvel_constr = x[-1, 4] == 0
    final_zvel_constr = x[-1, 5] == 0

    constraints.extend([initial_m_constr, final_m_constr])
    constraints.extend([initial_xpos_constr, initial_ypos_constr, initial_zpos_constr, initial_xvel_constr, initial_yvel_constr, initial_zvel_constr])
    constraints.extend([final_xpos_constr, final_xvel_constr, final_yvel_constr, final_zvel_constr])

    solution_constr = cp.norm(x[N, 1:3]) <= cp.norm(d_p3star)
    constraints.extend([solution_constr])

    # Solve the problem
    p4 = cp.Problem(cp.Maximize(objective), constraints)
    print("Solving the minimum fuel problem (P4)...")
    p4.solve(verbose=False)

    status = p4.status
    final_z = p4.value
    traj_hist = np.array(x.value)
    z_hist = np.array(z.value)
    u_hist = np.array(u.value)

    return status, final_z, traj_hist, z_hist, u_hist

if __name__ == "__main__":
    

    #state0 = [2000, 200, 100, -50, 2, -3]
    state0 = [999.582, 82.422, -89.20907, -33.9954, 1.01036, -3.031]

    n_T = [-1.0, 0.0, 0.0]

    m_wet = 2000.0 # 35600
    m_fuel = 1200.0 # 10000

    T_min = 0.3*20000 #0.3*5000 #0.4*411000
    T_max = 20000.0 # 411000
    alpha = 1/(300*9.807)
    gamma_gs = 20*(np.pi/180)
    theta = 75*(np.pi/180)
    Vmax = 150.0

    g_mag = 9.807

    tf = 75.0
    dt = 0.5

    status_flag = 1
    status, miss_dist, traj_hist, z_hist, u_hist = solve_p3(state0, n_T, m_wet, m_fuel, T_min, T_max, alpha, gamma_gs, theta, Vmax, g_mag, tf, dt)
    print(f"Problem 3 Status: '{status}'")
    if "optimal" not in status:
        print("Problem 3 solution not found. Consider relaxing powered descent constraints.")
        status_flag = 0

    d_p3star = traj_hist[-1, 1:3] # Last index is exclusive

    status, final_z, traj_hist, z_hist, u_hist = solve_p4(state0, n_T, m_wet, m_fuel, T_min, T_max, alpha, gamma_gs, theta, Vmax, g_mag, tf, dt, d_p3star)
    print(f"Problem 4 Status: '{status}'")
    if "optimal" not in status:
        print("Problem 4 solution not found. Consider relaxing powered descent constraints.")
        status_flag = 0

    time = np.linspace(0, int(tf), int(tf/dt)+1)

    fig1 = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(traj_hist[:,2],traj_hist[:,1],traj_hist[:,0])
    ax.set_title('3D Trajectory')

    fig2, (ax1,ax2,ax3) = plt.subplots(3,1)
    ax1.plot(time,traj_hist[:,0])
    ax2.plot(time,traj_hist[:,1])
    ax3.plot(time,traj_hist[:,2])
    fig2.suptitle('Position')

    fig3, (ax1,ax2,ax3) = plt.subplots(3,1)
    ax1.plot(time,traj_hist[:,3])
    ax2.plot(time,traj_hist[:,4])
    ax3.plot(time,traj_hist[:,5])
    fig3.suptitle('Velocity')

    fig4, (ax1,ax2,ax3) = plt.subplots(3,1)
    ax1.plot(time,u_hist[:,0])
    ax2.plot(time,u_hist[:,1])
    ax3.plot(time,u_hist[:,2])
    fig4.suptitle('Thrust/Mass')

    plt.show()


