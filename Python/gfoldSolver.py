import numpy as np
import cvxpy as cp

def commTest(passedVals):
    x = otherTest(passedVals)
    return x

def otherTest(passedVals):
    x = [passedVals[0], passedVals[1], passedVals[2]]
    return x

def solveGfoldOptim(passedVals):
    

    state0 = [passedVals[0], passedVals[1], passedVals[2], passedVals[3], passedVals[4], passedVals[5]]

    n_T = [passedVals[6], passedVals[7], passedVals[8]]

    m_wet = passedVals[9]
    m_fuel = passedVals[10]

    T_min = passedVals[11]
    T_max = passedVals[12]
    alpha = passedVals[13]
    gamma_gs = passedVals[14]
    theta = passedVals[15]
    Vmax = passedVals[16]

    g_mag = passedVals[17]

    tf = passedVals[18]
    dt = passedVals[19]

    status_flag = 1
    status, miss_dist, traj_hist, z_hist, u_hist = solve_p3(state0, n_T, m_wet, m_fuel, T_min, T_max, alpha, gamma_gs, theta, Vmax, g_mag, tf, dt)
    print(f"Problem 3 Status: '{status}'")
    if "optimal" not in status:
        print("Problem 3 solution not found. Consider relaxing powered descent constraints.")
        status_flag = 0
        return [status_flag, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

    d_p3star = traj_hist[-1, 1:3] # Last index is exclusive
    status, final_z, traj_hist, z_hist, u_hist = solve_p4(state0, n_T, m_wet, m_fuel, T_min, T_max, alpha, gamma_gs, theta, Vmax, g_mag, tf, dt, d_p3star)
    print(f"Problem 4 Status: '{status}'")
    if "optimal" not in status:
        print("Problem 4 solution not found. Consider relaxing powered descent constraints.")
        status_flag = 0
        return [status_flag, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

    numel_z = z_hist.size

    print("Synthesizing trajectory optimization results...")
    rx_hist = traj_hist[:,0]
    ry_hist = traj_hist[:,1]
    rz_hist = traj_hist[:,2]
    vx_hist = traj_hist[:,3]
    vy_hist = traj_hist[:,4]
    vz_hist = traj_hist[:,5]
    ux_hist = u_hist[:,0]
    uy_hist = u_hist[:,1]
    uz_hist = u_hist[:,2]

    rx_hist = np.squeeze(rx_hist.reshape((1, numel_z)))
    ry_hist = np.squeeze(ry_hist.reshape((1, numel_z)))
    rz_hist = np.squeeze(rz_hist.reshape((1, numel_z)))
    vx_hist = np.squeeze(vx_hist.reshape((1, numel_z)))
    vy_hist = np.squeeze(vy_hist.reshape((1, numel_z)))
    vz_hist = np.squeeze(vz_hist.reshape((1, numel_z)))

    ux_hist = np.squeeze(ux_hist.reshape((1, numel_z)))
    uy_hist = np.squeeze(uy_hist.reshape((1, numel_z)))
    uz_hist = np.squeeze(uz_hist.reshape((1, numel_z)))

    z_hist = np.squeeze(z_hist.reshape(1, numel_z))

    master_array = np.concatenate((rx_hist, ry_hist, rz_hist, vx_hist, vy_hist, vz_hist, z_hist, ux_hist, uy_hist, uz_hist), axis=0)
    master_array = np.insert(master_array, 0, status_flag)
    master_list = list(master_array)

    return master_list


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
        vmax_constr = cp.norm(x[i, 3:6]) <= Vmax # Check this for indexing
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