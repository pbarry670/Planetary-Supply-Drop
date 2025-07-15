import numpy as np
import cvxpy as cp

def solveGfoldOptim(passedVals):
    

    state0 = [passedVals(0), passedVals(1), passedVals(2), passedVals(3), passedVals(4), passedVals(5)]

    n_T = [passedVals(6), passedVals(7), passedVals(8)]

    m_wet = passedVals(9)
    m_fuel = passedVals(10)

    T_min = passedVals(11)
    T_max = passedVals(12)
    alpha = passedVals(13)
    gamma_gs = passedVals(14)
    theta = passedVals(15)
    Vmax = passedVals(16)

    g_mag = passedVals(17)

    tf = passedVals(18)
    dt = passedVals(19)

    status, miss_dist, traj_hist, m_hist, T_hist = solve_p3(state0, n_T, m_wet, m_fuel, T_min, T_max, alpha, gamma_gs, theta, Vmax, g_mag, tf, dt)
    
    d_p3star = traj_hist[-1, 1:2]

    status, miss_dist, traj_hist, m_hist, T_hist = solve_p4(state0, n_T, m_wet, m_fuel, T_min, T_max, alpha, gamma_gs, theta, Vmax, g_mag, tf, dt, d_p3star)



    return [1.1, 2.2, 3.3]



def solve_p3(state0, n_T, m_wet, m_fuel, T_min, T_max, alpha, gamma_gs, theta, Vmax, g_mag, tf, dt):

    # Problem Data
    x0 = state0[0]
    y0 = state0[1]
    z0 = state0[2]
    xdot0 = state0[3]
    ydot0 = state0[4]
    zdot0 = state0[5]

    n_Tx = n_T[0]
    n_Ty = n_T[1]
    n_Tz = n_T[2]

    g = np.array([-1*g_mag, 0, 0])

    N = int(tf/dt)

    # Construct the problem
    x = cp.Variable((N+1, 6))
    m = cp.Variable((N+1))
    T = cp.Variable((N+1, 3))
    Gamma = cp.Variable(N+1)

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

    objective = cp.norm(x[N, 1:2]) # Minimize the yz-distance from the origin, which is the target

    # Construct the constraints
    constraints = []
    for i in range(N):
        dyn_constr = x[i+1, :] == A@x[i, :] + B@(cp.sum(g, T[i, :]))
        m_constr = m[i+1] == -alpha*Gamma[i]
        vmax_constr = cp.norm(x[i, 3:5]) <= Vmax
        glideslope_constr = cp.norm( (x[i, 1] - x[N, 1]) + (x[i,2] - x[N,2]) ) - (1/np.tan(gamma_gs))*(x[i,0] - x[N,0]) <= 0
        Tnorm_constr = cp.norm(T[i, :]) <= Gamma[i]
        Tmin_constr = Gamma[i] >= T_min
        Tmax_constr = Gamma[i] <= T_max
        Tpoint_constr = n_Tx*T[i,0] + n_Ty*T[i,1] + n_Tz*T[i,2] <= np.cos(theta)*Gamma[i]

        constraints.extend(dyn_constr, m_constr, vmax_constr, glideslope_constr, Tnorm_constr, Tmin_constr, Tmax_constr, Tpoint_constr)

    initial_m_constr = m[0] == m_wet
    final_m_constr = m[N] >= m_wet - m_fuel
    initial_xpos_constr = x[0,0] == x0
    initial_ypos_constr = x[0,1] == y0
    initial_zpos_constr = x[0,2] == z0
    initial_xvel_constr = x[0,3] == xdot0
    initial_yvel_constr = x[0,4] == ydot0
    initial_zvel_constr = x[0,5] == zdot0
    final_xpos_constr = x[N, 0] == 0
    final_xvel_constr = x[N, 3] == 0
    final_yvel_constr = x[N, 4] == 0
    final_zvel_constr = x[N, 5] == 0

    constraints.extend(initial_m_constr, final_m_constr)
    constraints.extend(initial_xpos_constr, initial_ypos_constr, initial_zpos_constr, initial_xvel_constr, initial_yvel_constr, initial_zvel_constr)
    constraints.extend(final_xpos_constr, final_xvel_constr, final_yvel_constr, final_zvel_constr)

    # Solve the problem
    p3 = cp.Problem(cp.Minimize(objective), constraints)
    p3.solve()

    status = p3.status
    miss_dist = p3.value
    traj_hist = np.array(x.value)
    m_hist = np.array(m.value)
    T_hist = np.array(T.value)

    return status, miss_dist, traj_hist, m_hist, T_hist



def solve_p4(state0, n_T, m_wet, m_fuel, T_min, T_max, alpha, gamma_gs, theta, Vmax, g_mag, tf, dt, d_p3star):

    # Problem Data
    x0 = state0[0]
    y0 = state0[1]
    z0 = state0[2]
    xdot0 = state0[3]
    ydot0 = state0[4]
    zdot0 = state0[5]

    n_Tx = n_T[0]
    n_Ty = n_T[1]
    n_Tz = n_T[2]

    g = np.array([-1*g_mag, 0, 0])

    N = int(tf/dt)

    # Construct the problem
    x = cp.Variable((N+1, 6))
    m = cp.Variable((N+1))
    T = cp.Variable((N+1, 3))
    Gamma = cp.Variable(N+1)

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
    
    objective = cp.sum(dt*Gamma[0:N]) # Minimize the thrust (and therefore fuel usage) across the trajectory

    # Construct the constraints
    constraints = []
    for i in range(N):
        dyn_constr = x[i+1, :] == A@x[i, :] + B@(cp.sum(g, T[i, :]))
        m_constr = m[i+1] == -alpha*Gamma[i]
        vmax_constr = cp.norm(x[i, 3:5]) <= Vmax
        glideslope_constr = cp.norm( (x[i, 1] - x[N, 1]) + (x[i,2] - x[N,2]) ) - (1/np.tan(gamma_gs))*(x[i,0] - x[N,0]) <= 0
        Tnorm_constr = cp.norm(T[i, :]) <= Gamma[i]
        Tmin_constr = Gamma[i] >= T_min
        Tmax_constr = Gamma[i] <= T_max
        Tpoint_constr = n_Tx*T[i,0] + n_Ty*T[i,1] + n_Tz*T[i,2] <= np.cos(theta)*Gamma[i]

        constraints.extend(dyn_constr, m_constr, vmax_constr, glideslope_constr, Tnorm_constr, Tmin_constr, Tmax_constr, Tpoint_constr)

    initial_m_constr = m[0] == m_wet
    final_m_constr = m[N] >= m_wet - m_fuel
    initial_xpos_constr = x[0,0] == x0
    initial_ypos_constr = x[0,1] == y0
    initial_zpos_constr = x[0,2] == z0
    initial_xvel_constr = x[0,3] == xdot0
    initial_yvel_constr = x[0,4] == ydot0
    initial_zvel_constr = x[0,5] == zdot0
    final_xpos_constr = x[N, 0] == 0
    final_xvel_constr = x[N, 3] == 0
    final_yvel_constr = x[N, 4] == 0
    final_zvel_constr = x[N, 5] == 0

    constraints.extend(initial_m_constr, final_m_constr)
    constraints.extend(initial_xpos_constr, initial_ypos_constr, initial_zpos_constr, initial_xvel_constr, initial_yvel_constr, initial_zvel_constr)
    constraints.extend(final_xpos_constr, final_xvel_constr, final_yvel_constr, final_zvel_constr)

    solution_constr = cp.norm(x[N, 1:2]) <= cp.norm(d_p3star)
    constraints.extend(solution_constr)

    # Solve the problem
    p4 = cp.Problem(cp.Minimize(objective), constraints)
    p4.solve()

    status = p4.status
    miss_dist = p4.value
    traj_hist = np.array(x.value)
    m_hist = np.array(m.value)
    T_hist = np.array(T.value)

    return status, miss_dist, traj_hist, m_hist, T_hist


