import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt

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
    initial_xpos_constr = x[0,0] == x0
    initial_ypos_constr = x[0,1] == y0
    initial_zpos_constr = x[0,2] == z0
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
    p3.solve(verbose=True)

    status = p3.status
    miss_dist = p3.value
    traj_hist = np.array(x.value)
    z_hist = np.array(z.value)
    u_hist = np.array(u.value)

    return status, miss_dist, traj_hist, z_hist, u_hist



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
        dyn_constr = x[i+1, :] == A@x[i, :] + B@g + B@(T[i,:]/m[i])
        m_constr = m[i+1] == -alpha*Gamma[i]
        vmax_constr = cp.norm(x[i, 3:5]) <= Vmax
        glideslope_constr = cp.norm( (x[i, 1] - x[N, 1]) + (x[i,2] - x[N,2]) ) - (1/np.tan(gamma_gs))*(x[i,0] - x[N,0]) <= 0
        Tnorm_constr = cp.norm(T[i, :]) <= Gamma[i]
        Tmin_constr = Gamma[i] >= T_min
        Tmax_constr = Gamma[i] <= T_max
        Tpoint_constr = n_Tx*T[i,0] + n_Ty*T[i,1] + n_Tz*T[i,2] <= np.cos(theta)*Gamma[i]

        constraints.extend([dyn_constr, m_constr, vmax_constr, glideslope_constr, Tnorm_constr, Tmin_constr, Tmax_constr, Tpoint_constr])

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

    constraints.extend([initial_m_constr, final_m_constr])
    constraints.extend([initial_xpos_constr, initial_ypos_constr, initial_zpos_constr, initial_xvel_constr, initial_yvel_constr, initial_zvel_constr])
    constraints.extend([final_xpos_constr, final_xvel_constr, final_yvel_constr, final_zvel_constr])

    solution_constr = cp.norm(x[N, 1:2]) <= cp.norm(d_p3star)
    constraints.extend([solution_constr])

    # Solve the problem
    p4 = cp.Problem(cp.Minimize(objective), constraints)
    p4.solve()

    status = p4.status
    miss_dist = p4.value
    traj_hist = np.array(x.value)
    m_hist = np.array(m.value)
    T_hist = np.array(T.value)

    return status, miss_dist, traj_hist, m_hist, T_hist

if __name__ == "__main__":
    

    #state0 = [2000, 200, 100, -50, 2, -3]
    state0 = [1000, 400, -300, -34, 10, -15]

    n_T = [-1, 0, 0]

    m_wet = 2000 # 35600
    m_fuel = 1200 # 10000

    T_min = 0.3*20000 #0.3*5000 #0.4*411000
    T_max = 20000 # 411000
    alpha = 1/(311*9.807)
    gamma_gs = 20*(np.pi/180)
    theta = 75*(np.pi/180)
    Vmax = 100

    g_mag = 9.807

    tf = 75 # 50
    dt = 0.1

    status, miss_dist, traj_hist, z_hist, u_hist = solve_p3(state0, n_T, m_wet, m_fuel, T_min, T_max, alpha, gamma_gs, theta, Vmax, g_mag, tf, dt)
    print(status)
    print(miss_dist)
    print(traj_hist[-1, 1])
    print(traj_hist[-1, 2])


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

    #status, miss_dist, traj_hist, m_hist, T_hist = solve_p4(state0, n_T, m_wet, m_fuel, T_min, T_max, alpha, gamma_gs, theta, Vmax, g_mag, tf, dt, d_p3star)


