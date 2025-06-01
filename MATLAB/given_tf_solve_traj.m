function [generated_traj]= given_tf_solve_traj(tf, dt, gvec, sc_params)
    %% Solver
    tf = tf;
    dt = dt; 
    g = gvec;
    m_wet = sc_params.m_wet; % kg
    Isp = sc_params.Isp; % s
    alpha = sc_params.alpha;
    rho_1 = sc_params.Tmin; % N
    rho_2 = sc_params.Tmax; % N
    ts = 0:dt:tf;
    r0 = sc_params.r0;
    r0dot = sc_params.r0dot;
    theta_alt = sc_params.theta_alt; % 87 deg
    S = [0, 1, 0, 0, 0, 0;
         0, 0, 1, 0, 0, 0];
    c = [-tand(theta_alt); 0; 0; 0; 0; 0];
    
    N = tf/dt; % number of control inputs; number of temporal nodes minus 1
    % There are N+1 temporal nodes
    
    cvx_begin
        variable u(3,N)
        variable sigm(N)
        expression x(6,N+1)
        expression z(N+1)
        expression z0(N+1)
        expression mu1(N+1)
        expression mu2(N+1)
        for k = 1:N+1
            z0(k) = log(m_wet - alpha*rho_2*ts(k));
            mu_1(k) = rho_1*exp(-z0(k));
            mu_2(k) = rho_2*exp(-z0(k));
        end
        z(1) = log(m_wet); % initial mass
        x(1:3, 1) = r0; % initial position
        x(4:6, 1) = r0dot; % initial velocity

        minimize dt*sum(sigm);
        subject to
            x(:, end) = [0;0;0;0;0;0]; % final position/velocity
            for k = 1:N
                x(1, k+1) = x(1, k) + x(4,k)*dt + 0.5*(dt^2)*(u(1,k) + g(1)); % Dynamics
                x(2, k+1) = x(2, k) + x(5,k)*dt + 0.5*(dt^2)*(u(2,k) + g(2));
                x(3, k+1) = x(3, k) + x(6,k)*dt + 0.5*(dt^2)*(u(3,k) + g(3));
                x(4, k+1) = x(4,k) + dt*(u(1,k) + g(1));
                x(5, k+1) = x(5,k) + dt*(u(2,k) + g(2));
                x(6, k+1) = x(6,k) + dt*(u(3,k) + g(3));
                z(k+1) = z(k) + dt*(-alpha*sigm(k)); % Change in z
    
                norm(u(:, k)) <= sigm(k); % thrust constraint
                sigm(k) >= mu_1(k)*(1 - (z(k) - z0(k)) + ((z(k) - z0(k))^2/2)); % thrust lower bound
                sigm(k) <= mu_2(k)*(1 - (z(k) - z0(k))); % thrust upper bound
                z(k) >= log(m_wet - alpha*rho_2*ts(k)); % z lower bound
                z(k) <= log(m_wet - alpha*rho_1*ts(k)); % z upper bound
                x(1, k) >= 0; % subsurface constraint
                norm(S*x(:, k)) + c'*x(:,k) <= 0; % glide slope constraint
            end
    cvx_end

    generated_traj.x = x;
    generated_traj.z = z;
    generated_traj.u = u;
    generated_traj.sigm = sigm;
    generated_traj.ts = ts;
    generated_traj.objval = cvx_optval;
end