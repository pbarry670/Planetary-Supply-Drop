clear; clc; close all;
set(0, 'DefaultFigureWindowStyle','docked')
set(0, 'DefaultLineLineWidth', 2);

dt = 0.05;

% Capsule Parameters
capsule.Rp = 6378*1000; % Earth radius, mm
capsule.h = 500*1000; % orbital radius, m
capsule.m = 2530; % capsule mass, kg
capsule.m_heatshield = 500; % heatshield mass, kg
capsule.m_parachute = 30; % parachute mass, kg
capsule.m_fuel = 900; % kg
capsule.m_wet = capsule.m-capsule.m_heatshield-capsule.m_parachute;
capsule.CD = 1.1; % capsule drag coeff
capsule.A = 2.5; % m^2
capsule.beta = capsule.m/(capsule.CD*capsule.A); % capsule ballistic coefficient
tspan_ballistic = 0:dt:5000;

mu_Earth = 3.986e14; %m^3/s^2
V_0_ballistic = sqrt(mu_Earth/(capsule.Rp+capsule.h)); % initial velocity, m/s
gamma_0_ballistic = 5 * (pi/180); % initial flight path angle, radians
h_0_ballistic = capsule.h; % starting altitude, m

x0_ballistic = [V_0_ballistic; gamma_0_ballistic; h_0_ballistic];

options_ballistic = odeset('Events', @ballistic_event);
[ts_ballistic, xs_ballistic] = ode45(@(t,x)ballistic_dynamics(x,capsule), tspan_ballistic, x0_ballistic, options_ballistic);

capsule.m = capsule.m - capsule.m_heatshield;
V_0_chute = xs_ballistic(end, 1);
h_0_chute = xs_ballistic(end, 3);
vy0 = 0; % m/s
py0 = 0; % m
vz0 = 0; % m/s
pz0 = 0; % m
x0_chute = [V_0_chute, h_0_chute, vy0, py0, vz0, pz0];
tspan_chute = ts_ballistic(end) + dt:dt:2000;
options_chute = odeset('Events', @chute_event);
[ts_chute, xs_chute] = ode45(@(t,x)chute_dynamics(x, capsule), tspan_chute, x0_chute, options_chute);

edl_ts = [ts_ballistic; ts_chute];
edl_Vs = [xs_ballistic(:, 1); xs_chute(:, 1)];
edl_hs = [xs_ballistic(:,3); xs_chute(:,2)];
edl_gammas = [xs_ballistic(:,2); (pi/2)*ones(length(ts_chute), 1)];
chute_deploy_t = ts_ballistic(end);

dthetas = [];
for i = 1:length(ts_ballistic)
    gamma_i = xs_ballistic(i,2); % flight path angle
    V_i = xs_ballistic(i,1); % velocity
    r_i = 6378*1000 + xs_ballistic(i,3); % distance from center of Earth (m)

    dtheta_i = (V_i*cos(gamma_i)*dt)/r_i;
    % This comes from V_theta (component of velo parallel to horizon) being
    % such that V_theta = V*cos(gamma). Then the arc length covered in that
    % time step is V_theta*dt, which is equal to r*dtheta, where r is the
    % distance from the center of the Earth to that point.

    dthetas = [dthetas, dtheta_i];

end
DTheta = sum(dthetas);
disp(['Total Delta-Theta For Ballistic Phase: ', num2str(DTheta*(180/pi)), ' degrees.' ])
disp(['Total Time for Ballistic/Chute Phase: ', num2str(ts_chute(end)), ' seconds.'])

figure(1); clf
plot(edl_ts, edl_Vs)
xlabel('Time (s)', 'FontSize', 20)
ylabel('Descent Speed (m/s)', 'FontSize', 20)
title('Descent Speed Vs Time', 'FontSize', 20)
hold on
xline(chute_deploy_t)
grid on

figure(2); clf
plot(edl_ts, edl_gammas.*(180/pi))
xlabel('Time (s)', 'FontSize', 20)
ylabel('Flight Path Angle (deg)', 'FontSize', 20)
title('Flight Path Angle Vs Time', 'FontSize', 20)
hold on
xline(chute_deploy_t)
grid on

figure(3); clf
plot(edl_ts, edl_hs)
xlabel('Time (s)', 'FontSize', 20)
ylabel('Altitude (m)', 'FontSize', 20)
title('Altitude vs Time', 'FontSize', 20)
hold on
xline(chute_deploy_t)
grid on


%% Now, Powered Descent Phase
dt = 0.5;
gvec = [-9.807; 0; 0]; % For Earth

capsule.Isp = 315; % s
capsule.Tmin = 1000; % N
capsule.Tmax =  5000; % N
capsule.r0 = [xs_chute(end, 2); xs_chute(end, 4); xs_chute(end, 6)]; % m
capsule.r0dot = [-xs_chute(end,1); xs_chute(end,3); xs_chute(end,5)]; % m/s
capsule.theta_alt = 87;
capsule.alpha = 1/(capsule.Isp*9.807);

tfmin = ((capsule.m_wet - capsule.m_fuel)*norm(capsule.r0dot))/capsule.Tmax;
tfmax = capsule.m_fuel/(capsule.alpha*capsule.Tmin);
tf = 100;

traj_star = given_tf_solve_traj(tf, dt, gvec, capsule);

figure(4); clf
plot3(traj_star.x(2,:), traj_star.x(3,:), traj_star.x(1,:))
xlabel('Y Position (m)', 'FontSize', 20)
ylabel('Z Position (m)', 'FontSize', 20)
zlabel('X Position (m)', 'FontSize', 20)
grid on
box on
title('Spacecraft Position', 'FontSize', 20)
hold on
plot3(0, 0, 0, 'o','Color','b','MarkerSize',10,...
    'MarkerFaceColor','#D9FFFF')
plot3(traj_star.x(2,1), traj_star.x(3,1), traj_star.x(1,1), 'o', 'Color', 'g', 'MarkerSize', 10, 'MarkerFaceColor', 'g')
plot3(traj_star.x(2,end), traj_star.x(3,end), traj_star.x(1,end), 'o', 'Color', 'r', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
legend('Trajectory', 'Goal', 'Start', 'Finish', 'FontSize', 15)


%% Functions
function [xdot]= ballistic_dynamics(x, capsule)

    V = x(1);
    gamma = x(2);
    h = x(3);
    
    rho = 1.225*exp(-h/8420);
    % 1.225 kg/m^3 is sea level air density
    % 8420 m is the "effective height of atmosphere"
    % This is a crude atmospheric model
    
    g_actual = 9.807*(capsule.Rp/(capsule.Rp + h))^2;
    % gravitational model
    
    Vdot = -(rho*(V^2))/(2*capsule.beta) + g_actual*sin(gamma);
    gammadot = ((-((V^2)*cos(gamma))/(capsule.Rp + h)) + g_actual*cos(gamma))/V;
    hdot = -V*sin(gamma);
    
    xdot = [Vdot; gammadot; hdot];

end

function [value, isterminal, direction] = ballistic_event(t,x)
    
    value = x(3) - 5000; % event is when h = 5000 m
    % 5000 m is parachute deployment altitude
    isterminal = 1; % halt integration
    direction =  -1; % value of h - desired_parachute_deploy_alt is decreasing

end

function [value, isterminal, direction] = chute_event(t,x)

    value = x(2) - 1000; % event is when h = 1000 m
    % 1000 m is end of chute descent and beginning of powered descent
    isterminal = 1;
    direction = -1; % value of h - desired_pwr_descent_start_alt is decreasing

end

function [xdot]= chute_dynamics(x, capsule)

    V = x(1);
    h = x(2);
    vy = x(3);
    vz = x(5);

    rho = 1.225*exp(-h/8420);
    % 1.225 kg/m^3 is sea level air density
    % 8420 m is the "effective height of atmosphere"
    % This is a crude atmospheric model

    d = 7.5; % parachute diameter, m
    A = pi*(d/2)^2; % parachute area, m^2 
    CD = 0.8;  % capsule + parachute drag coefficient
    F_D = 0.5 * rho * (V^2) * CD * A; % parachute drag force
    g = 9.807; % m/s^2

    Vdot = -F_D/capsule.m + g;
    hdot = -V;

    vydot = randn+0.02; % velocity perturbation
    pydot = vy + 0.05*vydot; % position rate of change
    vzdot = randn-0.04; % velocity perturbation
    pzdot = vz + 0.05*vzdot; % position rate of change

    xdot = [Vdot; hdot; vydot; pydot; vzdot; pzdot];

end


