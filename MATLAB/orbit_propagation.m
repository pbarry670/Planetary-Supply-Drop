clear; clc; close all

set(0, 'DefaultFigureWindowStyle','docked')
set(0, 'DefaultLineLineWidth', 2);

%% Satellite Orbit Parameters
altitude_orbit = 500; % km
r_orbit = altitude_orbit*1000 + 6378*1000; % m
mu_Earth = 3.986e14; % m^3/s^2
speed_orbit = sqrt(mu_Earth/r_orbit); % m/s
m_Sat = 3000; % kg
T_orbit = 2*pi*sqrt((r_orbit^3)/mu_Earth); % s

desired_drop_location_lat = 38.836587; % deg
desired_drop_location_lon = -77.196691; % deg
desired_drop_location_alt = 0; % m; on Earth's surface
latlon_tolerance = 0.05; % deg

alpha_0 = 0; % rad. Assumes ECEF and ECI are aligned at beginning of integration
alpha = 0; % rad
we = (2*pi)/86164; % rad/s

%% Integration Conditions
x0_Earth = [-27590008.8805254;
            143923575.304665;
            19210.6298233792;
            -29.7767856127157;
            -5.53623265806239;
            -0.000194504010588403]*1000; % Earth state at JAN 01 2025 00:00:00.00

% Obtain x0_Sat from orbital elements
a_Sat = r_orbit; % semimajor axis (m)
e_Sat = 0; % eccentricity (Note: controller fails for e > 0.06 since assumption is circ orbit)
i_Sat = 0; % inclination (deg)
O_Sat = 0; % RAAN (deg)
w_Sat = 0; % argument of periapsis (deg)
f_Sat = 0; % true anomaly (deg)

[r_Sat, v_Sat] = elements_to_rv(a_Sat, e_Sat, i_Sat*pi/180, O_Sat*pi/180, w_Sat*pi/180, f_Sat*pi/180, mu_Earth);
x0_Sat = [r_Sat; v_Sat];
dt = 1;
tspan = 0:dt:1000000;

xs_Earth = zeros(6, length(tspan));
xs_SatAct = zeros(6,length(tspan));
xs_SatRef = zeros(6, length(tspan));
us_ECI = zeros(3, length(tspan));
us_LVLH = zeros(3, length(tspan));
xs_LVLH = zeros(6, length(tspan));
alphas = zeros(1, length(tspan));
xs_LLA = zeros(3, length(tspan));
latlon_error = zeros(2, length(tspan));

xs_Earth(:, 1) = x0_Earth;
xs_SatAct(:, 1) = x0_Sat;
xs_SatRef(:, 1) = x0_Sat;
us_ECI(:,1) = zeros(3,1);
us_LVLH(:,1) = zeros(3,1);
xs_LVLH(:,1) = zeros(6,1);
alphas(1) = 0; 
xs_LLA(1) = compute_LLA(r_Sat);
latlon_error(:, 1) = [desired_drop_location_lat; desired_drop_location_lon];

%% Controller Derivation
% Controller acts in the LVLH frame, with states being the pos/vel relative
% to the desired orbit. The desired orbit is integrated in tandem with the
% actual orbit, so the controller can act on the actual orbit and reduce
% its error.
n = sqrt(mu_Earth/(r_orbit^3));
A = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     3*(n^2) 0 0 0 2*n 0;
     0 0 0 -2*n 0 0;
     0 0 -(n^2) 0 0 0];
B = [0 0 0;
     0 0 0;
     0 0 0;
     1/m_Sat 0 0;
     0 1/m_Sat 0;
     0 0 1/m_Sat];

Q = [10 0 0 0 0 0;
     0 10 0 0 0 0;
     0 0 10 0 0 0;
     0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1];

R = [1 0 0;
     0 1 0;
     0 0 1];

K = lqr(A, B, Q, R);



%% Simulation
for i = 1:1:2*T_orbit
    [xE_next, xSatAct_next, xSatRef_next] = rk4_orbit_integ(xs_Earth(:,i), xs_SatAct(:,i), xs_SatRef(:, i), us_ECI(:,i), dt);
    xs_Earth(:,i+1) = xE_next;
    xs_SatAct(:,i+1) = xSatAct_next;
    xs_SatRef(:,i+1) = xSatRef_next;

    [aRef, eRef, iRef, ORef, TRef] = rv_to_elements(xs_SatRef(1:3, i+1), xs_SatRef(4:6, i+1), mu_Earth);
    dx_ECI = xs_SatAct(:,i+1) - xs_SatRef(:, i+1);
    R_ECI_2_LVLH = compute_rotation_ECI_2_LVLH(iRef, ORef, TRef);
    x_LVLH = [R_ECI_2_LVLH*dx_ECI(1:3);
                R_ECI_2_LVLH*dx_ECI(4:6)];
    u_LVLH = -K*x_LVLH;
    
    R_LVLH_2_ECI = inv(R_ECI_2_LVLH);

    us_LVLH(:,i+1) = u_LVLH;
    us_ECI(:,i+1) = R_LVLH_2_ECI*u_LVLH;
    xs_LVLH(:,i+1) = x_LVLH;
    
    alpha = we*tspan(i+1) + alpha_0;
    alphas(i+1) = alpha;

    R_ECI_2_ECEF = compute_R_ECI_2_ECEF(alpha);
    r_ECEF = R_ECI_2_ECEF*xs_SatAct(1:3, i+1);
    [lat, lon, alt] = compute_LLA(r_ECEF); % Returns lat/lon in deg
    xs_LLA(:,i+1) = [lat; lon; alt];
    latlon_error(:, i+1) = [desired_drop_location_lat - lat; desired_drop_location_lon - lon];

    seconds_iterated = i;
end

% Desired inclination
iDesired = desired_drop_location_lat*(pi/180);

% Find DV for actual satellite
vCurrentAct = xs_SatAct(4:6, seconds_iterated+1);
[aAct, eAct, iAct, OAct, TAct] = rv_to_elements(xs_SatAct(1:3, seconds_iterated+1), xs_SatAct(4:6, seconds_iterated+1), mu_Earth);
[rDesiredAct, vDesiredAct] = elements_to_rv(aAct, eAct, iDesired, OAct, TAct, 0, mu_Earth);
DeltaVAct = vDesiredAct - vCurrentAct;

% Find DV for reference satellite
vCurrentRef = xs_SatRef(4:6, seconds_iterated+1);
[aRef, eRef, iRef, ORef, TRef] = rv_to_elements(xs_SatRef(1:3, seconds_iterated+1), xs_SatRef(4:6, seconds_iterated+1), mu_Earth);
[rDesiredRef, vDesiredRef] = elements_to_rv(aRef, eRef, iDesired, ORef, TRef, 0, mu_Earth);
DeltaVRef = vDesiredRef - vCurrentRef;

% Apply DeltaVs
xs_SatAct(4:6, seconds_iterated+1) = vCurrentAct + DeltaVAct;
xs_SatRef(4:6, seconds_iterated+1) = vCurrentRef + DeltaVRef;

% Parameters to determine drop point
hasOrbitedOnce = false;
hasOrbitedTwice = false;
readyToDrop = false;
phaseShiftDetermined = false;
lonAtDesiredLatFirst = 0;
lonAtDesiredLatSecond = 0;
lonPhaseShiftPerOrbit = 0;
latTolerance = 0.01; % deg
lonTolerance = 0.05; % deg
timeOfFirstPass = 0;
timeOfSecondPass = 0;
timeOfFinalPass = 0;
DTheta = 52.0775; % deg; from EDL script
DThetaTolerance = 0.05; % deg
timeToDropFromFinalPass = ((360-DTheta)/360)*T_orbit; % s
hasDeployed = false;
timeToLand = ceil(992.6651); % comes from EDL script


next_start_t = seconds_iterated + 1;
for i = next_start_t:1:(length(tspan)-1)
    [xE_next, xSatAct_next, xSatRef_next] = rk4_orbit_integ(xs_Earth(:,i), xs_SatAct(:,i), xs_SatRef(:, i), us_ECI(:,i), dt);
    xs_Earth(:,i+1) = xE_next;
    xs_SatAct(:,i+1) = xSatAct_next;
    xs_SatRef(:,i+1) = xSatRef_next;

    [aRef, eRef, iRef, ORef, TRef] = rv_to_elements(xs_SatRef(1:3, i+1), xs_SatRef(4:6, i+1), mu_Earth);
    dx_ECI = xs_SatAct(:,i+1) - xs_SatRef(:, i+1);
    R_ECI_2_LVLH = compute_rotation_ECI_2_LVLH(iRef, ORef, TRef);
    x_LVLH = [R_ECI_2_LVLH*dx_ECI(1:3);
                R_ECI_2_LVLH*dx_ECI(4:6)];
    u_LVLH = -K*x_LVLH;
    
    R_LVLH_2_ECI = inv(R_ECI_2_LVLH);

    us_LVLH(:,i+1) = u_LVLH;
    us_ECI(:,i+1) = R_LVLH_2_ECI*u_LVLH;
    xs_LVLH(:,i+1) = x_LVLH;

    alpha = we*tspan(i+1) + alpha_0;
    alphas(i+1) = alpha;

    R_ECI_2_ECEF = compute_R_ECI_2_ECEF(alpha);
    r_ECEF = R_ECI_2_ECEF*xs_SatAct(1:3, i+1);
    [lat, lon, alt] = compute_LLA(r_ECEF); % Returns lat/lon in deg
    xs_LLA(:,i+1) = [lat; lon; alt];
    latlon_error(:, i+1) = [desired_drop_location_lat - lat; desired_drop_location_lon - lon];
    
    if abs(desired_drop_location_lat - lat) <= latTolerance
        if ~hasOrbitedOnce
            lonAtDesiredLatFirst = lon;
            hasOrbitedOnce = true;
            timeOfFirstPass = tspan(i);
        elseif ~hasOrbitedTwice && (tspan(i) - timeOfFirstPass) >= 0.5*T_orbit
            lonAtDesiredLatSecond = lon;
            hasOrbitedTwice = true;
            timeOfSecondPass = tspan(i);
        end
    end

    if hasOrbitedOnce && hasOrbitedTwice && ~phaseShiftDetermined
        lonPhaseShiftPerOrbit = lonAtDesiredLatSecond - lonAtDesiredLatFirst; % deg
        phaseShiftDetermined = true;
    end

    if phaseShiftDetermined && abs(desired_drop_location_lat - lat) <= latTolerance
        if abs(lon + lonPhaseShiftPerOrbit) <= lonTolerance
            readyToDrop = true;
            timeOfFinalPass = tspan(i);
        end
    end

    if readyToDrop && tspan(i) - timeOfFinalPass >= timeToDropFromFinalPass && ~hasDeployed
        hasDeployed = true;
        disp(['Capsule deployed from satellite at t = ', num2str(tspan(i)), ' seconds.'])
        disp(['Deployment Latitude: ', num2str(lat), ' degrees.'])
        disp(['Deployment Longitude: ', num2str(lon), ' degrees.'])
        r_DropPoint_ECI = xs_SatAct(1:3, i+1);
        timeOfDeployment = tspan(i);
    end
    
    if hasDeployed && tspan(i) == timeOfDeployment + timeToLand
        r_Target_ECEF = compute_ECEF_from_LLA(desired_drop_location_lat, desired_drop_location_lon, desired_drop_location_alt);
        R_ECEF_2_ECI = compute_R_ECEF_2_ECI(alpha);
        r_Target_ECI = R_ECEF_2_ECI*r_Target_ECEF;
    end
    % if norm(latlon_error(:,i+1)) <= latlon_tolerance
    %         disp(['Satellite passes within ', ... 
    %             num2str(latlon_tolerance), ' deg of target surface location at ' ...
    %             't = ' num2str(tspan(i+1)) ' seconds.'])
    %         disp(['Error: ', num2str(norm(latlon_error(:,i+1))), ' degrees.'])
    %         disp('')
    % end
    
    seconds_iterated = i;
end


%% Plots
% Heliocentric plot of Earth orbit
figure(1); clf
plot3(xs_Earth(1,:), xs_Earth(2,:), xs_Earth(3,:));
hold on
plot3(0,0,0, 'o', 'Color', 'k', 'MarkerSize', 10, 'MarkerFaceColor', 'y')
xlabel('X Position (m)', 'FontSize', 15)
ylabel('Y Position (m)', 'FontSize', 15)
zlabel('Z Position (m)', 'FontSize', 15)
title('Heliocentric View of Earth Orbit', 'FontSize', 15)

% Geocentric plot of satellite orbit, actual and reference
figure(2); clf
plot3(xs_SatAct(1,:), xs_SatAct(2,:), xs_SatAct(3,:));
hold on
plot3(xs_SatRef(1,:), xs_SatRef(2,:), xs_SatRef(3,:));
plot3(0,0,0, 'o', 'Color', 'k', 'MarkerSize', 10, 'MarkerFaceColor', 'g')
%plot3(r_DropPoint_ECI(1), r_DropPoint_ECI(2), r_DropPoint_ECI(3), 'd', 'Color', 'm', 'MarkerSize', 10)
%plot3(r_Target_ECI(1), r_Target_ECI(2), r_Target_ECI(3), 'd', 'Color', 'c','MarkerSize', 10)
xlabel('X Position (m)', 'FontSize', 15)
ylabel('Y Position (m)', 'FontSize', 15)
zlabel('Z Position (m)', 'FontSize', 15)
legend('Actual Orbit', 'Reference Orbit', 'Earth', 'Location', 'best')
title('Geocentric View of Satellite Orbit', 'FontSize', 15)
axis equal

% Satellite altitude as a function of time
figure(3); clf
plot(tspan, vecnorm(xs_SatAct(1:3, :)) - 6378*1000)
xlabel('Time (s)', 'FontSize', 15)
ylabel('Satellite Altitude', 'FontSize', 15)
title('Satellite Altitude Over Time', 'FontSize', 15)

% Heliocentric plot of both orbits
figure(4); clf
plot3(xs_Earth(1,:), xs_Earth(2,:), xs_Earth(3,:), 'Color', 'g');
hold on
plot3(xs_Earth(1,:) + xs_SatAct(1,:), xs_Earth(2,:) + xs_SatAct(2,:), xs_Earth(3,:) + xs_SatAct(3,:), 'Color', 'k');
plot3(0,0,0, 'o', 'Color', 'k', 'MarkerSize', 10, 'MarkerFaceColor', 'y')
legend('Earth', 'Satellite', 'Sun', 'FontSize', 13, 'Location', 'best')
title('Heliocentric View of Earth and Satellite Orbits', 'FontSize', 15)
%legend('Earth', 'Satellite', 'FontSize', 13, 'Location', 'best')

% LVLH State (Satellite Error State)
figure(5); clf
plot(tspan, xs_LVLH(1,:))
hold on
plot(tspan, xs_LVLH(2,:))
plot(tspan, xs_LVLH(3,:))
plot(tspan, xs_LVLH(4,:))
plot(tspan, xs_LVLH(5,:))
plot(tspan, xs_LVLH(6,:))
xlabel('Time (s)', 'FontSize', 15)
ylabel('Satellite Error States', 'FontSize', 15)
title('Satellite Error States Relative to Desired Reference Orbit', 'FontSize', 15)
legend('e_x', 'e_y', 'e_z', 'e_{vx}', 'e_{vy}', 'e_{vz}', 'FontSize', 15)

% Control Inputs
figure(6); clf
plot(tspan, us_LVLH(1,:))
hold on
plot(tspan, us_LVLH(2,:))
plot(tspan, us_LVLH(3,:))
xlabel('Time (s)', 'FontSize', 15)
ylabel('Control Thruster Forces', 'FontSize', 15)
title('Control Inputs in LVLH Frame', 'FontSize', 15)
legend('F_x', 'F_y', 'F_z', 'FontSize', 13, 'Location', 'best')

% Latitude/Longitude
figure(7); clf
plot(tspan, xs_LLA(1,:))
hold on
plot(tspan, xs_LLA(2,:))
xlabel('Time (s)', 'FontSize', 15)
ylabel('ECEF Position', 'FontSize', 15)
title('Satellite Latitude and Longitude Over Time', 'FontSize', 15)
legend('Latitude', 'Longitude', 'FontSize', 13, 'Location', 'best')

figure(8); clf
plot(tspan, latlon_error(1,:))
hold on
plot(tspan, latlon_error(2,:))
xlabel('Time (s)', 'FontSize', 15)
ylabel('Geodetic Error (deg)', 'FontSize', 15)
title('Satellite Latitude and Longitude Error Over Time', 'FontSize', 15)
legend('Latitude Error', 'Longitude Error', 'FontSize', 13, 'Location', 'best')

function[xdot]= Earth_actual_orbit_prop(x, u, r_S2E)

    % r_S2E = vector from Sun to Earth

    px = x(1);
    py = x(2);
    pz = x(3);
    vx = x(4);
    vy = x(5);
    vz = x(6);

    r = sqrt(px^2 + py^2 + pz^2); % m
    V = sqrt(vx^2 + vy^2 + vz^2); % m/s
    Re = 6378*1000; % m
    J2 = 0.0010827;
    mu = 3.986e14; %m^3/s^2

    J2_accel_vec = [(1-5*(pz/r)^2)*(px/r); 
                    (1-5*(pz/r)^2)*(py/r); 
                    (3-5*(pz/r)^2)*(pz/r)];

    f_J2 = (-3/2)*J2*(mu/(r^2))*((Re/r)^2)*J2_accel_vec;

    h = r - Re;
    rho = 1.225*exp(-h/8420); % kg/m^3
    CD = 2.3;
    A = 4; % m^2
    m = 3000; % kg

    f_drag = -(1/2)*rho*((CD*A)/m)*V*[vx; vy; vz];

    P = 4.56e-6; % solar momentum flux, N/m^2
    theta = acos(dot(r_S2E, [px py pz])/(r*norm(r_S2E)));
    if theta < pi/2 - acos(Re/r)
        v = 0; % light indicator: =0 if in darkness
    else
        v = 1;
    end

    CR = 1; % reflectivity coefficient
    r_S2Sat = [r_S2E(1) + px; r_S2E(2)+py; r_S2E(3)+pz];
    u_S2Sat = r_S2Sat/norm(r_S2Sat);

    f_SRP = -P*((v*A)/m)*CR*u_S2Sat;
    f_grav = (-mu/(r^3))*[px;py;pz];

    pxdot = vx;
    pydot = vy;
    pzdot = vz;
    vxdot = f_grav(1) + f_J2(1) + f_drag(1) + f_SRP(1) + u(1)/m;
    vydot = f_grav(2) + f_J2(2) + f_drag(2) + f_SRP(2) + u(2)/m;
    vzdot = f_grav(3) + f_J2(3) + f_drag(3) + f_SRP(3) + u(3)/m;

    xdot = [pxdot; pydot; pzdot; vxdot; vydot; vzdot];
    
end

function [xdot]= Earth_reference_orbit_prop(x)

    px = x(1);
    py = x(2);
    pz = x(3);
    r = sqrt(px^2 + py^2 + pz^2);

    mu = 3.986e14; %m^3/s^2
    f_grav = (-mu/(r^3))*[px;py;pz];

    xdot = [x(4);
            x(5);
            x(6);
            f_grav(1);
            f_grav(2);
            f_grav(3)];

end

function [xdot]= Sun_orbit_prop(x)

    px = x(1);
    py = x(2);
    pz = x(3);
    r = sqrt(px^2 + py^2 + pz^2);

    mu = 1.327e20; %m^3/s^2
    f_grav = (-mu/(r^3))*[px;py;pz];

    xdot = [x(4);
            x(5);
            x(6);
            f_grav(1);
            f_grav(2);
            f_grav(3)];

end

function[xE_next, xSatAct_next, xSatRef_next]= rk4_orbit_integ(xE, xSatAct, xSatRef, u_ECI, dt)

        k1E = Sun_orbit_prop(xE);
        k2E = Sun_orbit_prop(xE + k1E*(dt/2));
        k3E = Sun_orbit_prop(xE + k2E*(dt/2));
        k4E = Sun_orbit_prop(xE + k3E*dt);

        xE_next = xE + dt*((1/6)*k1E + (1/3)*k2E + (1/3)*k3E + (1/6)*k4E);

        k1SatAct = Earth_actual_orbit_prop(xSatAct, u_ECI, xE(1:3));
        k2SatAct = Earth_actual_orbit_prop(xSatAct + k1SatAct*(dt/2), u_ECI, xE(1:3));
        k3SatAct = Earth_actual_orbit_prop(xSatAct + k2SatAct*(dt/2), u_ECI, xE(1:3));
        k4SatAct = Earth_actual_orbit_prop(xSatAct + k3SatAct*dt, u_ECI, xE(1:3));

        xSatAct_next = xSatAct + dt*((1/6)*k1SatAct + (1/3)*k2SatAct + (1/3)*k3SatAct + (1/6)*k4SatAct);

        k1SatRef = Earth_reference_orbit_prop(xSatRef);
        k2SatRef = Earth_reference_orbit_prop(xSatRef + k1SatRef*(dt/2));
        k3SatRef = Earth_reference_orbit_prop(xSatRef + k2SatRef*(dt/2));
        k4SatRef = Earth_reference_orbit_prop(xSatRef + k3SatRef*dt);

        xSatRef_next = xSatRef + dt*((1/6)*k1SatRef + (1/3)*k2SatRef + (1/3)*k3SatRef + (1/6)*k4SatRef);

end

function[r, v]= elements_to_rv(a, e, i, O, w, f, mu)

    % Input all element angles in radians
    % Input mu in SI units, no km stuff
   
    
    r_norm = (a*(1-e^2))/(1+e*cos(f));
    r_peri = [r_norm*cos(f); r_norm*sin(f); 0];
    
    p = a*(1-e^2);
    v_peri = [sqrt(mu/p)*-1*sin(f); sqrt(mu/p)*(e+cos(f)); 0];
    
    % Rotation matrix terms
    R_11 = cos(O)*cos(w) - sin(O)*sin(w)*cos(i);
    R_12 = -cos(O)*sin(w) - sin(O)*cos(w)*cos(i);
    R_13 = sin(O)*sin(i);
    
    R_21 = sin(O)*cos(w) + cos(O)*sin(w)*cos(i);
    R_22 = -sin(O)*sin(w) + cos(O)*cos(w)*cos(i);
    R_23 = -cos(O)*sin(i);
    
    R_31 = sin(w)*sin(i);
    R_32 = cos(w)*sin(i);
    R_33 = cos(i);
    
    R = [R_11, R_12, R_13;
         R_21, R_22, R_23;
         R_31, R_32, R_33];
    
    r = R*r_peri;
    v = R*v_peri;

end

function[a, e, i, O, T]= rv_to_elements(r, v, mu)
    
    % Input r, v, mu in SI units (m, s)
    % Output a (m), e, i/O/w/f (rads)

    h = cross(r, v);
    eps = (norm(v)^2)/2 - mu/norm(r);

    a = -mu/(2*eps); % semimajor axis, m

    if 1 + 2*eps*norm(h)*norm(h)/(mu^2) <= 0
        e = 0.001;
    else
        e = sqrt(1 + 2*eps*norm(h)*norm(h)/(mu^2)); % eccentricity
    end

    p = norm(h)*norm(h)/mu;
    cosf = (1/(norm(r)*e))*(p-norm(r));
    sinf = (p/(norm(h)*e))*(dot(r, v)/norm(r));
    
    f = atan2(sinf, cosf);
    
    % Find i
    cosi = h(3)/norm(h);
    i = acos(cosi);

    % Find O
    if h(2) == 0 && h(1) == 0
        O = 0;
    else
        sinO = h(1)/sqrt(h(1)^2 + h(2)^2);
        cosO = -h(2)/sqrt(h(1)^2 + h(2)^2);
        O = atan2(sinO, cosO);
    end

    % Find T = w+f
    n_vec = [cos(O);
             sin(O);
             0];
    cosT = dot(n_vec, r)/norm(r);
    T = acos(cosT);
    if dot(r,v) < 0
        T = 2*pi - T;
    end
end

function [R_LVLH_2_ECI]= compute_rotation_LVLH_2_ECI(i, O, T)

    R = zeros(3,3);
    R(1,1) = cos(O)*cos(T)*cos(i)*cos(i) - sin(O)*sin(T)*cos(i) + cos(O)*cos(T)*sin(i)*sin(i);
    R(1,2) = -(cos(O)*sin(T)*cos(i)*cos(i) + cos(T)*sin(O)*cos(i) + cos(O)*sin(T)*sin(i)*sin(i));
    R(1,3) = sin(i)*sin(O);
    R(2,1) = cos(T)*sin(O)*cos(i)*cos(i) + cos(O)*sin(T)*cos(i) + cos(T)*sin(O)*sin(i)*sin(i);
    R(2,2) = -(sin(O)*sin(T)*cos(i)*cos(i) - cos(O)*cos(T)*cos(i) + sin(O)*sin(T)*sin(i)*sin(i));
    R(2,3) = -(cos(O)*sin(i));
    R(3,1) =  sin(i)*sin(T);
    R(3,2) = cos(T)*sin(i);
    R(3,3) = cos(i);

    R_LVLH_2_ECI = R;

end

function [R_ECI_2_LVLH]= compute_rotation_ECI_2_LVLH(i, O, T)

    R = zeros(3,3);
    R(1,1) = cos(O)*cos(T) - sin(O)*sin(T)*cos(i);
    R(1,2) = sin(O)*cos(T) + cos(O)*sin(T)*cos(i);
    R(1,3) = sin(T)*sin(i);
    R(2,1) = -cos(O)*sin(T) - sin(O)*cos(T)*cos(i);
    R(2,2) = -sin(O)*sin(T) + cos(O)*cos(T)*cos(i);
    R(2,3) = cos(T)*sin(i);
    R(3,1) = sin(O)*sin(i);
    R(3,2) = -cos(O)*sin(i);
    R(3,3) = cos(i);

    R_ECI_2_LVLH = R;

end

function[R]= compute_R_ECEF_2_ECI(alpha)   
    R = [cos(alpha), -sin(alpha), 0;
         sin(alpha), cos(alpha), 0;
         0, 0, 1];
end

function[R]= compute_R_ECI_2_ECEF(alpha)
    R = [cos(alpha), sin(alpha), 0;
         -sin(alpha), cos(alpha), 0;
         0, 0, 1];
end

function[lat, lon, alt]= compute_LLA(r_ECEF)   
    lat = asin(r_ECEF(3)/norm(r_ECEF));
    
    coslon = r_ECEF(1)/(norm(r_ECEF)*cos(lat));
    sinlon = r_ECEF(2)/(norm(r_ECEF)*cos(lat));
    lon = atan2(sinlon, coslon);

    lat = lat*(180/pi);
    lon = lon*(180/pi);
    alt = norm(r_ECEF) - 6378*1000;
end

function [r_ECEF] = compute_ECEF_from_LLA(lat, lon, alt)
    lat = lat*(pi/180);
    lon = lon*(pi/180);

    x = (6378*1000 + alt)*cos(lat)*cos(lon);
    y = (6378*1000 + alt)*cos(lat)*sin(lon);
    z = (6378*1000 + alt)*sin(lat);

    r_ECEF = [x; y; z];
end