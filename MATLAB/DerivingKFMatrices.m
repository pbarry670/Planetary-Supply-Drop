syms x y z xdot ydot zdot mu dt

F = [xdot;
     ydot;
     zdot;
     -mu*x/((sqrt(x^2 + y^2 + z^2))^3);
     -mu*y/((sqrt(x^2 + y^2 + z^2))^3);
     -mu*z/((sqrt(x^2 + y^2 + z^2))^3)];

A = jacobian(F, [x y z xdot ydot zdot]);


syms xs ys zs xsdot ysdot zsdot

% nonlinear measurement model: range and range-rate
G = [sqrt((x-xs)^2 + (y-ys)^2 + (z-zs)^2);
     ( (x-xs)*(xdot-xsdot) + (y-ys)*(ydot-ysdot) + (z-zs)*(zdot-zsdot) )/(sqrt((x-xs)^2 + (y-ys)^2 + (z-zs)^2))];

H = jacobian(G, [x y z xdot ydot zdot]);