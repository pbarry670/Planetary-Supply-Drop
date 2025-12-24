syms xPrior1 xPrior2 xPrior3 xPrior4 xPrior5 xPrior6 mu dt

F = [xPrior4;
     xPrior5;
     xPrior6;
     -mu*xPrior1/((sqrt(xPrior1^2 + xPrior2^2 + xPrior3^2))^3);
     -mu*xPrior2/((sqrt(xPrior1^2 + xPrior2^2 + xPrior3^2))^3);
     -mu*xPrior3/((sqrt(xPrior1^2 + xPrior2^2 + xPrior3^2))^3)];

A = jacobian(F, [xPrior1 xPrior2 xPrior3 xPrior4 xPrior5 xPrior6])


syms station_r1 station_r2 station_r3 station_v1 station_v2 station_v3

% nonlinear measurement model: range and range-rate
G = [sqrt((xPrior1-station_r1)^2 + (xPrior2-station_r2)^2 + (xPrior3-station_r3)^2);
     ( (xPrior1-station_r1)*(xPrior4-station_v1) + (xPrior2-station_r2)*(xPrior5-station_v2) + (xPrior3-station_r3)*(xPrior6-station_v3) )/(sqrt((xPrior1-station_r1)^2 + (xPrior2-station_r2)^2 + (xPrior3-station_r3)^2))];

H = jacobian(G, [xPrior1 xPrior2 xPrior3 xPrior4 xPrior5 xPrior6])