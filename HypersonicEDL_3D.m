function xp = HypersonicEDL_3D(t, x)

global N
global R_0;     % mean equatorial radius of Mars
global B_c;     % ballistic coefficient
global v_c;     % normalizing velocity constant
global rho_0;   % reference-level density
global g;       % acceleration due to gravity
global C_LD;    % lift-to-drag ratio
global h_1 h_2;
% global mu_Mars; % Standard gravitational parameter of Mars

xp = zeros(N, 1);

% Martian atmospheric density
rho = rho_0*exp((h_2 - x(1)*R_0)/h_1);
% % acceleration due to gravity
% g = mu_Mars/(R_0 + x(1))^2;
% % normalizing velocity constant
% v_c = sqrt(g*R_0);


% Altidude
xp(1) = x(2)*sin(x(3));
% Velocity
xp(2) = -rho*R_0/(2*B_c)*x(2)^2 - g*R_0/v_c^2*sin(x(3));
% Flight-path angle
xp(3) = rho*R_0/(2*B_c)*C_LD*x(2) + g*R_0/v_c^2*cos(x(3))*(x(2)/(1 + x(1)) - 1/x(2));


end