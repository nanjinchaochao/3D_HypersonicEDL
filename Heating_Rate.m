function Q_dot = Heating_Rate(x)
% MSL (Design): peak 200(W/cm^2) and Diameter 4.5m Highest aerothermodynamic
% (margined at ~250 W/cm2)
global h_1 h_2;
global rho_0;
C_f = 1;         % the skin-friction coefficient of the exterior surface area of the entry capsule
D = 2.5;         % Diameter of aeroshell 
S = (D/2)^2*pi;           % m^2: exterior surface area
h_real = x(1);
V_real = x(2);
V = V_real*1000; % m/s

rho = rho_0*exp((h_2 - h_real)/h_1)*10^-9;      % kg/m^3;
Q_dot = C_f*rho*V^3*S/4/10000;

% % Sutton Graves Equation
% k= 1.9027e-4;
% R_n = D/4;          % bluntness
% Q_dot2 = k*sqrt(rho/R_n)*V^3/10000;