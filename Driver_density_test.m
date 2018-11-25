%%% Driver: test the density

clear; close; clc;


h_1 = 9.8;                  % km
h_2 = 20;                   % km
R_0 = 3397;                 % km
rho_0 = 0.0019*10^9;        % kg/km^3;

h = linspace(0, 30, 100);
for i = 1 : 100
    rho(i) = rho_0*exp((h_2 - h(i))/h_1);
end

figure(1)
semilogx(rho.*10^-9, h, 'ro-');
hold on;