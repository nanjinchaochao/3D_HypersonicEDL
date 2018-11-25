%%% Driver: Test the dynamics of hypersonic EDL modeled through Vinh's equations
%%% Aa three state model where the dynamics is assumed to be purely longitudinal
%%% is employed
clc; clear; close all;
global N
global R_0;     % mean equatorial radius of Mars
global B_c;     % ballistic coefficient
global v_c;     % normalizing velocity constant
global rho_0;   % reference-level density
global mu_Mars; % standard gravitational parameter of Mars
global C_LD;    % lift-to-drag ratio
global h_1 h_2;
global g;       % acceleration due to gravity

N = 3;

%% System parameters:
R_0 = 3397;                 % km
h_1 = 9.8;                  % km
h_2 = 20;                   % km
B_c = 72.8*10^6;            % kg/km^2
C_LD = 0.3;
rho_0 = 0.0019*10^9;        % kg/km^3;
mu_Mars = 4.282837*10^4;    % km^3/s^2;
g = mu_Mars/R_0^2;          % km/s^2;
v_c = sqrt(g*R_0);          % km/s;
NonD = [R_0 v_c 180/pi];    % nondimensional parameters for three states

% Time
t_f = 0.5;                  % forward propagation for 3 days
% tLEN = length(t_f);
t_step = 0.01;              % fixed time-step numerical propagation of geosynchronous orbits
t_N = R_0/v_c;              % nondimensional parameter for time
NSAM = 2000;

% Initial uncertainty defined on the position and velocity
mu0 = [80 3.5 -2];
% Covar0 = diag(abs(mu0)*0.1); 
Covar0 = diag((abs(mu0)*0.1).^2); 
sam_0 = mvnrnd(mu0, Covar0, NSAM);
w_0 = mvnpdf(sam_0, mu0, Covar0);

% mu0 = [80 3.5 -2];
% a = mu0 - mu0*0.05;
% b = mu0 + mu0*0.05;
% sam_0 = repmat(a, NSAM, 1) + repmat((b - a), NSAM, 1).*rand(NSAM, N);
% w_0 = repmat(1/abs(prod(b - a)), NSAM, 1);

figure(1)
plot3(sam_0(:, 1), sam_0(:, 2), sam_0(:, 3), 'r.')
hold on;
grid on;
xlabel('h [km]'); ylabel('V [km/s]');zlabel('FPA [degrees]');
set(gca,'FontSize',18, 'fontweight','bold')

options = odeset('RelTol', 1e-9, 'AbsTol', 1e-10);
% k = 1;
for ntr = 1 : NSAM
    ic = [sam_0(ntr, :)./NonD w_0(ntr)];
    [t_ode, x_ode] = ode45(@HypersonicEDL3D_SLE, 0 : t_step: t_f, ic, options);
    sam(ntr, :) = x_ode(end, 1 : N).*NonD;
    w(ntr, 1) = x_ode(end, N + 1); 
    [M(ntr, 1) , q(ntr, 1)] = ChuteDeployment(sam(ntr, :));
    Q_dot(ntr, 1) = Heating_Rate(sam(ntr, :));
%     if idx == 1
%         sam_arise0(k, :) = sam_0(ntr, :);
%         w_arise0(k, 1) = w_0(ntr);
%         sam_arise(k, :) = sam(ntr, :);
%         w_arise(k, 1) = w(ntr);
%         M_arise(k, 1) = M(ntr, 1);
%         q_arise(k, 1) = q(ntr, 1);
%         k = k + 1;
%     end
end

figure(2)
plot3(sam(:, 1), sam(:, 2), sam(:, 3), 'r.');
hold on;
grid on;
xlabel('h [km]'); ylabel('V [km/s]');zlabel('FPA [degrees]');
set(gca,'FontSize',18, 'fontweight','bold')

figure(3)
plot(M, q, 'r.');
hold on;
grid on;
xlabel('Mach'); ylabel('Pressume [Pa]');
set(gca,'FontSize',18, 'fontweight','bold')

figure(4)
histogram(Q_dot);

% if ~isempty(M_arise)
%     figure(1)
%     plot3(sam_arise0(:, 1), sam_arise0(:, 2), sam_arise0(:, 3), 'bs');
%     hold on;
% 
%     figure(2)
%     plot3(sam_arise(:, 1), sam_arise(:, 2), sam_arise(:, 3), 'bs');
%     hold on;
% 
%     figure(3)
%     plot(M_arise, q_arise, 'bs');
%     hold on;
% 
%     figure(4)
%     h1 = histogram(w);
%     hold on;
%     h2 = histogram(w_arise);
%     hold on;
%     xlabel('pdf');
%     legend('All','Bad');
%     set(gca,'FontSize',18, 'fontweight','bold');
% end
