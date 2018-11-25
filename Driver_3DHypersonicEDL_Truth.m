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


t_f = linspace(0, 0.5, 11);     
tLEN = length(t_f);
t_step = 0.01;            
t_N = R_0/v_c;    
NSAM = 500000;

% Initial uncertainty defined on the position and velocity
mu0 = [80 3.5 -2];
% Covar0 = diag(abs(mu0)*0.1); 
Covar0 = diag((abs(mu0)*0.1).^2); 
isample = mvnrnd(mu0, Covar0, NSAM);
w_0 = mvnpdf(isample, mu0, Covar0);

for i = 1 : NSAM
    isample_ND(i, :) = isample(i, :)./NonD;
    [isample_M(i, 1), isample_q(i, 1)] = ChuteDeployment(isample(i, :));
    isample_Q(i, 1) = Heating_Rate(isample(i, :));
    
end


if t_f(1) == 0
    MCSAM(1).sam_T = isample;
    MCSAM(1).samND_T = isample_ND;
    MCSAM(1).w_T = w_0;
    MCSAM(1).M_T = isample_M;
    MCSAM(1).P_T = isample_q;
    MCSAM(1).Q_T = isample_Q;
    
    Mean_T(1, :) = mean(MCSAM(1).sam_T);
    Std_T(1, :) = std(MCSAM(1).sam_T);
    Moment3_T(1, :) = moment(MCSAM(1).sam_T, 3);
    Moment4_T(1, :) = moment(MCSAM(1).sam_T, 4);
    
    M_Mean_T(1, :) = mean(MCSAM(1).M_T);
    M_Std_T(1, :) = std(MCSAM(1).M_T);
    M_Moment3_T(1, :) = moment(MCSAM(1).M_T, 3);
    M_Moment4_T(1, :) = moment(MCSAM(1).M_T, 4);
    
    P_Mean_T(1, :) = mean(MCSAM(1).P_T);
    P_Std_T(1, :) = std(MCSAM(1).P_T);
    P_Moment3_T(1, :) = moment(MCSAM(1).P_T, 3);
    P_Moment4_T(1, :) = moment(MCSAM(1).P_T, 4);
    
    Q_Mean_T(1, :) = mean(MCSAM(1).Q_T);
    Q_Std_T(1, :) = std(MCSAM(1).Q_T);
    Q_Moment3_T(1, :) = moment(MCSAM(1).Q_T, 3);
    Q_Moment4_T(1, :) = moment(MCSAM(1).Q_T, 4);
    
    tctstart = 2;
    
else
    tctstart = 1;
end

% mu0 = [80 3.5 -2];
% a = mu0 - mu0*0.05;
% b = mu0 + mu0*0.05;
% sam_0 = repmat(a, NSAM, 1) + repmat((b - a), NSAM, 1).*rand(NSAM, N);
% w_0 = repmat(1/abs(prod(b - a)), NSAM, 1);

% figure(1)
% plot3(sam_0(:, 1), sam_0(:, 2), sam_0(:, 3), 'r.')
% hold on;
% grid on;
% xlabel('h [km]'); ylabel('V [km/s]');zlabel('FPA [degrees]');
% set(gca,'FontSize',18, 'fontweight','bold')

options = odeset('RelTol', 1e-9, 'AbsTol', 1e-10);

for tctr = tctstart : tLEN
    fprintf('Currently running %d/%d time instance, %f ahead \n', tctr, tLEN, t_f(tctr));
    %%%
    % Forward Propagation (make it as a module, and how to generalize it when noise is involved?)
    %%%
    for ctr = 1 : NSAM
        ic = [MCSAM(tctr - 1).samND_T(ctr,:) MCSAM(tctr - 1).w_T(ctr,1)];                                            % Pick up the particle for forward propagation
        [t_ode, x_ode] = ode45(@HypersonicEDL3D_SLE, t_f(tctr - 1): t_step: t_f(tctr), ic, options); 
        MCSAM(tctr).samND_T(ctr, :) = x_ode(end, 1 : N);                                                                 % Collect the propagated particle
        MCSAM(tctr).sam_T(ctr, :) = MCSAM(tctr).samND_T(ctr, :).*NonD;                 % Collect its associated pdf value at the current time
        MCSAM(tctr).w_T(ctr, 1) = x_ode(end, N + 1);
        [MCSAM(tctr).M_T(ctr, 1) , MCSAM(tctr).P_T(ctr, 1)] = ChuteDeployment(MCSAM(tctr).sam_T(ctr, :));
        MCSAM(tctr).Q_T(ctr, 1) = Heating_Rate(MCSAM(tctr).sam_T(ctr, :));
    end
    
    
    
    % Moment estimation
    Mean_T(tctr, :) = mean(MCSAM(tctr).sam_T);
    Std_T(tctr, :) = std(MCSAM(tctr).sam_T);
    Moment3_T(tctr, :) = moment(MCSAM(tctr).sam_T, 3);
    Moment4_T(tctr, :) = moment(MCSAM(tctr).sam_T, 4);
    
    M_Mean_T(tctr, :) = mean(MCSAM(tctr).M_T);
    M_Std_T(tctr, :) = std(MCSAM(tctr).M_T);
    M_Moment3_T(tctr, :) = moment(MCSAM(tctr).M_T, 3);
    M_Moment4_T(tctr, :) = moment(MCSAM(tctr).M_T, 4);
    
    P_Mean_T(tctr, :) = mean(MCSAM(tctr).P_T);
    P_Std_T(tctr, :) = std(MCSAM(tctr).P_T);
    P_Moment3_T(tctr, :) = moment(MCSAM(tctr).P_T, 3);
    P_Moment4_T(tctr, :) = moment(MCSAM(tctr).P_T, 4);
    
    Q_Mean_T(tctr, :) = mean(MCSAM(tctr).Q_T);
    Q_Std_T(tctr, :) = std(MCSAM(tctr).Q_T);
    Q_Moment3_T(tctr, :) = moment(MCSAM(tctr).Q_T, 3);
    Q_Moment4_T(tctr, :) = moment(MCSAM(tctr).Q_T, 4);
end


