% Driver: Adaptive FMC for uncertainty characterization in 3D Mars EDL.

clear; clc; close;

global N;
global R_0;     % mean equatorial radius of Mars
global B_c;     % ballistic coefficient
global v_c;     % normalizing velocity constant
global rho_0;   % reference-level density
global g;       % acceleration due to gravity
global C_LD;    % lift-to-drag ratio
global h_1 h_2;
global NonD;
global t_step;

%% System parameters:
N = 3;
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



%% Adaptive MC parameters:
global t_f;
global Num_resam_bootstrap;
global Accuracy_check_frequency;
global Acc_UB Acc_LB P_Acc;
% global Acc_UB_q Acc_LB_q Acc_q;
global alpha;
global Num_candidate;
global P_removal; 
global mu0 Covar0;
global Num_removal;
global options;

alpha = 0.5;                                  
Num_candidate = 5000; 
P_removal = struct();                               
Num_removal = 0;
Num_resam_bootstrap = 5000;                                                               
Accuracy_check_frequency = 100;
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-8);

% Accuracy level for estimating Mach number and dynamics pressure
%%% Mach: 1.1 - 2.2
% P_Acc = 0.05;
% Acc_UB = 0.1;       
% Acc_LB = 0.01;
P_Acc = 0.005;
Acc_UB = 0.01;       
Acc_LB = 0.001;

%%% Dynamics pressure
% P_Acc = 0.5;
% Acc_UB = 0.1;       
% Acc_LB = 0.01;
% % P_Acc = 0.05;
% % Acc_UB = 0.01;       
% % Acc_LB = 0.001;


% Time
t_f = linspace(0, 0.5, 11);     
tLEN = length(t_f);
t_step = 0.01;            
t_N = R_0/v_c;            


% Initial uncertainty defined on the position and velocity
mu0 = [80 3.5 -2];
Covar0 = diag((abs(mu0)*0.1).^2); 


Num_start = 10;                                     % # of particles to start the simulation 
                                                    % (user-defined input or set as a default, e.g., 10)
a = zeros(N, 1);                            
b = ones(N, 1);  
U = lhsu(a, b, Num_start);                          % Generate uniform distributed particles by LHS

% Inverse transform:
for num = 1 : Num_start
    isample(num, :) = icdf('Normal', U(num, :), mu0, sqrt(diag(Covar0))');
    [QoI(num, 1) , ~] = ChuteDeployment(isample(num, :));
    isample_ND(num, :) = isample(num, :)./NonD;
end

% Accuracy Estimation: Boostrap 
Acc_cur = BootstrapMean(QoI, Num_resam_bootstrap);  
% Acc2 = BootstrapApproximation(QoI, isample);
Acc_all_Old(1, 1) = Acc_cur;   


if t_f(1) == 0
    MCSAM(1).U = U;
    MCSAM(1).sam = isample;
    MCSAM(1).sam_ND = isample_ND;
    MCSAM(1).w = mvnpdf(isample, mu0, Covar0);
    MCSAM(1).QoI = QoI;
    MCSAM(1).Acc = Acc_cur;
    
    tctstart = 2;
    
else
    tctstart = 1;
end


%% Particle-addition scheme is applied until its performance reaches the desired level at the initial time 
[MCSAM(1).sam, MCSAM(1).sam_ND, MCSAM(1).QoI, MCSAM(1).w, MCSAM(1).U, MCSAM(1).Acc] ...
= Particle_Addition(MCSAM(1).sam, MCSAM(1).sam_ND, MCSAM(1).QoI, MCSAM(1).w, MCSAM(1).U, MCSAM(1).Acc, t_f(1));

Num_current_ensemble(1, 1) = size(MCSAM(1).sam, 1);              % # of the current ensemble
Acc_all_New(1, 1) = MCSAM(1).Acc;                            % The transient performance of the current ensemble 

% % Moment estimation
% Mean(1, :) = mean(MCSAM(1).sam);
% Std(1, :) = std(MCSAM(1).sam);
% Moment3(1, :) = moment(MCSAM(1).sam, 3);
% Moment4(1, :) = moment(MCSAM(1).sam, 4);

% %% Plot
% figure(1)
% plot(t_f*t_N/60, repmat(Acc_UB, tLEN, 1), 'r-', 'linewidth', 2);
% hold on;
% plot(t_f*t_N/60, repmat(Acc_LB, tLEN, 1), 'r-', 'linewidth', 2);
% hold on;
% plot(t_f(1)*t_N/60, MCSAM(1).Acc, 'kd','MarkerSize', 8, 'linewidth', 2);
% hold on;
% xlabel('Time [min]'); ylabel('Error');
tic;
%% Adaptive FMC
for tctr = tctstart : tLEN
%     fprintf('Currently running %d/%d time instance, %f ahead \n', tctr, tLEN, t_f(tctr));
    %%%
    % Forward Propagation (make it as a module, and how to generalize it when noise is involved?)
    %%%
    for ctr = 1 : size(MCSAM(tctr - 1).sam, 1)
        ic = [MCSAM(tctr - 1).sam_ND(ctr,:) MCSAM(tctr - 1).w(ctr,1)];                                            % Pick up the particle for forward propagation
        [t_ode, x_ode] = ode45(@HypersonicEDL3D_SLE, t_f(tctr - 1): t_step: t_f(tctr), ic, options); 
        MCSAM(tctr).sam_ND(ctr, :) = x_ode(end, 1 : N);                                                                 % Collect the propagated particle
        MCSAM(tctr).sam(ctr, :) = MCSAM(tctr).sam_ND(ctr, :).*NonD;                 % Collect its associated pdf value at the current time
        MCSAM(tctr).w(ctr, 1) = x_ode(end, N + 1);
        [MCSAM(tctr).QoI(ctr, 1) , ~] = ChuteDeployment(MCSAM(tctr).sam(ctr, :));
%         if isnan(MCSAM(tctr).QoI(ctr, 1))
%             aaa = 1;
%         end
    end
    
    %%% 
    % Accuracy Estimation: Boostrap 
    %%%  
    Acc_cur = BootstrapMean(MCSAM(tctr).QoI, Num_resam_bootstrap);
    Acc_all_Old(tctr, 1) = Acc_cur;
%     figure(1)
%     plot(t_f(tctr)*t_N/60, Acc_cur, 'bo', 'MarkerSize', 8, 'linewidth', 2);
%     hold on;
    
    %%%   
    % Adaptive scheme
    %%%  
    % Particle_Addition when the performance of the current ensemble is greater than the prescribed upper error bound
    % (make it as a module)
    if Acc_cur > Acc_UB
        [MCSAM(tctr).sam, MCSAM(tctr).sam_ND, MCSAM(tctr).QoI, MCSAM(tctr).w, MCSAM(tctr).U, MCSAM(tctr).Acc] ...
         = Particle_Addition(MCSAM(tctr).sam, MCSAM(tctr).sam_ND, MCSAM(tctr).QoI, MCSAM(tctr).w, MCSAM(tctr - 1).U, Acc_cur, t_f(tctr));
    
    % Particle_Removal when the performance of the current ensemble is smaller than the prescribed lower error bound
    % (make it as a module, and how to let the user have the flexibility to select this?)
    elseif Acc_cur < Acc_LB
        [MCSAM(tctr).sam, MCSAM(tctr).sam_ND, MCSAM(tctr).w, MCSAM(tctr).U, MCSAM(tctr).Acc]...
         = Particle_Removal(MCSAM(tctr).sam, MCSAM(tctr).sam_ND, MCSAM(tctr).w, MCSAM(tctr - 1).U, Acc_cur, t_f(tctr));
    
    % Within the prescribed accuracy bounds
    else
        MCSAM(tctr).Acc = Acc_cur;
        MCSAM(tctr).U = MCSAM(tctr - 1).U;
    end
    
    Num_current_ensemble(tctr, 1) = size(MCSAM(tctr).sam, 1);
    Acc_all_New(tctr, 1) = MCSAM(tctr).Acc;
%     figure(1)
%     plot(t_f(tctr)*t_N/60, MCSAM(tctr).Acc, 'kd', 'MarkerSize', 8, 'linewidth', 2);
%     hold on;
    
%     % Moment estimation
%     Mean(tctr, :) = mean(MCSAM(tctr).sam);
%     Std(tctr, :) = std(MCSAM(tctr).sam);
%     Moment3(tctr, :) = moment(MCSAM(tctr).sam, 3);
%     Moment4(tctr, :) = moment(MCSAM(tctr).sam, 4);
end
t_AMC = toc;
fprintf('Adaption process finished \n');
%% Generate accuracy histroy in terms of min and max number of samples
tic;
Num_min = min(Num_current_ensemble(:));
% MCSAM(1).pseudo_min = mvnrnd(mu0, Covar0, Num_min);
% MCSAM(1).pseudo_min_w = mvnpdf(MCSAM(1).pseudo_min, mu0, Covar0);
MCSAM(1).pseudo_min_isample = rand([Num_min, N]);
for ntr = 1 : Num_min
    MCSAM(1).pseudo_min(ntr, :) = icdf('Normal', MCSAM(1).pseudo_min_isample(ntr, :), mu0, sqrt(diag(Covar0))');
    MCSAM(1).pseudo_min_w(ntr, :) = mvnpdf(MCSAM(1).pseudo_min(ntr, :), mu0, Covar0);
    
    [MCSAM(1).pseudo_min_QoI(ntr, :), ~] = ChuteDeployment(MCSAM(1).pseudo_min(ntr, :));
    MCSAM(1).pseudo_min_ND(ntr, :) = MCSAM(1).pseudo_min(ntr, :)./NonD;
end
Acc_pseudo_min(1) = BootstrapMean(MCSAM(1).pseudo_min_QoI, Num_resam_bootstrap);
% % Moment estimation
% Mean_min(1, :) = mean(MCSAM(1).pseudo_min);
% Std_min(1, :) = std(MCSAM(1).pseudo_min);
% Moment3_min(1, :) = moment(MCSAM(1).pseudo_min, 3);
% Moment4_min(1, :) = moment(MCSAM(1).pseudo_min, 4);

for tctr = tctstart : tLEN
        for ctr = 1 : Num_min
            ic = [MCSAM(tctr - 1).pseudo_min_ND(ctr,:) MCSAM(tctr - 1).pseudo_min_w(ctr,1)];  
            [t_ode, x_ode] = ode45(@HypersonicEDL3D_SLE, t_f(tctr - 1): t_step: t_f(tctr), ic, options);
            MCSAM(tctr).pseudo_min_ND(ctr,:) = x_ode(end, 1 : N);  
            MCSAM(tctr).pseudo_min_w(ctr, 1) = x_ode(end, N + 1); 
            MCSAM(tctr).pseudo_min(ctr, :) = MCSAM(tctr).pseudo_min_ND(ctr,:).*NonD;
            [MCSAM(tctr).pseudo_min_QoI(ctr, :), ~] = ChuteDeployment(MCSAM(tctr).pseudo_min(ctr, :));
        end
%         MCSAM(tctr).pseudo_minN = MCSAM(tctr).pseudo_min./repmat(NonDim, Num_min, 1);
%         bootstat_pseudo_min = bootstrp(nboot, @mean, MCSAM(tctr).pseudo_minN);
%         MeanSD_pseudo_min = std(bootstat_pseudo_min)/sqrt(Num_min);
%         Acc_pseudo_min(tctr) = sqrt(sum(MeanSD_pseudo_min.^2)/N);
        Acc_pseudo_min(tctr) = BootstrapMean(MCSAM(tctr).pseudo_min_QoI, Num_resam_bootstrap);
%         % Moment estimation
%         Mean_min(tctr, :) = mean(MCSAM(tctr).pseudo_min);
%         Std_min(tctr, :) = std(MCSAM(tctr).pseudo_min);
%         Moment3_min(tctr, :) = moment(MCSAM(tctr).pseudo_min, 3);
%         Moment4_min(tctr, :) = moment(MCSAM(tctr).pseudo_min, 4);
    
end
t_min = toc;

tic;
Num_max = max(Num_current_ensemble(:));
% MCSAM(1).pseudo_max = mvnrnd(mu0, Covar0, Num_max);
% MCSAM(1).pseudo_max_w = mvnpdf(MCSAM(1).pseudo_max, mu0, Covar0);
MCSAM(1).pseudo_max_isample = rand([Num_max, N]);
for ntr = 1 : Num_max
    MCSAM(1).pseudo_max(ntr, :) = icdf('Normal', MCSAM(1).pseudo_max_isample(ntr, :), mu0, sqrt(diag(Covar0))');
    MCSAM(1).pseudo_max_w(ntr, :) = mvnpdf(MCSAM(1).pseudo_max(ntr, :), mu0, Covar0);
    
    [MCSAM(1).pseudo_max_QoI(ntr, :), ~] = ChuteDeployment(MCSAM(1).pseudo_max(ntr, :));
    MCSAM(1).pseudo_max_ND(ntr, :) = MCSAM(1).pseudo_max(ntr, :)./NonD;
end
Acc_pseudo_max(1) = BootstrapMean(MCSAM(1).pseudo_max_QoI, Num_resam_bootstrap);
% Moment estimation
% Mean_max(1, :) = mean(MCSAM(1).pseudo_max);
% Std_max(1, :) = std(MCSAM(1).pseudo_max);
% Moment3_max(1, :) = moment(MCSAM(1).pseudo_max, 3);
% Moment4_max(1, :) = moment(MCSAM(1).pseudo_max, 4);

for tctr = tctstart : tLEN
        for ctr = 1 : Num_max
            ic = [MCSAM(tctr - 1).pseudo_max_ND(ctr,:) MCSAM(tctr - 1).pseudo_max_w(ctr,1)];  
            [t_ode, x_ode] = ode45(@HypersonicEDL3D_SLE, t_f(tctr - 1): t_step: t_f(tctr), ic, options);
            MCSAM(tctr).pseudo_max_ND(ctr,:) = x_ode(end, 1 : N);  
            MCSAM(tctr).pseudo_max_w(ctr, 1) = x_ode(end, N + 1); 
            MCSAM(tctr).pseudo_max(ctr, :) = MCSAM(tctr).pseudo_max_ND(ctr,:).*NonD;
            [MCSAM(tctr).pseudo_max_QoI(ctr, :), ~] = ChuteDeployment(MCSAM(tctr).pseudo_max(ctr, :));
        end
%         MCSAM(tctr).pseudo_maxN = MCSAM(tctr).pseudo_max./repmat(NonDim, Num_max , 1);
%         bootstat_pseudo_max = bootstrp(nboot, @mean, MCSAM(tctr).pseudo_maxN);
%         MeanSD_pseudo_max = std(bootstat_pseudo_max)/sqrt(Num_max);
%         Acc_pseudo_max(tctr) = sqrt(sum(MeanSD_pseudo_max.^2)/N);
        Acc_pseudo_max(tctr) = BootstrapMean(MCSAM(tctr).pseudo_max_QoI, Num_resam_bootstrap);
%         % Moment estimation
%         Mean_max(tctr, :) = mean(MCSAM(tctr).pseudo_max);
%         Std_max(tctr, :) = std(MCSAM(tctr).pseudo_max);
%         Moment3_max(tctr, :) = moment(MCSAM(tctr).pseudo_max, 3);
%         Moment4_max(tctr, :) = moment(MCSAM(tctr).pseudo_max, 4);
    
end
t_max = toc;

% figure(1)
% plot(t_f*t_N./60, Acc_pseudo_min,'g^--', 'MarkerSize', 8, 'linewidth', 2);
% hold on;
% plot(t_f*t_N./60, Acc_pseudo_max,'mv--', 'MarkerSize', 8, 'linewidth', 2);
% set(gca,'FontSize',18, 'fontweight','bold');
