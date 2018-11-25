%{
Particle-addition scheme is applied when the transient performance of the
current ensemble exceeds the prescribed upper error bound. Then, Optimally
selected particles are sequentially introduced in batches to the initial 
ensemble, and then forward propagated to join the current ensemble, until
it reaches the prescribed level.
%}

function [Sam, Sam_ND, QoI, w, U, Acc] = Particle_Addition...
    (Sam_Old, Sam_ND_Old, QoI_Old, w_Old, U_Old, Acc_Old, Current_time)
global Acc_UB P_Acc;
global N;
global Num_resam_bootstrap;
global Accuracy_check_frequency;
global NonD;
global P_removal; 
global Num_removal;
global mu0 Covar0;
global t_step;
global options;

Sam = Sam_Old;
Sam_ND = Sam_ND_Old;
QoI = QoI_Old;
w = w_Old;
U = U_Old;
Acc = Acc_Old;
%% At the initial time, i.e., Current_time = 0
if Current_time == 0
   
       while Acc > P_Acc
           % Check the transient performance every Accuracy_check_frequency times
           for i = 1 : Accuracy_check_frequency
               New_U = ParticleGeneration_SamplingEfficiency(U);
               U = [U; New_U];
               New_sam = icdf('Normal', New_U, mu0, sqrt(diag(Covar0))');
               Sam = [Sam; New_sam];
               New_sam_ND = New_sam./NonD;
               Sam_ND = [Sam_ND; New_sam_ND];
               New_w = mvnpdf(New_sam, mu0, Covar0);
               w = [w; New_w];
               [M , ~] = ChuteDeployment(New_sam);
               QoI = [QoI; M];
%                if isnan(M)
%                     aaa = 1;
%                end
           end

          % Accuracy Estimation: Boostrap 
          Acc = BootstrapMean(QoI, Num_resam_bootstrap);
       end
else
%% After forward propagation, i.e., Current_time > 0    
    while Acc > Acc_UB
        % Check if there is any previously "halted" particles (how to make it better?)
        if Num_removal ~= 0
            % Adding the most recently "halted" particles
            for k = 1 : size(P_removal(Num_removal).sam_ND, 1)                 % # of "halted" samples                                                                   
                ic = [P_removal(Num_removal).sam_ND(k, :) P_removal(Num_removal).w(k)];
                [te, x] = ode45(@HypersonicEDL3D_SLE, P_removal(Num_removal).T : t_step: Current_time, ic, options);
                Sam_ND = [Sam_ND; x(end, 1 : N)];
                Sam = [Sam; x(end, 1 : N).*NonD];
                w = [w; x(end, N + 1)];
                U = [U; P_removal(Num_removal).U(k, :)];
                [M , ~] = ChuteDeployment(Sam(end, :));
                QoI = [QoI M];
                % Check if the current performance needs to be estimated
                if mod(k, Accuracy_check_frequency) == 0
                  Acc = BootstrapMean(QoI, Num_resam_bootstrap); 
                  % The prescribed accuracy satisfied
                  if Acc <= Acc_UB
                    % Reset: collect the un-used "halted" particles
                    P_removal(Num_removal).w = P_removal(Num_removal).w(k + 1 : end);
                    P_removal(Num_removal).sam_ND = P_removal(Num_removal).sam_ND(k + 1 : end, :);
                    P_removal(Num_removal).sam = P_removal(Num_removal).sam(k + 1 : end, :);
                    P_removal(Num_removal).U  = P_removal(Num_removal).U(k + 1 : end, :);
                    return;                                                 % jump out of the function
                  end  
                end  
            end
            % Reset: all the "halted" particles at time Num_removal are used
            P_removal(Num_removal) = [];                                     % delete all the "halted" particles at time Num_removal
            Num_removal = Num_removal - 1;                                   % Go to another time instant where the particle-removal has been applied  
            Acc = BootstrapMean(QoI, Num_resam_bootstrap);  
            
        % No previous "halted" particles available, then generate new samples
        else
            for i = 1 : Accuracy_check_frequency
                New_U = ParticleGeneration_SamplingEfficiency(U);
                U = [U; New_U];
                New_sam_0 = icdf('Normal', New_U, mu0, sqrt(diag(Covar0))');
                New_w_0 = mvnpdf(New_sam_0, mu0, Covar0);
                New_sam_ND_0 = New_sam_0./NonD;
                % Forward Propagation
                ic = [New_sam_ND_0 New_w_0];
                [t_ode, x_ode] = ode45(@HypersonicEDL3D_SLE, 0 : t_step: Current_time, ic, options);
                New_sam_ND_t = x_ode(end, 1 : N);
                Sam_ND = [Sam_ND; New_sam_ND_t];
                New_sam_t = New_sam_ND_t.*NonD;
                Sam = [Sam; New_sam_t];
                New_w_t = x_ode(end, N + 1);
                w = [w; New_w_t];
                [M , ~] = ChuteDeployment(New_sam_t);
                QoI = [QoI; M];
%                 if isnan(M)
%                     aaa = 1;
%                 end
                
            end
            Acc = BootstrapMean(QoI, Num_resam_bootstrap);
        end    
    end
    
   
end