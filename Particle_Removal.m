%{
Particle-removal scheme is applied when the transient performance of the
current ensemble is under the prescribed lower error bound. Then, Particles 
are removed (halted) in the interest of alleviating computational load. 
Removal is based on ranking particles by their relative particle weightage, 
evaluated by solving the associated stochastic Liouville equation via the
method of characteristics.
%}

function [sam, sam_ND, QoI, w, U, Acc] = Particle_Removal(sam_Old, sam_ND_Old, QoI_Old, w_Old, U_Old, Acc_Old, Current_time)
global Acc_LB P_Acc;
global Num_resam_bootstrap;
global P_removal; 
global Num_removal;

sam = sam_Old;
sam_ND = sam_ND_Old;
w = w_Old;
Acc = Acc_Old;
U = U_Old;
QoI = QoI_Old;

%% Particle-Removal scheme
while Acc < Acc_LB
     % Estimate the number of samples needed to be remained
     num = ceil((Acc^2/P_Acc^2)*size(sam, 1));        
     % Normlization 
     w_Norm = w./sum(w);
     % Resample num samples in terms of their weights, w_Norm, without replacement
     [sam, idx]  = datasample(sam, num, 'Replace', false, 'Weights', w_Norm);  
     w = w(idx);   
     U = U(idx, :);
     sam_ND = sam_ND(idx, :);
     QoI = QoI(idx);
     % Accuracy Estimation: Boostrap mean 
     Acc = BootstrapMean(QoI, Num_resam_bootstrap);
    
end

Num_removal = Num_removal + 1;        %Store into a new field of the structure 
%% Find the "halted" samples 
Lia = ismember(w_Old, w, 'rows');
id = find(Lia == 0);
P_removal(Num_removal).sam = sam_Old(id, :);
P_removal(Num_removal).sam_ND = sam_ND_Old(id, :);
P_removal(Num_removal).w = w_Old(id);
P_removal(Num_removal).T = Current_time;
P_removal(Num_removal).U = U_Old(id, :);

end