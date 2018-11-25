%{
Bootstrap Performance Estimation: approximate the transient performance of 
the current ensemble, defined by the RMS value of the standard deviation of 
the MC estimation error in terms of the QoI over all the states.
%}

%% The propagated mean is used as the QoI ()
function Acc = BootstrapMean(sam, Num_resam_bootstrap)
global N;

% Sample mean
Mean_sam = mean(sam);

% Means obtained from the bootstrap resampling sets (how to generalize it to arbitrary QoI)
[bootstat, ~] = bootstrp(Num_resam_bootstrap, @mean, sam);

% The Variance of the bootstrap distribution used to approximate to the variance 
% of MC estimation error distribution in terms of the mean
Var_estimation_error = sum((bsxfun(@minus, bootstat, Mean_sam)).^2, 1)/Num_resam_bootstrap;

% The transient performance: the RMS vaule over all states for the SD of the 
% MC estimation error on the mean
Acc = sqrt(sum(Var_estimation_error)/N);

end