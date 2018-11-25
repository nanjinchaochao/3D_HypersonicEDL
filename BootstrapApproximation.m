function Acc = BootstrapApproximation(QoI, sam)
global N;
global Num_resam_bootstrap;

% Mean of QoI
Mean_QoI = mean(QoI);
NSAM = length(QoI);
for i = 1 : Num_resam_bootstrap
    idx = randi([1, NSAM], NSAM, 1);
    Bp_QoI = zeros(NSAM, 1);
    for j = 1 : NSAM
        Bp_QoI(j) =  Mach(sam(idx(j), :)); 
    end
    Bp_QoI_Mean(i) = mean(Bp_QoI);
end

% % Means obtained from the bootstrap resampling sets (how to generalize it to arbitrary QoI)
% [bootstat, bootsam] = bootstrp(Num_resam_bootstrap, @(samples) Mach(samples), sam);

Var_estimation_error = sum((bsxfun(@minus, Bp_QoI_Mean, Mean_QoI)).^2, 1)/Num_resam_bootstrap;
Acc = sqrt(sum(Var_estimation_error)/N);



    function M = Mach(x)
    R = 8.3145;         % universal gas constant: J/mol*K
    M = 44.01*10^-3;    % molar mass of the gas
    R_s = R/M;          % specific gas constant: J/kg*K
    h_real = x(1);
    V_real = x(2);
    % FPA = x(3);
    % h = h_real/R_0;
    % V = V_real/v_c;

    % if h_real > 7
    %     gamma = 1.4;
    %     T = -23.4 - 2.22*h_real + 273.15;   % K
    % else
    %     gamma = 1.22;
    %     T = -31 - 0.998*h_real + 273.15;
    % end


    % New Mars-GRAM temperature profiles
    H_profiles = [0 1 45 60 75 82 105 130];
    T_profiles = [240 225 140 160 125 130 115 180];
    % if h_real <= 0
    %     T = -31 - 0.998*h_real + 273.15;
    if h_real <= 1
        T = interp1(H_profiles(1:2), T_profiles(1:2), h_real, 'linear', 'extrap');  
    elseif h_real > 1 && h_real <= 45
        T = interp1(H_profiles(2:3), T_profiles(2:3), h_real, 'linear');  
    elseif h_real > 45 && h_real <= 60
        T = interp1(H_profiles(3:4), T_profiles(3:4), h_real, 'linear');  
    elseif h_real > 60 && h_real <= 75
        T = interp1(H_profiles(4:5), T_profiles(4:5), h_real, 'linear');  
    elseif h_real > 75 && h_real <= 82
        T = interp1(H_profiles(5:6), T_profiles(5:6), h_real, 'linear');  
    elseif h_real > 82 && h_real <= 105 
        T = interp1(H_profiles(6:7), T_profiles(6:7), h_real, 'linear');  
    elseif h_real > 105 
        T = interp1(H_profiles(7:8), T_profiles(7:8), h_real, 'linear', 'extrap');  
    end

    if T < 235
        gamma = 1.4;
    %     T = -23.4 - 2.22*h_real + 273.15;   % K
    else
        gamma = 1.3;
    %     T = -31 - 0.998*h_real + 273.15;
    end

    % Mach number
    M = (V_real*1000)/sqrt(gamma*R_s*T);       % sqrt(gamma*R_s*T): speed of sound
    % if ~isreal(M)
    %     idx = 1;
    % else
    %     idx = 0;
    % end
        
    end



end