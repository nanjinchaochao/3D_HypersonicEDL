%{
A two-tier particle selection approach in order to generate the next 
"optimal" particle by improving its sampling efficiency: a large number of 
candidate is first drawn from a uniform distribution using a pseudo-random 
generator at the initial time. The first level of selection employs a 
projective distance threshold function such that candidates that lie 
?too close? to the existing ensemble in terms of projective distance are 
ruled out. THen, the remaining, much smaller set of candidates is ranked 
in terms of discrepancy.
%}

function P_new = ParticleGeneration_SamplingEfficiency(sam)
global alpha;
global N; 
global Num_candidate;

NSAM = size(sam, 1);                                    % # number of the current samples
cand = rand([Num_candidate, N]);                      % Generate candidates from a uniform distribution in the "low dimensional" case
d_min = alpha/(NSAM + 1);                               % Minimum allowed projected distance
ntr = 1;                                                % Count the number of samples satisfied the non-collapsing threshold function
Discrepancy = [];              
cand_remain = [];

for ctr = 1 : Num_candidate
   %% Non-collapsing property: Minimun projected distance
   dis = abs(sam - repmat(cand(ctr, :), NSAM, 1));
   proj = min(dis, [], 2);
   proj_d = min(proj);
   if proj_d >= d_min     
       %% Space-filling property: centered L2-Discrepancy 
       % (how to generalize this such that the user has flexibility to select different criterion?)
       P_r = [sam; cand(ctr, :)];                       
       C_r = repmat(cand(ctr, :), NSAM + 1, 1);
       % First Part
       P1_r = 1 + abs(cand(ctr, :) - 0.5)/2 - (cand(ctr, :) - 0.5).^2/2;
       D1_r = (2/(NSAM + 1))*prod(P1_r, 2);
       % Second Part
       P2_r = 1 + 0.5*abs(P_r - 0.5) + 0.5*abs(C_r - 0.5) - 0.5*abs(P_r - C_r);
       P2_r = prod(P2_r, 2);
       P2_r(1 : end - 1) = 2*P2_r(1 : end - 1);
       D2_r = (1/(NSAM + 1)^2)*sum(P2_r);
       % Discrepancy
       Discrepancy(ntr, 1) = (13/12)^N - D1_r + D2_r;
       cand_remain(ntr, :) = cand(ctr, :);
       ntr = ntr + 1;

   end
end
    if isempty(Discrepancy) == 1
        fprintf('Did not find the next optimal simulation point.\n');
    else
        % Find the "optiaml" next particle with minimum discrepancy
        % (how to select particles in batches?)
        [~,I] = min(Discrepancy);                          
        P_new = cand_remain(I, :);     
    end
end 