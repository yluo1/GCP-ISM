%Compute weighted image-source contributions in R^N within distance sqrt(q)
%via weghted Gauss Circle Problem centered on receiver coordinate r,
%source coordinate s, and orthotope with dimensions l, reflection coefficients
%gamma_pos for wall along +axis, gamma_neg for wall along -axis.

%Geometric series solution for single entrant in LUT
%Jitters image-source coordinate by random offset

%Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input

%q:             Scalar, maximum squared distance (meters^2)
%s:             Scalar, source offset from origin (meters)
%r:             Scalar, receiver offset from origin (meters)
%l:             Scalar, orthotope dimensions (meters)

%gamma_pos:     Scalar, reflection coefficient for wall on +axis
%gamma_neg:     Scalar, reflection coefficient for wall on -axis

%rand_off_LUT:  [M x 1] random image-source coordinate offset LUT per lattice coordinate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%sq:            Scalar, image-source volume for squared distances (meters) <= q

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage

% gamma_pos = 0.5;
% gamma_neg = 0.31;
% Q = 10000;
% s = 0.49;
% r = 0.1;
% l = 1;
% 
% sq_direct      = GCP_ISM_direct(Q, s, r, l, gamma_pos, gamma_neg);
% sq_geometric   = GCP_ISM_1D_geom_count_asym(Q, Q, s, r, l, gamma_pos, gamma_neg);
% 
% err = norm(sq_direct - sq_geometric)

function sq = GCP_ISM_1D_geom_count_rand(q, s, r, l, gamma_pos, gamma_neg, rand_off_LUT)

M = numel(rand_off_LUT);

sqrt_q = sqrt(q);
[bq] = g_upper((r + sqrt_q) / l, s / l);
[aq] = g_lower((r - sqrt_q) / l, s / l);

[bq] = g_upper((r + sqrt_q + rand_off_LUT(mod(bq, M) + 1) ) / l, s / l);
[aq] = g_lower((r - sqrt_q + rand_off_LUT(mod(aq, M) + 1) ) / l, s / l);


% %Check
% u_aq = u_vector(aq, s, r, l);
% u_bq = u_vector(bq, s, r, l);
% all(diff([-sqrt_q, u_aq, u_bq, sqrt_q]) >= 0) %Is monotonic increasing

if bq >= aq

   gamma_prod = gamma_pos * gamma_neg;

   if gamma_prod == 1
       sq = bq - aq + 1;
       
   else
       sq = ((1 - gamma_prod^(ceil(abs((bq + 1)/2)) ) )   + ...            %Even, nonneg
             (1 - gamma_prod^(ceil(abs(bq/2)) )) * (gamma_pos)  + ...    %Odd, nonneg
             (1 - gamma_prod^(ceil(abs((aq - 1)/2)) ))   + ...            %Even, neg
             (1 - gamma_prod^(ceil(abs(aq/2)) )) * (gamma_neg)    ...    %Odd, neg
             ) /  ...
             (1 - gamma_prod) ...
             - 1;
   end
   
else 
    sq = 0;
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function s = geom_count(gamma_pos, gamma_neg, up, lo)
% 
% gamma_prod = gamma_pos * gamma_neg;
% 
% %Pass
% % s = ((1 - gamma_prod^(floor(abs(up/2)) + 1 ))   + ...            %Even, nonneg
% %      (1 - gamma_prod^(floor(abs(up/2)) + rem(abs(up), 2) )) * (gamma_pos)  + ...    %Odd, nonneg
% %      (1 - gamma_prod^(floor(abs(lo/2)) + 1 ))  + ...            %Even, neg
% %      (1 - gamma_prod^(floor(abs(lo/2)) + rem(abs(lo), 2) )) * (gamma_neg)    ...    %Odd, neg
% %      ) /  ...
% %      (1 - gamma_prod) ...
% %      - 1;
% 
% % %Pass
% % s = ((1 - gamma_prod^(floor(abs(up/2)) + 1 ) )   + ...            %Even, nonneg
% %      (1 - gamma_prod^(ceil(abs(up/2)) )) * (gamma_pos)  + ...    %Odd, nonneg
% %      (1 - gamma_prod^(floor(abs(lo/2)) + 1 ))   + ...            %Even, neg
% %      (1 - gamma_prod^(ceil(abs(lo/2)) )) * (gamma_neg)    ...    %Odd, neg
% %      ) /  ...
% %      (1 - gamma_prod) ...
% %      - 1;
% 
% %Pass
% s = ((1 - gamma_prod^(ceil(abs((up + 1)/2)) ) )   + ...            %Even, nonneg
%      (1 - gamma_prod^(ceil(abs(up/2)) )) * (gamma_pos)  + ...    %Odd, nonneg
%      (1 - gamma_prod^(ceil(abs((lo - 1)/2)) ))   + ...            %Even, neg
%      (1 - gamma_prod^(ceil(abs(lo/2)) )) * (gamma_neg)    ...    %Odd, neg
%      ) /  ...
%      (1 - gamma_prod) ...
%      - 1;

