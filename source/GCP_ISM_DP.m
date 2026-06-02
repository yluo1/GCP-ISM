%Compute weighted image-source contributions in R^N within distance sqrt(q)
%via weghted Gauss Circle Problem centered on receiver coordinate r,
%source coordinate s, and orthotope with dimensions l, reflection coefficients
%gamma_pos for wall along +axis, gamma_neg for wall along -axis.

%Dynamic programming method: Compute S(q, N) from sub-solutions of weighted S(m, N-1).
%Requires integer Q, s, r, l

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input

%Q:             Scalar, maximum squared distance (meters^2)

%s:             [1 x N] Source offset from origin (meters)
%r:             [1 x N] Receiver offset from origin (Meters)
%l:             [1 x N] Orthotope dimensions (meters)

%gamma_pos:     [1 x N] Reflection coefficient for wall on +axis
%gamma_neg:     [1 x N] Reflection coefficient for wall on -axis

%jitter_coord_bnd:   [1 x 2]    Jitter the image source coordinates by 
%                               unifrnd(min(jitter_coord_bnd), max(jitter_coord_bnd))
%                               (Default = [0, 0] is disabled)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%S:     [(Q+1) x N] LUT of image-source volumes S(q+1, N) for squared distances (meters) <= q

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:

% gamma_pos = [0.5, 0.3, 0.9];
% gamma_neg = [0.4, 0.2, 0.8];
% Q = 100;
% s = [1, 2, 1];
% r = [2, 1, 1];
% l = [5, 6, 3];
% 
% S_DP = GCP_ISM_DP(Q, s, r, l, gamma_pos, gamma_neg);

% N = numel(s);
% S_direct = zeros(Q+1, N);
% for n = 1:N
%     for q = 0:Q
%         S_direct(q+1, n) = GCP_ISM_direct(q, s(1:n), r(1:n), l(1:n), gamma_pos(1:n), gamma_neg(1:n));
%     end
% end
% err = norm(S_DP - S_direct)

function S = GCP_ISM_DP(Q, s, r, l, gamma_pos, gamma_neg, jitter_coord_bnd)

if nargin < 7 || isempty(jitter_coord_bnd)
    jitter_coord_bnd = [0, 0];
end

N = numel(s);
S = zeros([(Q+1), N]);

if isequal(jitter_coord_bnd, [0, 0]) %No jitter

    %Compute LUT basecase
    S(:, 1) = GCP_ISM_basecase(Q, s(1), r(1), l(1), gamma_pos(1), gamma_neg(1), 'geometric');
    
    %Iterate over higher dimensions
    for n = 2:N    
        for q = 0:Q

            %Compute lattice coordinate bounds
            sqrt_q = sqrt(q);
            [bq] = g_upper((r(n) + sqrt_q) / l(n), s(n) / l(n));
            [aq] = g_lower((r(n) - sqrt_q) / l(n), s(n) / l(n));
    
            if bq >= aq
                for m = aq:bq
                    um = u_vector(m, s(n), r(n), l(n));

                    S(q+1, n) = S(q+1, n)    + ...
                                S(q-um^2 + 1, n-1)  * gamma_pos(n)^abs(floor((m+1)/2)) * gamma_neg(n)^abs(floor((1-m)/2));   
                end
            end
        end
    end

else %With jitter

    %Compute LUT basecase
    S(:, 1) = GCP_ISM_basecase(Q, s(1), r(1), l(1), gamma_pos(1), gamma_neg(1), 'geometric', jitter_coord_bnd);
    
    %Iterate over higher dimensions
    for n = 2:N

        sqrt_Q = sqrt(Q);
        [bq] = g_upper((r(n) + sqrt_Q) / l(n), s(n) / l(n));
        [aq] = g_lower((r(n) - sqrt_Q) / l(n), s(n) / l(n));
        
        rand_off_LUT = unifrnd(min(jitter_coord_bnd), max(jitter_coord_bnd), abs(bq-aq) + 1, 1);
        rand_off_LUT(1) = 0; %Exclude orignal source

        M = numel(rand_off_LUT);

        for q = 0:Q

            %Compute lattice coordinate bounds
            sqrt_q = sqrt(q);
            [bq] = g_upper((r(n) + sqrt_q) / l(n), s(n) / l(n));
            [aq] = g_lower((r(n) - sqrt_q) / l(n), s(n) / l(n));

            [bq] = g_upper((r(n) + sqrt_q + rand_off_LUT(mod(bq, M) + 1) ) / l(n), s(n) / l(n));
            [aq] = g_lower((r(n) - sqrt_q + rand_off_LUT(mod(aq, M) + 1) ) / l(n), s(n) / l(n));
        
            if bq >= aq
                for m = aq:bq
                    um = u_vector(m, s(n), r(n), l(n)) + rand_off_LUT(mod(m, M) + 1);
                    
                    S(q+1, n) = S(q+1, n)    + ...
                                S(min(max(1, q-round(um^2) + 1), Q+1), n-1)  * gamma_pos(n)^abs(floor((m+1)/2)) * gamma_neg(n)^abs(floor((1-m)/2));

                end
            end
        end
    end

end

