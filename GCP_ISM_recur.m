%Compute weighted image-source contributions in R^N within distance sqrt(q)
%via weghted Gauss Circle Problem centered on receiver coordinate r,
%source coordinate s, and orthotope with dimensions l, reflection coefficients
%gamma_pos for wall along +axis, gamma_neg for wall along -axis.

%Recursive method: Compute S(q, N) from sub-solutions of weighted S(q-u(m)^2, N-1).

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input

%q:             Scalar, squared distance (meters^2)

%s:             [1 x N] Source offset from origin (meters)
%r:             [1 x N] Receiver offset from origin (Meters)
%l:             [1 x N] Orthotope dimensions (meters)

%gamma_pos:     [1 x N] Reflection coefficient for wall on +axis
%gamma_neg:     [1 x N] Reflection coefficient for wall on -axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%sq:        Scalar, weighted image-source contributions within distance sqrt(q)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:

% gamma_pos = [0.5, 0.3, 0.9];
% gamma_neg = [0.4, 0.2, 0.8];
% Q = 100;
% s = [1, 2, 1];
% r = [2, 1, 1];
% l = [5, 6, 3];
% 
% tic; sq_recur = GCP_ISM_recur(Q, s, r, l, gamma_pos, gamma_neg); toc
% tic; sq_direct = GCP_ISM_direct(Q, s, r, l, gamma_pos, gamma_neg); toc

% err = norm(sq_recur - sq_direct)

function sq = GCP_ISM_recur(q, s, r, l, gamma_pos, gamma_neg)

N = numel(s);

if N == 1

    sq = GCP_ISM_1D_geom_count(q, s, r, l, gamma_pos, gamma_neg);

elseif N > 1
    
    sqrt_q = sqrt(q);
    [bq] = g_upper((r(N) + sqrt_q) / l(N), s(N) / l(N));
    [aq] = g_lower((r(N) - sqrt_q) / l(N), s(N) / l(N));

    sq = 0;
    nidx = 1:(N-1); %Lower dimension indices

    % aq = aq - 1;
    % bq = bq + 1;

    if bq >= aq
        for m = aq:bq
            um = u_vector(m, s(N), r(N), l(N));
            sq = sq + GCP_ISM_recur(q-um^2, s(nidx), r(nidx), l(nidx), gamma_pos(nidx), gamma_neg(nidx)) ...
                * gamma_pos(N)^abs(floor((m+1)/2)) * gamma_neg(N)^abs(floor((1-m)/2));
        end
    end
    
else
    sq = 0;
end
