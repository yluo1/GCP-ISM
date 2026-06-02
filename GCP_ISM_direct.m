%Compute weighted image-source contributions in R^N within distance sqrt(q)
%via weghted Gauss Circle Problem centered on receiver coordinate r,
%source coordinate s, and orthotope with dimensions l, reflection coefficients
%gamma_pos for wall along +axis, gamma_neg for wall along -axis.

%Direct method: Finds all image-source coordinates within sqrt(q) Eucildean norm
%take the weighted sum

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
% Q = 1000;
% s = [1, 2, 1];
% r = [2, 1, 1];
% l = [5, 6, 3];
% 
% [sq] = GCP_ISM_direct(Q, s, r, l, gamma_pos, gamma_neg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Lamba scaling equivalences

% gamma_pos = [0.5, 0.3, 0.9];
% gamma_neg = [0.4, 0.2, 0.8];
% Q = 10*pi;
% s = [1, 2, 1];
% r = [2, 1, 1];
% l = [5, 6, 3];
% 
% lambda_list = [1, 2, 3];
% N_lambda = numel(lambda_list);
% sq_lambda = zeros(N_lambda, 1);
% for n = 1:N_lambda
%   lambda = lambda_list(n);
%   [sq_lambda(n)] = GCP_ISM_direct(lambda^2 * Q, lambda * s, lambda * r, lambda * l, gamma_pos, gamma_neg);
% end
% sq_lambda

function [sq] = GCP_ISM_direct(q, s, r, l, gamma_pos, gamma_neg)

N = numel(s);

%Compute bounds
ub = zeros(1, N);
lb = zeros(1, N);
x_list = cell(1, N);
for n = 1:N
    while u_vector(ub(n), s(n), r(n), s(n))^2 < q
        ub(n) = ub(n) + 1;
    end
    while u_vector(lb(n), s(n), r(n), s(n))^2 < q
        lb(n) = lb(n) - 1;
    end
    x_list{n} = lb(n):ub(n);
end

%Compute lattice coordinates
V_set = table2array(combinations(x_list{:})); %[M x N]

U_set = u_vector(V_set, s, r, l);   %[M x N]
R_sq = sum(U_set.^2, 2);            %[M x 1]

mask_valid = (R_sq <= q); %[M x 1]

R_set = double(mask_valid);
for n = 1:N
    R_set(mask_valid) = R_set(mask_valid) .* ... 
        gamma_pos(n).^abs(floor((V_set(mask_valid, n) + 1)/2) )  .* ...
        gamma_neg(n).^abs(floor((1 - V_set(mask_valid, n))/2) );
end
sq = sum(R_set);

%V = -10:10
%gamma_pos_pow = abs(floor((V + 1)/2) )
%gamma_neg_pow = abs(floor((1 - V)/2) )