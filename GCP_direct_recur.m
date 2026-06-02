%Find number of integers x in Z^D s.t. ||x||_2 <= k
%Use the cardinality recurrence relation without memoization
%Requires error correction due to rounding errors in numerical precision

%Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%k          [1 x M] Radius
%N:         Number of dimensions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%c:         [1 x M] Count

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:

%GCP_direct_recur(0, 1)       %1
%GCP_direct_recur(1, 1)       %3
%GCP_direct_recur(1, 2)       %5
%GCP_direct_recur(1, 3)       %7

%GCP_direct_recur(3, 4)       %425
%GCP_direct_recur(8, 4)       %20185

%GCP_direct_recur(14, 3)       %11513
%GCP_direct_recur(30, 3)       %113081

%GCP_direct_recur(11, 5)      %853399

function c = GCP_direct_recur(k, N)

M = numel(k);
c = zeros(size(k));

for m = 1:M
    c(m) = GCP_direct_recur_func(k(m), N);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = GCP_direct_recur_func(k, N)

k_floor = floor(k);

tol = 1e-10;
if k - k_floor > 1-tol %Numerical correction
    k_floor = ceil(k);
end

if N == 1
    c = 1 + 2 * k_floor;
else
    c = GCP_direct_recur(k, N - 1);
    for m = 1:k_floor
        c = c + 2 * GCP_direct_recur_func(sqrt(k^2 - m^2), N - 1);
    end
end