%Find number of integers x in Z^D s.t. ||x||_2 <= k(m) for integer radius k(m)
%Dynamic programming method. 

%Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%k          [1 x M] Radius, integer
%N:         Number of dimensions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%c:         [1 x M] Count

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:

%GCP_DP(0, 1)       %1
%GCP_DP(1, 1)       %3
%GCP_DP(1, 2)       %5
%GCP_DP(1, 3)       %7

%GCP_DP(3, 4)       %425
%GCP_DP(8, 4)       %20185

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare with  direct computation

% k = 0:8;
% N = 5;
%tic; c_direct = GCP_direct(k, N); toc
%tic; c_DP = GCP_DP(k, N); toc
%err = norm(c_direct - c_DP)


function c = GCP_DP(k, N)

M = numel(k);

k_ceil = ceil(max(k));

Q = k_ceil^2 + 1;

C = zeros(Q, N);

%Base case n = 1
for q = 0:(Q-1)
    C(q + 1, 1) = 1 + 2 * floor(sqrt(q));
end
%Higher dimension cases n > 1
for n = 2:N
    for q = 0:(Q-1)
        % for m = -floor(sqrt(q)) : floor(sqrt(q))
        %     C(q + 1, n) = C(q + 1, n) + C((q - m^2) + 1, n - 1);
        % end

        %Unrolled
        C(q + 1, n) = C(q + 1, n) + C(q + 1, n - 1); %m = 0
        for m = 1:floor(sqrt(q))
             C(q + 1, n) = C(q + 1, n) + 2 * C((q - m^2) + 1, n - 1);
        end
    end
end

%Output
c = zeros(size(k));
for m = 1:M
    c(m) = C((k(m)^2) + 1, N);
end