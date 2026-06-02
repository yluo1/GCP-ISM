%Find number of integers x in Z^D s.t. ||x||_2 <= k(m)
%Direct method generates all x in ||x||_1 <= k(m) and checks above condition

%Optionally unrolls loops for N <= 6.
%Otherwise compute and memory cost ~ (1+2*k(m))^N

%Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%k                  [1 x M] Radius
%N:                 Number of dimensions
%unroll:            Logical, if true, use the unrolled solutions for N<=6

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%c:         [1 x M] Count

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:

%GCP_direct(0, 1)       %1
%GCP_direct(1, 1)       %3
%GCP_direct(1, 2)       %5
%GCP_direct(1, 3)       %7

%GCP_direct(3, 4)       %425
%GCP_direct(8, 4)       %20185

%GCP_direct(sqrt(2), 2)       %9
%GCP_direct(sqrt(2)-eps, 2)   %5

function c = GCP_direct(k, N, unroll)

if nargin < 3 || isempty(unroll)
    unroll = false;
end

M = numel(k);

k_ceil = ceil(max(k));
c = zeros(size(k));

if N <= 6 && unroll

    if  N == 1  %Unrolled
    
        v = zeros(1, N);
        idx = -k_ceil:k_ceil;
        for i = idx
            v(1) = i;
    
            R_sq = sum(v.^2, 2);
            for m = 1:M
                c(m) = c(m) + (R_sq <= k(m)^2);
            end
        end
        
    elseif N == 2 %Unrolled
    
        v = zeros(1, N);
        idx = -k_ceil:k_ceil;
        for i = idx
            v(1) = i;
            for j = idx
                v(2) = j;    
    
                R_sq = sum(v.^2, 2);
                for m = 1:M
                    c(m) = c(m) + (R_sq <= k(m)^2);
                end
            end
        end
    
    elseif N == 3 %Unrolled
    
        v = zeros(1, N);
        idx = -k_ceil:k_ceil;
        for i = idx
            v(1) = i;
            for j = idx
                v(2) = j;      
                for kk = idx
                    v(3) = kk;
    
                    R_sq = sum(v.^2, 2);
                    for m = 1:M
                        c(m) = c(m) + (R_sq <= k(m)^2);
                    end
                end
            end
        end
    
    elseif N == 4 %Unrolled
    
        v = zeros(1, N);
        idx = -k_ceil:k_ceil;
        for i = idx
            v(1) = i;
            for j = idx
                v(2) = j;      
                for kk = idx
                    v(3) = kk;
                    for ll = idx
                        v(4) = ll;
    
                        R_sq = sum(v.^2, 2);
                        for m = 1:M
                            c(m) = c(m) + (R_sq <= k(m)^2);
                        end
                    end
                end
            end
        end
    
    elseif N == 5 %Unrolled
    
        v = zeros(1, N);
        idx = -k_ceil:k_ceil;
        for i = idx
            v(1) = i;
            for j = idx
                v(2) = j;      
                for kk = idx
                    v(3) = kk;
                    for ll = idx
                        v(4) = ll;
                        for mm = idx
                            v(5) = mm;
    
                            R_sq = sum(v.^2, 2);
                            for m = 1:M
                                c(m) = c(m) + (R_sq <= k(m)^2);
                            end
                        end
                    end
                end
            end
        end
    
    elseif N == 6 %Unrolled
    
        v = zeros(1, N);
        idx = -k_ceil:k_ceil;
        for i = idx
            v(1) = i;
            for j = idx
                v(2) = j;      
                for kk = idx
                    v(3) = kk;
                    for ll = idx
                        v(4) = ll;
                        for mm = idx
                            v(5) = mm;
                            for nn = idx
                                v(6) = nn;
        
                                R_sq = sum(v.^2, 2);
                                for m = 1:M
                                    c(m) = c(m) + (R_sq <= k(m)^2);
                                end
                             end
                        end
                    end
                end
            end
        end
    end

else %General case

    x_list = repmat({-k_ceil:k_ceil}, 1, N);
    V_comb = combinations(x_list{:});
    V_set = table2array(V_comb);
    
    R_sq = sum(V_set.^2, 2);
    
    for m = 1:M
        c(m) = sum(R_sq <= k(m)^2);
    end

end

