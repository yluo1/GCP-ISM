%Estimate echo density of RIR over sliding window
%Echo ~ number of samples outside sliding window's std

%Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%h:             [N x 1] filter
%Fs:            Sampling rate
%win_size:      Window size (odd)
%enable_disp:   Logical. If true, plot spectrogram and density curves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%dst:           [N x 1] Density profile (0 = sparse, > 1 dense) over time at t
%krt:           [N x 1] Kurtosis of window at t
%stdw:          [N x 1] standard deviation of window at t
%t:             [N x 1] time per second in seconds

function [dst, krt, stdw, t] = estimate_RIR_density(h, Fs, win_size)

if nargin < 3 || isempty(win_size)
%       win_size = 64;
        win_size = 128;
%       win_size = 256;
%       win_size = 512;
        win_size = win_size + 1;
end

if rem(win_size, 2) == 0 %Force odd
    win_size = win_size + 1;
end

h = h(:);
N = numel(h);      

t = (0:N-1) / Fs;
dst = zeros(N, 1);   %Density
krt = zeros(N, 1);   %Excess Kurtosis profile
stdw = zeros(N, 1);    %Standard deviation profile

%Sliding window centered on sample n
win_size_half = (win_size - 1) / 2;
for n = 1:N    
    n_start = max(n - win_size_half, 1);
    n_end   = min(n + win_size_half, N);
    x = h(n_start:n_end);
    
%    stdw(n) = std(x);
%    dst(n) = sum((abs(x-mean(x)) > stdw(n)) / (n_end - n_start + 1) ); %Fraction of samples outside 1 std.
 
    stdw(n) = sqrt(sum(x.^2) / (n_end - n_start + 1) );
    dst(n) = sum((abs(x) > stdw(n)) / (n_end - n_start + 1) ); %Fraction of samples outside 1 std.

    krt(n) = kurtosis(x);
end
dst = dst / erfc(1/sqrt(2));  %Denominator is expected number samples lying outside 1 std
