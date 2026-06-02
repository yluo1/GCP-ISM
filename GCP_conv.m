%Find number of integers x in Z^D s.t. ||x||_2 <= k(m) for integer radius k(m)
%Convolution method, frequency domain

%Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%k          [1 x M] Radius, integer
%N:         Number of dimensions
%mode:      String,     'padded_FFT'     %O(N k^2) memory, O(N k^2 log(k)) compute
%                       'trunc_FFT'      %O(2 k^2) memory, O(N k^2 log(k)) compute
%                       'conv'           %O(N k^2) memory, O(N k^3)        compute

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%c:         [1 x M] Count

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:

%GCP_conv(0, 1)       %1
%GCP_conv(1, 1)       %3
%GCP_conv(1, 2)       %5
%GCP_conv(1, 3)       %7

%GCP_conv(3, 4)       %425
%GCP_conv(8, 4)       %20185

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare with direct computation

% k = 0:12;
% N = 5;
%tic; c_direct = GCP_direct(k, N); toc
%tic; c_conv = GCP_conv(k, N); toc
%err = norm(c_direct - c_conv)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Compare with DP

%mode = 'padded_FFT';
%%mode = 'trunc_FFT';
%%mode = 'conv';

% k = 0:100;
% N = 6;
%tic; c_DP = GCP_DP(k, N); toc
%tic; c_conv = GCP_conv(k, N, mode); toc
%err = norm(c_DP - c_conv)


function c = GCP_conv(k, N, mode)

if nargin < 3 || isempty(mode)
    mode = 'padded_FFT';
end

M = numel(k);

k_ceil = ceil(max(k));

Q = k_ceil^2 + 1;

C1 = zeros(Q, 1);

%Base case n = 1
for q = 0:(Q-1)
    C1(q + 1) = 1 + 2 * floor(sqrt(q));
end

%Perfect square kernel
f = zeros(Q, 1);
f(1) = 1;
f( ((1:k_ceil).^2) + 1) = 2;

if strcmp(mode, 'padded_FFT')

    %Convolution in frequency-domain after collapsing recurrence relation
    %O(N k^2) memory, O(N k^2 log(k)) compute
    C = ifft((fft(f, N * Q).^(N - 1)) .* fft(C1, N * Q));
    C = real(round(C));

elseif strcmp(mode, 'trunc_FFT')

    %Convolution in frequency-domain per dimension
    %O(2 k^2) memory, O(N k^2 log(k)) compute

    C = C1;
    f_spec = fft(f, 2 * Q);
    for n = 2:N
        C = ifft(f_spec .* fft(C, 2 * Q));
        C = round(C(1:Q));
    end

elseif strcmp(mode, 'conv')

    %Convolution in time-domain per dimension
    %O(N k^2) memory, O(N k^3) compute

    C = C1;
    for n = 2:N
        C = conv(f, C);
    end

else
    error('Unsupported mode');
end
%Output
c = zeros(size(k));
for m = 1:M
    c(m) = C((k(m)^2) + 1);
end