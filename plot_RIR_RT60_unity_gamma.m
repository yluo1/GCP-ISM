%Plot GCP_ISM RIR with various gamma under various phase-inversion

%Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_RIR_RT60_unity_gamma 

%clear all

[s_full, r_full, l_full, gamma_pos_full, gamma_neg_full, T] = get_default_RIR_params;

N = 3;

ndims = 1:N

s = s_full(ndims);
r = r_full(ndims);
l = l_full(ndims);

T = 1.0;

%T60 = T * ones(1, N);
%Alternating phase
gamma_pos_mask = (-1).^(1:N);
gamma_neg_mask = (-1).^((1:N) + 1);

%Negatve gamma
% gamma_pos_mask = -ones(1, N);
% gamma_neg_mask = -ones(1, N);

%One sided walls phase invert
% gamma_pos_mask = ones(1, N);
% gamma_neg_mask = -ones(1, N);

%In-phase
% gamma_pos_mask = (1).^(1:N);
% gamma_neg_mask = (1).^((1:N) + 1);

%xi = 4;
%gamma = compute_RT60_gain_from_dim(l(1:N), T60, xi);
gamma = ones(1, N);

gamma_pos = gamma .* gamma_pos_mask;
gamma_neg = gamma .* gamma_neg_mask;


h = RIR_GCP_ISM_LUT(T, s(1:N), r(1:N), l(1:N), gamma_pos, gamma_neg, 'lambda', [1], 'direction', 'inverse',...
                'enable_disp', true, 'disp_EDC_fig', true);



Fs = 48000;
c = 343; %Speed of sound
M_h = ceil(T * Fs);  %Number of taps

Ts = 1 / Fs;
c_Ts = c * Ts;


k_max = M_h / Fs * c;
q_max = k_max^2;

lambda = 1;
Q = ceil(q_max);
Q_lambda = Q * lambda^2;

S = GCP_ISM_conv(Q_lambda, lambda * s, lambda * r, lambda * l, gamma_pos, gamma_neg);  
%S = GCP_ISM_DP(Q_lambda, lambda * s, lambda * r, lambda * l, gamma_pos, gamma_neg);  

figure; plot(S(:, N)); grid on; axis tight;