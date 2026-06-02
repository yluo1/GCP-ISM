%Generate and plot Forward GCP-ISM

%Author: Yuancheng Luo, 2026

function plot_GCP_ISM_inverse

[s_full, r_full, l_full, gamma_pos_full, gamma_neg_full, T] = get_default_RIR_params;

N = numel(s_full);

ndims = 1:3

s = s_full(ndims);
r = r_full(ndims);
l = l_full(ndims);
gamma_pos = gamma_pos_full(ndims);
gamma_neg = gamma_neg_full(ndims);

h = RIR_GCP_ISM_LUT(T, s, r, l, gamma_pos, gamma_neg, 'lambda', [1], 'direction', 'inverse',...
    'enable_disp', true);

exportgraphics(gcf, 'inverse_GCP.png');