%Generate RIRs along each sub-dimension
%Author: Yuancheng Luo, 2026

function plot_GCP_ISM_unrolled

[s_full, r_full, l_full, gamma_pos_full, gamma_neg_full, T] = get_default_RIR_params;

N = numel(s_full);

%for n = 1:N
%for n = 1
%for n = 2
%for n = 3
%for n = 4
%for n = 5
%for n = 6

for n = 4:6

    ndims = 1:n

    s = s_full(ndims);
    r = r_full(ndims);
    l = l_full(ndims);
    gamma_pos = gamma_pos_full(ndims);
    gamma_neg = gamma_neg_full(ndims);
    
    h = RIR_GCP_ISM_LUT(T, s, r, l, gamma_pos, gamma_neg, 'lambda', [1], 'direction', 'inverse', ...
        'enable_disp', true, 'clim', [-160, -40], 'fig_size', [900 900] * (3/4), ...
        'font_size', 22);

    exportgraphics(gcf, ['inverse_GCP_N_', num2str(n), '.png']);

end

;
