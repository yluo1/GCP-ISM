%Generate RIRs along each sub-dimension
%Phase-flip reflections
%Note that impedance's phase should not invert unless medium is reflecting off of boundary with lower impedance
%e.g. air in tube (high impedance) phase-flip when reflecting off of tube's opening exposed to open-air

%Air reflecting off walls should not phase-flip

%Author: Yuancheng Luo, 2026

function plot_GCP_ISM_refl_phaseflip

[s_full, r_full, l_full, gamma_pos_full, gamma_neg_full, T] = get_default_RIR_params;

%T = 1;
climits = [-160, -40];


RT60_dB_hi = -20;
RT60_dB_lo = -40;

n = 6;

ndims = 1:n

s = s_full(ndims);
r = r_full(ndims);
l = l_full(ndims);

%No phase-flip
gamma_pos = gamma_pos_full(ndims);
gamma_neg = gamma_neg_full(ndims);


h = RIR_GCP_ISM_LUT(T, s, r, l, gamma_pos, gamma_neg, 'lambda', [1], 'direction', 'inverse', ...
    'enable_disp', true, 'clim', climits, 'fig_size', [900 900] * (3/4), ...
    'font_size', 22, 'legend_location', 'east', 'RT60_dB_hi', RT60_dB_hi, 'RT60_dB_lo', RT60_dB_lo);


%Phase-flip both +- facing walls
gamma_pos = -gamma_pos_full(ndims);
gamma_neg = -gamma_neg_full(ndims);


h = RIR_GCP_ISM_LUT(T, s, r, l, gamma_pos, gamma_neg, 'lambda', [1], 'direction', 'inverse', ...
    'enable_disp', true, 'clim', climits, 'fig_size', [900 900] * (3/4), ...
    'font_size', 22, 'legend_location', 'southeast', 'RT60_dB_hi', RT60_dB_hi, 'RT60_dB_lo', RT60_dB_lo);

exportgraphics(gcf, ['inverse_GCP_N_', num2str(n), '_phaseflip_posneg.png']);


% %Phase-flip only - facing walls
% gamma_pos = gamma_pos_full(ndims);
% gamma_neg = -gamma_neg_full(ndims);
% 
% 
% h = RIR_GCP_ISM_LUT(T, s, r, l, gamma_pos, gamma_neg, 'lambda', [1], 'direction', 'inverse', ...
%     'enable_disp', true, 'clim', clim, 'fig_size', [900 900] * (3/4), ...
%     'font_size', 22, 'legend_location', 'southeast');
% 
% exportgraphics(gcf, ['inverse_GCP_N_', num2str(n), '_phaseflip_neg.png']);


%Phase-flip only + facing walls
gamma_pos = -gamma_pos_full(ndims);
gamma_neg = gamma_neg_full(ndims);


h = RIR_GCP_ISM_LUT(T, s, r, l, gamma_pos, gamma_neg, 'lambda', [1], 'direction', 'inverse', ...
    'enable_disp', true, 'clim', climits, 'fig_size', [900 900] * (3/4), ...
    'font_size', 22, 'legend_location', 'southeast', 'RT60_dB_hi', RT60_dB_hi, 'RT60_dB_lo', RT60_dB_lo);

exportgraphics(gcf, ['inverse_GCP_N_', num2str(n), '_phaseflip_pos.png']);


%Phase-flip alternating walls
gamma_pos = gamma_pos_full(ndims) .* (-1).^((ndims) + 0);
gamma_neg = gamma_neg_full(ndims) .* (-1).^((ndims) + 1);

% gamma_pos = gamma_pos_full(ndims) .* (-1).^((ndims) + 1);
% gamma_neg = gamma_neg_full(ndims) .* (-1).^((ndims) + 0);


h = RIR_GCP_ISM_LUT(T, s, r, l, gamma_pos, gamma_neg, 'lambda', [1], 'direction', 'inverse', ...
    'enable_disp', true, 'clim', climits, 'fig_size', [900 900] * (3/4), ...
    'font_size', 22, 'legend_location', 'southeast', 'RT60_dB_hi', RT60_dB_hi, 'RT60_dB_lo', RT60_dB_lo);

exportgraphics(gcf, ['inverse_GCP_N_', num2str(n), '_phaseflip_alt.png']);



