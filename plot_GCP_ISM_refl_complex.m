%Generate RIRs along each sub-dimension
%Phase-rotate reflection coefficients

%Author: Yuancheng Luo, 2026

function plot_GCP_ISM_refl_complex

[s_full, r_full, l_full, gamma_pos_full, gamma_neg_full, T] = get_default_RIR_params;

%Number of frequency points from DC to Nyquist 
%P = 8;
%P = 9;
P = 17; %OG
%P = 33;
%P = 64;

if P == 1
    w = 0;
else 
    w = linspace(0, pi, P)';
end

%Design material reflection filter
h_refl = filter_two_tap_FIR(0, -2.5, false)
H_refl = freqz(h_refl, 1, [w; pi]); H_refl = H_refl(1:end-1);  %[P x 1]

%T = 1;

n = 3;
%n = 6;

ndims = 1:n

s = s_full(ndims);
r = r_full(ndims);
l = l_full(ndims);

% %Cascade original reflection cofficients with filter
gamma_pos = H_refl * gamma_pos_full(ndims);
gamma_neg = H_refl * gamma_neg_full(ndims);

%Cascade original reflection cofficients with filter
% gamma_pos = H_refl * gamma_pos_mag(ndims) .* (-1).^((ndims) + 0);;
% gamma_neg = H_refl * gamma_neg_mag(ndims) .* (-1).^((ndims) + 1); ;


%Generate RIR
h = RIR_GCP_ISM_LUT_freq(T, s, r, l, gamma_pos, gamma_neg, 'lambda', [1], 'mode', 'ifft', ...
    'enable_disp', true, 'clim', [-120, -40], 'legend_location', 'east');

exportgraphics(gcf, ['inverse_GCP_N_', num2str(n), '_complex.png']);
