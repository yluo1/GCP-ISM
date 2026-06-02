%Test GCP RT60 controlled RIR

%Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_RIR_RT60_controlled

%clear all

[s_full, r_full, l_full, gamma_pos_full, gamma_neg_full, T] = get_default_RIR_params;

N = numel(s_full);



ndims = 1:6
%ndims = 1:3
%ndims = 1:2

s = s_full(ndims);
r = r_full(ndims);
l = l_full(ndims);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_dim_list = [1:6];
%max_dim_list = [3:6];

N_max_dims = numel(max_dim_list);

%Vary T
T_list = [0.5, 1, 2, 3, 4];
N_T = numel(T_list);

time_sec = zeros(N_T, N_max_dims);

h_inv_list = cell(N_T, N_max_dims);

%max_iter = 10;
max_iter = 1;

%Iterate over max dims
for n = 1:N_max_dims %Iterate over N
    N = max_dim_list(n);

    for i = 1:N_T %Iterate over T

        %Increase time and reflection coefficients
        T = T_list(i);

        T60 = T * ones(1, N);
        gamma_pos_mask = (-1).^(1:N);
        gamma_neg_mask = (-1).^((1:N) + 1);
           
        xi = 4;
        gamma = compute_RT60_gain_from_dim(l(1:N), T60, xi);
        gamma_pos = gamma .* gamma_pos_mask;
        gamma_neg = gamma .* gamma_neg_mask;

        
        for j = 1:max_iter
            tic;
            h_inv_list{i, n} = RIR_GCP_ISM_LUT(T, s(1:N), r(1:N), l(1:N), gamma_pos, gamma_neg, 'lambda', [1], 'direction', 'inverse',...
                'enable_disp', false);
            time_sec(i, n) = time_sec(i, n) + toc;
        end
        time_sec(i, n) = time_sec(i, n) / max_iter
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot RIR

options = [];
options.enable_disp = false;
options.clim = [-120, -40];
options.fig_size = [900 600] * (3/4);
options.font_size = 16;
options.legend_location = 'southeast';
options.Fs = 48000;

idx_t = 5;
idx_n = 6;
h_fig_sample = plot_ISM_RIR(h_inv_list{idx_t, idx_n}, 'Fs', options.Fs, 'name', ['Inverse GCP-ISM: N = ', num2str(max_dim_list(idx_n))], 'clim', options.clim, ...
    'fig_size', options.fig_size, 'font_size', options.font_size, 'legend_location', options.legend_location);

exportgraphics(h_fig_sample, 'RIR_RT60_controlled_sample.png');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fontsize = 18;
markersize = 12;
linewidth = 1.25;
marker_list = {'o', '*', 's', 'd', 'x', '+'};

h_fig = figure;
h_fig.Position = [100, 100, [600, 300] * 0.75]
legend_cell = {};
for i = 1:N_max_dims
    plot(T_list, time_sec(:, i), 'linewidth', linewidth, 'Marker', marker_list{i}, 'MarkerSize', markersize); hold on;
    legend_cell{end+1} = ['$N = ', num2str(max_dim_list(i)), '$'];
end
grid on; axis tight;
xlabel('RIR Length T (Seconds)', 'fontsize', fontsize-2, 'interpreter', 'latex');
ylabel('Runtime (Seconds)', 'fontsize', fontsize-2, 'interpreter', 'latex');
title('RIR Generation Time for Increasing Lengths $T$', 'fontsize', fontsize-2, 'interpreter', 'latex');
set(gca, 'fontsize', fontsize - 3);
xticks(T_list);
yticks([0, 2, 4, 6, 8, 10, 12]);

legend(legend_cell,  'location', 'best', 'interpreter', 'latex', 'NumColumns', 2);

%Export
exportgraphics(h_fig, 'RIR_RT60_controlled.png');
