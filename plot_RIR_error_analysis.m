%Test GCP errors for various lambda scaling

%Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_RIR_error_analysis

%clear all

[s_full, r_full, l_full, gamma_pos_full, gamma_neg_full, T] = get_default_RIR_params;

N = numel(s_full);

%Increase time and reflection coefficients
T = 0.5;
%T = 0.1;

gamma_pos_full = sqrt(gamma_pos_full);
gamma_neg_full = sqrt(gamma_neg_full);

ndims = 1:4
%ndims = 1:3
%ndims = 1:2

s = s_full(ndims);
r = r_full(ndims);
l = l_full(ndims);
gamma_pos = gamma_pos_full(ndims);
gamma_neg = gamma_neg_full(ndims);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_dim_list = [3, 4];
N_max_dims = numel(max_dim_list);

%Vary lambda
lambda_list = [1 2 4 8 16];
N_lambda = numel(lambda_list);

NMSE = zeros(N_lambda, N_max_dims);
time_sec = zeros(N_lambda, N_max_dims);

h_ref_list = cell(1, N_max_dims);
h_inv_list = cell(N_lambda, N_max_dims);

max_iter = 10;
%max_iter = 1;

%Iterate over max dims
for n = 1:N_max_dims
    N = max_dim_list(n);

    %Reference target
    h_ref_list{n} = RIR_ISM_direct(T, s(1:N), r(1:N), l(1:N), gamma_pos(1:N), gamma_neg(1:N), 'enable_disp', true); 
    
    for i = 1:N_lambda
    
        for j = 1:max_iter
            tic;
            h_inv_list{i, n} = RIR_GCP_ISM_LUT(T, s(1:N), r(1:N), l(1:N), gamma_pos(1:N), gamma_neg(1:N), 'lambda', [lambda_list(i)], 'direction', 'inverse',...
                'enable_disp', false);
            time_sec(i, n) = time_sec(i, n) + toc;
        end
        time_sec(i, n) = time_sec(i, n) / max_iter

        NMSE(i, n) = sum((h_ref_list{n} - h_inv_list{i, n} ).^2) / sum(h_ref_list{n}.^2)
    end


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fontsize = 18;
markersize = 12;
linewidth = 1.25;

h_fig = figure;
h_fig.Position = [100, 100, [600, 300] * 0.75]
yyaxis left;
semilogx(lambda_list, 10 * log10(NMSE(:, 1)), 'linewidth', linewidth, 'Marker', 'o', 'MarkerSize', markersize); hold on;
semilogx(lambda_list, 10 * log10(NMSE(:, 2)), 'linewidth', linewidth, 'Marker', '*', 'MarkerSize', markersize); hold on;
grid on; axis tight;
ylabel('NMSE (dB)', 'fontsize', fontsize-2, 'interpreter', 'latex');
title('RIR Generation Time and Error for $\lambda$ Scaling', 'fontsize', fontsize-2, 'interpreter', 'latex');
set(gca, 'fontsize', fontsize - 3);
xticks(lambda_list);


yyaxis right;
loglog(lambda_list, time_sec(:, 1), 'linewidth', linewidth, 'Marker', 'o','MarkerSize', markersize); hold on;
loglog(lambda_list, time_sec(:, 2), 'linewidth', linewidth, 'Marker', '*','MarkerSize', markersize); hold on;
grid on; axis tight;
xlabel('$\lambda$ Scaling', 'fontsize', fontsize, 'interpreter', 'latex');
ylabel('Runtime (Seconds)', 'fontsize', fontsize, 'interpreter', 'latex');
%yticks([1000 10000]);

legend({'$N = 3$', '$N = 4$', '$N = 3$', '$N = 4$'}, 'location', 'west', 'interpreter', 'latex');

%Export
exportgraphics(h_fig, 'RIR_error_runtime.png');
