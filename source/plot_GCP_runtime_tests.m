%Test GCP runtimes for cases of k, N

%Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_GCP_runtime_tests
%clear all

%N_list = 1:12;             %N dimensions
N_list = 1:12;             %N dimensions

%k = floor(48000 / 2);     %k distance
%k = 1024;     %k distance

%k_list = 1024;
%k_list = 512;
%k_list = 256;
%k_list = 32;

k_list = [256, 1024]; %Original
%k_list = [256, 512]; 
%k_list = [16, 32]; %Test


M_N = numel(N_list);
M_K = numel(k_list);

zero_table = nan([M_K, M_N]);

GCP_direct_sol      = zero_table;
GCP_direct_times    = zero_table;

GCP_direct_recur_sol      = zero_table;
GCP_direct_recur_times    = zero_table;

GCP_DP_sol      = zero_table;
GCP_DP_times    = zero_table;

GCP_conv_sol      = zero_table;
GCP_conv_times    = zero_table;

%Solve
%max_iter = 5;
max_iter = 10;
%max_iter = 1;

for i = 1:M_N %Iterate over N
    N = N_list(i);
  
    for j = 1:M_K %Iterate over K
        k = k_list(j);

        disp(['N = ', num2str(N), ', k = ', num2str(k)])

        if  N <= 3 && k <= 256
            disp('direct');
            tic
            for iter = 1:max_iter       
                GCP_direct_sol(j, i) = GCP_direct(k, N, true);
            end
            GCP_direct_times(j, i)  = toc;
        end

        %Direct recurrence solution
        if N <= 3 || (N <= 4 && k <= 1024) || (N <= 5 && k <= 256)
            disp('direct recur');
            tic
            for iter = 1:max_iter       
                GCP_direct_recur_sol(j, i)  = GCP_direct_recur(k, N);
            end
            GCP_direct_recur_times(j, i) =  toc;
        end
    
        %Dynamic programming solution        
        disp('DP');
        tic
        for iter = 1:max_iter      
           GCP_DP_sol(j, i) = GCP_DP(k, N);
        end
        GCP_DP_times(j, i) = toc;
    
        %Convolution solution
        disp('conv');
        tic
        for iter = 1:max_iter   
            GCP_conv_sol(j, i) = GCP_conv(k, N, 'padded_FFT');
        end
        GCP_conv_times(j, i) = toc;


    end

end


GCP_direct_sol
direct_ms = GCP_direct_times * 1000 / max_iter

GCP_direct_recur_sol;
direct_recur_ms = GCP_direct_recur_times * 1000 / max_iter
isequal(GCP_direct_sol, GCP_direct_recur_sol)

GCP_DP_sol;
DP_ms = GCP_DP_times * 1000 / max_iter
isequal(GCP_direct_sol, GCP_DP_sol)

GCP_conv_sol;
conv_ms = GCP_conv_times * 1000 / max_iter
isequal(GCP_direct_sol, GCP_conv_sol)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fontsize = 20;
marker_list = {'o', '*', '+', 'x', 's', 'd', '^', 'v', '>', '<'};
linestyle_list = {'-', '--', '-.', ':'};

h_fig = figure;
h_fig.Position = [100, 100, [900, 450]*0.75];

legend_str = {};
linewidth = 1.5;
markersize = 12;

for j = 1:M_K

    if ~all(isnan(direct_ms(j, :)))
        loglog( N_list, direct_ms(j, :), 'r', 'linestyle', linestyle_list{j}, 'Marker', marker_list{j}, 'linewidth', linewidth, 'MarkerSize', markersize);    hold on; 
        legend_str{end+1} =  ['$k:', num2str(k_list(j)) , '$ ', 'DirectGrid'];
    end
    if  ~all(isnan(direct_recur_ms(j, :)))
        loglog( N_list, direct_recur_ms(j, :), 'b', 'linestyle', linestyle_list{j}, 'Marker', marker_list{j}, 'linewidth', linewidth, 'MarkerSize', markersize);    hold on;
        legend_str{end+1} =  ['$k:', num2str(k_list(j)) , '$ ', 'DirectRecur.'];
    end
    if  ~all(isnan(DP_ms(j, :)))
        loglog( N_list, DP_ms(j, :), 'g',   'linestyle', linestyle_list{j},    'Marker', marker_list{j}, 'linewidth', linewidth, 'MarkerSize', markersize);    hold on;
        legend_str{end+1} =  ['$k:', num2str(k_list(j)) , '$ ', 'LUTSum'];
    end
    if  ~all(isnan(conv_ms(j, :)))
        loglog( N_list, conv_ms(j, :), 'k', 'linestyle', linestyle_list{j},     'Marker', marker_list{j}, 'linewidth', linewidth, 'MarkerSize', markersize);    hold on;
        legend_str{end+1} =  ['$k:', num2str(k_list(j)) , '$ ', 'LUTConv.'];
    end
end
  

xticks(N_list);
yticks([1, 10, 100, 1000, 10000]);
grid on; axis tight;
xlabel('$N$ Dimensions', 'fontsize', fontsize + 4, 'interpreter', 'latex');
ylabel('Time (ms)', 'interpreter', 'latex', 'fontsize', fontsize + 4);


title(['Gauss Circle Problem Runtimes'], ...
    'fontsize', fontsize + 5, 'interpreter', 'latex');


set(gca, 'fontsize', fontsize - 4);

h_lg = legend(legend_str{:}, ...
    'location', 'southeast', 'Interpreter', 'latex', 'NumColumns', 2);

set(h_lg, 'fontsize', fontsize - 5);

h_lg.BackgroundAlpha = 0.7; 

%Save
exportgraphics(gcf, 'GCP_runtimes.png');
