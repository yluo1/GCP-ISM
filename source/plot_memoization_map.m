%Plot memoization density map

%Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_memoization_map

fontsize = 18;
markersize = 12;
linewidth = 1.5;

q = 1:16; 
k = sqrt(q);

k_lambda_2 = sqrt(q/4);

k_lambda_3 = sqrt(q/9);

k_lambda_4 = sqrt(q/16);

%%%%%%%%%%%
h_fig = figure;
h_fig.Position = [100, 100, 1200/2, 600/2];

loglog(q, k, 'o-', q, k_lambda_2, 'd--', q, k_lambda_3, 's-.',  q, k_lambda_4, '*:', ...
    'MarkerSize', markersize, 'linewidth', linewidth);

xlabel('Area $q$ (Square Meters)', 'fontsize', fontsize, 'interpreter', 'latex');
ylabel('$k$-Radius (Meters)', 'fontsize', fontsize, 'interpreter', 'latex');
title('Memoization Map Under Scaling', 'fontsize', fontsize + 1, 'interpreter', 'latex');

grid on; axis tight;
set(gca, 'fontsize', fontsize - 1);

yticks([floor(min(k)) : ceil(max(k))]);

xticks(2.^(log2(floor(min(q))) : log2(ceil(max(q))) ));

h_lg = legend('1', '2', '3', '4',  'location', 'best', 'Orientation', 'horizontal');
h_lg.Title.String = 'Scale Factor \lambda';
set(h_lg, 'fontsize', fontsize - 1);

exportgraphics(h_fig, 'memoization_map.png');