%Plot image-source model RIR and spectrogram

%Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Input
%h:     [M x 1] RIR 
%options: struct

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%h_fig:     Handle to figure
%h_edc:     Handle to EDC curve

function [h_fig, h_edc] = plot_ISM_RIR(h, options)

arguments
    h = [1];

    options.Fs = 48000;
    options.enable_disp = false;
    options.win_size = 512;
    options.N_FFT = 512;
    options.clim = [-120, -40];
    options.spectrogram_scale {mustBeMember(options.spectrogram_scale, {'linear', 'log'})}  = 'linear';
    options.fig_size = [900 600] * (3/4);
    options.font_size = 16;

    options.name = '';
    options.colormap = parula;
    options.legend_location = 'east';

    options.RT60_dB_hi = -10;
    options.RT60_dB_lo = -30;

    options.disp_EDC_fig = false;

end

h_real = real(h);
h_imag = imag(h);

if all(isreal(h)) %RIR h is all reals
    num_cols = 1;

    name_real = options.name;
    name_imag = options.name;
    
else
    num_cols = 2;

    name_real = [options.name, ' Real'];
    name_imag = [options.name, ' Imag'];
end

%Plotting
fontsize = options.font_size;

%Time-domain
h_fig = figure;
h_fig.Position = [100, 100, options.fig_size(1) * num_cols, options.fig_size(2)];
tiledlayout(2, num_cols,  'TileIndexing', 'columnmajor'); 


for nc = 1:num_cols
    if nc == 1
        h = h_real;
        name = name_real;
    else
        h = h_imag;
        name = name_imag;
    end

    %RIR
    nexttile;
    [dst, krt, stdw] = estimate_RIR_density(h, options.Fs, options.win_size);
    edc = flipud(cumsum(flipud(h).^2)/sum(h.^2));
    edc_dB = 10 * log10(edc);

    idx_p50 = find(edc <= 0.5, 1);

    N_h = numel(h);
    % idx_lo = ceil(N_h * options.T_frac_lo);
    % idx_hi = ceil(N_h * options.T_frac_hi);

    idx_hi = find(edc_dB <= options.RT60_dB_hi, 1);
    idx_lo = find(edc_dB <= options.RT60_dB_lo, 1);   

    t = (0:numel(h)-1)' / options.Fs;
    t_ms = t * 1000;

    % slope = (options.RT60_dB_hi - options.RT60_dB_lo) / (t(idx_hi) - t(idx_lo));
    % RT60 = -60 / slope;

    [p, stats] = polyfit(t(idx_hi:idx_lo), edc_dB(idx_hi:idx_lo), 1);
    RT60 = -60 / p(1);
    
    
    yyaxis left;
    plot(t_ms, h, '-', t_ms, edc, '-.', t_ms(idx_p50), 0, 'r*', 'linewidth', 1.5);
    grid on;  axis tight;
    xlabel('Time (ms)', 'fontsize', fontsize);
    ylabel('Amplitude', 'fontsize', fontsize);
    ylim([-inf, max(h) * 1.1])
    
    yyaxis right;
    plot(t_ms, dst, 'linewidth', 1.5);
    ylim([0, max([1, max(dst)])])
    ylabel('Echo Density Profile', 'fontsize', fontsize);
    
    %Legend
    h_lg = legend('RIR', 'EDC', 'P50', 'Echo Density', 'location', options.legend_location, 'NumColumns', 2, BackgroundAlpha=.7);
    set(h_lg, 'fontsize', fontsize - 1);
    title(h_lg, ['RT60: ', num2str(RT60, 3), ' s, Rsq: ', num2str(stats.rsquared, 3)]);
    
    title([name, ' RIR'], 'fontsize', fontsize + 1);
    set(gca, 'fontsize', fontsize - 1);
    
    %Spectrogram
    nexttile;

    spectrogram(h, options.win_size, ceil(options.win_size * 0.9), options.N_FFT, options.Fs, "power", "yaxis");
    axis tight; 
    set(gca, 'YScale', options.spectrogram_scale);
    set(gca, 'fontsize', fontsize - 1);
    cb = colorbar; ylabel(cb, 'Power (dB)','FontSize', fontsize); title([name, ' Spectrogram'], 'fontsize', fontsize + 1);
    clim(options.clim); 
    
    if ~isempty(options.colormap)
        colormap(options.colormap)
    end



end

%EDC plot
if options.disp_EDC_fig

    h_edc = figure; 
    h_edc.Position = [300, 100, options.fig_size(1) * num_cols, options.fig_size(2)];
    tiledlayout(1, num_cols,  'TileIndexing', 'columnmajor'); 


    for nc = 1:num_cols
        if nc == 1
            h = h_real;
            name = name_real;
        else
            h = h_imag;
            name = name_imag;
        end

        %RIR
        nexttile;
        [dst, krt, stdw] = estimate_RIR_density(h, options.Fs, options.win_size);
        edc = flipud(cumsum(flipud(h).^2)/sum(h.^2));
        edc_dB = 10 * log10(edc);
    
        idx_p50 = find(edc <= 0.5, 1);
    
        N_h = numel(h);
        % idx_lo = ceil(N_h * options.T_frac_lo);
        % idx_hi = ceil(N_h * options.T_frac_hi);
    
        idx_hi = find(edc_dB <= options.RT60_dB_hi, 1);
        idx_lo = find(edc_dB <= options.RT60_dB_lo, 1);   
    
        t = (0:numel(h)-1)' / options.Fs;
        t_ms = t * 1000;
    
        plot(t_ms, 10*log10(edc), t_ms(idx_lo), 10*log10(edc(idx_lo)), 'ro', t_ms(idx_hi), 10*log10(edc(idx_hi)), 'ko', 'linewidth', 1.5); grid on; axis tight;
        grid on; axis tight;
        xlabel('Time (ms)', 'fontsize', fontsize);
        ylabel('dB', 'fontsize', fontsize);
        title([name, ' Energy Decay Curve'], 'fontsize', fontsize + 1);

    end

else
    h_edc = [];
end

