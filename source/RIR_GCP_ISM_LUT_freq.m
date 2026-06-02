%Compute Gauss circle problem image-source model RIR from N-dimensional lattice
%and complex reflection coefficients at uniform frequency-spaced w from DC to Nyquist

%Method: Lookup table evaluation of GCP-ISM volumes, 
%Inverse construction: q -> k
%Compute complex solutions H(i, w), 
%find time-domain g(i, p), 
%assemble RIR by summing delayed columns h(i) = sum_p^P g(i-p, p)

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%T:             Scalar, time (sec)

%s:             [1 x N] Source offset from origin (meters)
%r:             [1 x N] Receiver offset from origin (Meters)
%l:             [1 x N] Orthotope dimensions (meters)

%gamma_pos:     [P x N] Complex reflection coefficient for wall on +axis uniformly spanning DC to Nyquist (single-sided spectrum)
%gamma_neg:     [P x N] Complex reflection coefficient for wall on -axis uniformly spanning DC to Nyquist (single-sided spectrum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%h:        [ceil(T * Fs) x 1] RIR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage

% T = 0.3;
% s = [1, 2, 1];
% r = [2, 1, 1];
% l = [5, 6, 3];
% 
% ndims = 1:3
% %ndims = 1:2
% %ndims = 1
% s = s(ndims);
% r = r(ndims);
% l = l(ndims);


%%P = 1;
%%P = 2;
%%P = 4;
%%P = 16;
%P = 17;

%if P == 1
% w = 0;
% else 
%     w = linspace(0, pi, P)';  %Radians span DC to Nyquist [0, pi]
% end

%h_refl = filter_two_tap_FIR(-2, -40);
%H_refl = freqz(h_refl, 1, [w; pi]); H_refl = H_refl(1:end-1);  %[P x 1]

% gamma_pos = repmat(H_refl, [1, numel(ndims)]);
% gamma_neg = repmat(H_refl, [1, numel(ndims)]);

%h = RIR_GCP_ISM_LUT_freq(T, s, r, l, gamma_pos, gamma_neg, 'enable_disp', true);

function h = RIR_GCP_ISM_LUT_freq(T, s, r, l, gamma_pos, gamma_neg, options)

arguments
    T = 0.2;
    s = [1 2 1];
    r = [2 1 1];
    l = [5 6 3];
    gamma_pos = [0.93, 0.8, 0.9];
    gamma_neg = [0.72, 0.78, 0.8];

    %Options struct
    options.Fs = 48000;
    options.mode  {mustBeMember(options.mode, {'ifft', 'direct'})}  = 'ifft';
    options.lambda = 1;

    %Display options
    options.enable_disp = false;
    options.clim = [-120, -40];
    options.fig_size = [900 600] * (3/4);
    options.font_size = 16;
    options.legend_location = 'east';

end


P = size(gamma_pos, 1);

if P == 1
%if false

    h = RIR_GCP_ISM_LUT(T, s, r, l, gamma_pos, gamma_neg);

else

    %Assemble frequency response
    M_h = ceil(T * options.Fs);
    H = zeros(M_h, P);
    for p = 1:P
        H(:, p) = RIR_GCP_ISM_LUT(T, s, r, l, gamma_pos(p, :), gamma_neg(p, :), ...
            'lambda', options.lambda);
    end

    %Force real DC, Nyquist
    % H(:, 1)     = real(H(:, 1));    
    % H(:, end)   = real(H(:, end));
    H = [H, conj(fliplr(H(:, 2:end-1)))];

    if strcmp(options.mode, 'ifft')

        %IFFT version
        g = ifft(H, [], 2, 'symmetric');
        h = zeros(M_h + P - 1, 1);
        idx = (1:M_h)';    
        for p = 1:P
            h(idx + p - 1) = h(idx + p - 1) + g(:, p);
        end
        h = h(1:M_h);

    elseif strcmp(options.mode, 'direct')

        %Direct
        g = zeros(M_h, 1);
        P_twoside = 2 * (P - 1);
        for i = 1:M_h
            for n = 1:P_twoside
                if i - n > 0 && i - n + 1 <= M_h
                    for m = 1:P_twoside
                        if i - m >= 0
                            g(i) = g(i) + H(i - n + 1, m) * exp(1i * 2 * pi * (m-1) * (n-1) / P_twoside);
                        end
                    end
                end
            end
        end
        g = real(g) / P_twoside;
        h = g;
    end
    
end

name = 'Inverse GCP-ISM Freq.';

%Plotting
if options.enable_disp
    plot_ISM_RIR(h, 'Fs', options.Fs, 'name', name, 'clim', options.clim, ...
        'fig_size', options.fig_size, 'font_size', options.font_size, 'legend_location', options.legend_location);
end

