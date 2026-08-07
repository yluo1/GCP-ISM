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

%options:                   struct
%options.Fs:                Sample rate
%options.mode:              String, RIR assembly from TFs,  {'ifft', 'direct'}
%                                   'ifft'       IFFT and sum
%                                   'direct'     Direct computation
%options.lambda:            Scalar, scaling factor for coordinates, must be positive integer
%options.jitter_coord_bnd:  [1 x 2]    Jitter the image source coordinates by 
%                           unifrnd(min(jitter_coord_bnd), max(jitter_coord_bnd))
%                           (Default = [0, 0] is disabled)
%options.jitter_srand:      Random seed for jitter

%options.enable_disp:           Logical, if true, display RIR
%options.clim:                  [1 x 2] dB limits for color bar
%options.fig_size:              [1 x 2] figure [width, height] in pixels
%options.font_size:             Scalar, fontsize
%options.legend_location:       String, legend placement

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
    T (1,1) double {mustBeNonnegative} = 0.2;
    s (1,:) double = [1 2 1];
    r (1,:) double = [2 1 1];
    l (1,:) double = [5 6 3];
    gamma_pos (:,:) double = [0.93, 0.8, 0.9];
    gamma_neg (:,:) double = [0.72, 0.78, 0.8];

    %options struct
    options.Fs (1,1) double {mustBeNonnegative} = 48000;
    options.mode {mustBeMember(options.mode, {'ifft', 'direct'})}  = 'ifft';
    
    options.lambda (1,1) double {mustBeInteger, mustBePositive} = 1;
    options.jitter_coord_bnd (1,2) double = [0, 0]; 
    options.jitter_srand (1,1) double {mustBeInteger, mustBeNonnegative} = 6452;

    %Display options
    options.enable_disp (1,1) logical = false;
    options.clim (1,2) double = [-120, -40];
    options.fig_size (1,2) double {mustBePositive} = [900 600] * (3/4);
    options.font_size (1,1) double {mustBePositive} = 16;
    options.legend_location (1,:) char = 'east';

end


P = size(gamma_pos, 1);

if P == 1
    h = real(RIR_GCP_ISM_LUT(T, s, r, l, gamma_pos, gamma_neg, ...
            'lambda', options.lambda, 'jitter_coord_bnd', options.jitter_coord_bnd, 'jitter_srand', options.jitter_srand));
else

    %Assemble frequency response
    M_h = ceil(T * options.Fs);
    H = complex(zeros(M_h, P));    
    for p = 1:P
        H(:, p) = RIR_GCP_ISM_LUT(T, s, r, l, gamma_pos(p, :), gamma_neg(p, :), ...
            'lambda', options.lambda, 'jitter_coord_bnd', options.jitter_coord_bnd, 'jitter_srand', options.jitter_srand);
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
            h(idx + p - 1) = h(idx + p - 1) + real(g(:, p));
        end
        h = h(1:M_h);

    elseif strcmp(options.mode, 'direct')

        %Direct
        g = complex(zeros(M_h, 1));
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
        h = real(g) / P_twoside;
    end
    
end

name = 'Inverse GCP-ISM Freq.';

%Plotting
if options.enable_disp && coder.target('MATLAB')
    plot_ISM_RIR(h, 'Fs', options.Fs, 'name', name, 'clim', options.clim, ...
        'fig_size', options.fig_size, 'font_size', options.font_size, 'legend_location', options.legend_location);
end
