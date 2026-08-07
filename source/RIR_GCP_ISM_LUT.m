%Compute Gauss circle problem image-source model RIR from N-dimensional lattice

%Method: Lookup table evaluation of GCP-ISM volumes, 
%optional direct-evaluation for smaller distances,
%Forward construction: k -> q, take finite differences of volume function at distances between k, k+1, 
%                      place contributions at integer samples in h
%Inverse construction: q -> k, take finite differences of volume function at distances ebtween q, q+1
%                      center contributions on Lanzcos kernels with fractional delay

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input

%T:             Scalar, time (sec)

%s:             [1 x N] Source offset from origin (meters)
%r:             [1 x N] Receiver offset from origin (Meters)
%l:             [1 x N] Orthotope dimensions (meters)

%gamma_pos:     [1 x N] Reflection coefficient for wall on +axis
%gamma_neg:     [1 x N] Reflection coefficient for wall on -axis

%options:                   struct
%options.Fs:                Sample rate
%options.mode:              String, evaluation method,  {'DP', 'conv'}
%                                   'DP'       Dynamic programming memoization
%                                   'conv'     Convolution form
%options.lambda:            Scalar, scaling factor for coordinates, must be positive integer
%options.jitter_coord_bnd:  [1 x 2]    Jitter the image source coordinates by 
%                           unifrnd(min(jitter_coord_bnd), max(jitter_coord_bnd))
%                           (Default = [0, 0] is disabled)
%options.jitter_srand:      Random seed for jitter
%options.T_direct:          Scalar (seconds), perform direct computation for distances (seconds) <= T_direct
%options.direction:         String, RIR construction method, {'forward', 'inverse'}
%                               'forward'  Inter-sample finite difference
%                               'inverse'  Sum of kernels centered at squared-distances, weighted by finite difference 
%options.kernel_sample_width:   Half-window size of Lanczos kernel, must be non-negative integer

%options.error_check_samples:   Logical, if true, check for errors due to image coords falling between ceil(q_n) and q_n

%options.RT60_dB_hi:            dB upperbound of echo decay curve for computing RT60
%options.RT60_dB_lo:            dB lowerbound of echo decay curve for computing RT60

%options.disp_EDC_fig:         Logical, if true, disp echo decay curve

%options.enable_disp:           Logical, if true, display RIR
%options.clim:                  [1 x 2] dB limits for color bar [min, max]
%options.fig_size:              [1 x 2] figure [width, height] in pixels
%options.font_size:             Scalar, fontsize
%options.legend_location:       String, legend placement


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%h:        [ceil(T * Fs) x 1] RIR
%h_fig:    Handle to figure, [] if emable_disp is false

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage

% gamma_pos = [0.93, 0.8, 0.9];
% gamma_neg = [0.72, 0.78, 0.8];
% T = 0.2;
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
% gamma_pos = gamma_pos(ndims);
% gamma_neg = gamma_neg(ndims);
% 
% h_ref = RIR_ISM_direct(T, s, r, l, gamma_pos, gamma_neg, 'enable_disp', true); 
% h = RIR_GCP_ISM_LUT(T, s, r, l, gamma_pos, gamma_neg, 'direction', 'inverse', 'enable_disp', true); 
% norm(h_ref - h) 
% mag2db(norm(h_ref - h) / norm(h_ref))
% 
% h_lambda = RIR_GCP_ISM_LUT(T, s, r, l, gamma_pos, gamma_neg, 'direction', 'inverse', 'lambda', [8], 'enable_disp', true);
% norm(h_ref - h_lambda)
% mag2db(norm(h_ref - h_lambda) / norm(h_ref))  
% figure; plot((0:numel(h_ref)-1)/48000 * 1000, h_ref - h, 'o-', (0:numel(h_ref)-1)/48000 * 1000, h_ref - h_lambda, '*--', 'linewidth', 1.25); 
% grid on; axis tight; legend('Error Single Lambda', 'Error Multilambda'); xlabel('ms'); ylabel('error');


function [h, h_fig] = RIR_GCP_ISM_LUT(T, s, r, l, gamma_pos, gamma_neg, options)

arguments

    T (1,1) double {mustBeNonnegative} = 0.2;
    s (1,:) double = [1 2 1];
    r (1,:) double = [2 1 1];
    l (1,:) double = [5 6 3];
    gamma_pos (1,:) double = [0.93, 0.8, 0.9];
    gamma_neg (1,:) double = [0.72, 0.78, 0.8];

    %options struct
    options.Fs (1,1) double {mustBeNonnegative} = 48000;
    options.mode {mustBeMember(options.mode, {'DP', 'conv'})}  = 'conv';
    
    options.lambda (1,1) double {mustBeInteger, mustBePositive} = 1;
    options.jitter_coord_bnd (1,2) double = [0, 0]; 
    options.jitter_srand (1,1) double {mustBeInteger, mustBeNonnegative} = 6452;

    options.T_direct (1,1) double = -1;

    options.direction {mustBeMember(options.direction, {'forward', 'inverse'})} = 'inverse';

    options.kernel_sample_width (1,1) double {mustBeNonnegative}  = 10;

    options.error_check_samples (1,1) logical = false;
    
    options.RT60_dB_hi (1,1) double = -10;
    options.RT60_dB_lo (1,1) double = -30;
    options.disp_EDC_fig (1,1) logical = false;

    %Display options
    options.enable_disp (1,1) logical = false;
    options.clim (1,2) double = [-120, -40];
    options.fig_size (1,2) double {mustBePositive} = [900 600] * (3/4);
    options.font_size (1,1) double {mustBePositive} = 16;
    options.legend_location (1,:) char = 'east';

end

%Check inputs
[pass, pass_LUT] = check_ISM_inputs(s, r, l, options.lambda);
if ~pass
    error('source or receiver is out-of-bounds');
end
if ~pass_LUT
    error('coordinates are not integers');
end

N = numel(s);

c = 343; %Speed of sound
M_h = ceil(T * options.Fs);  %Number of taps

Ts = 1 / options.Fs;
c_Ts = c * Ts;

if coder.target('MATLAB')
    h = zeros(M_h, 1);
else
    h = complex(zeros(M_h, 1));
end

%Set name
name = 'GCP-ISM';

if strcmp(options.direction, 'forward')
    name = ['Forward ', name];
elseif strcmp(options.direction, 'inverse')
    name = ['Inverse ', name];
end

%Compute distances
k_max = M_h / options.Fs * c;
q_max = k_max^2;

Q = ceil(q_max);
k_direct = options.T_direct * c;

%Compute LUT
lambda = options.lambda;

%Q_lambda = Q; %Uniform resolution
Q_lambda = Q * lambda^2; %Max squared distance unadjusted

%Seed jitter random variable
if ~isequal(options.jitter_coord_bnd, [0, 0])
    rng(options.jitter_srand);
end

if strcmp(options.mode, 'DP')
    
    S = GCP_ISM_DP(Q_lambda, lambda * s, lambda * r, lambda * l, gamma_pos, gamma_neg, lambda * options.jitter_coord_bnd);

elseif strcmp(options.mode, 'conv')

    S = GCP_ISM_conv(Q_lambda, lambda * s, lambda * r, lambda * l, gamma_pos, gamma_neg, 'padded_FFT', lambda * options.jitter_coord_bnd);  

else
    error('Unknown mode');
end

if strcmp(options.direction, 'forward')

    %Iterate over sample distances k
    k = (0:M_h) *  c_Ts; %Meters
    q = k.^2; %Squared meters

    S_prev = 0;
    for n = 2:M_h+1
    
        q_n = lambda^2 * ceil(q(n));
        
        %Check for errors due to image falling between ceil(q_n) and q_n
        if options.error_check_samples
            error_s = GCP_ISM_recur(q(n), s, r, l, gamma_pos, gamma_neg) - S(q_n, N);
            if abs(error_s) > 1e-6            
                warning([num2str(error_s), ' error exceed threshold at n = ', num2str(n)]);
            end
        end
    
        if k(n) <= k_direct  %Direct computation
    
            S_curr = GCP_ISM_recur(q(n), s, r, l, gamma_pos, gamma_neg);
    %       S_curr = GCP_ISM_direct(q(n), s, r, l, gamma_pos, gamma_neg);
    
        else                %Use LUT
            S_curr = S(q_n, N);
        end
        
        %Finite difference
        h(n-1) = (S_curr - S_prev) / (k(n)^((N - 1)/2));
        S_prev = S_curr;
    end

elseif strcmp(options.direction, 'inverse') 

    %Iterate over area q
    a = options.kernel_sample_width; %kernel width

    for q = 1:Q_lambda
        
        k = sqrt(q / lambda^2);
        i = k / c * options.Fs; %Fractional index
        
        idx = (ceil(i - a) : floor(i + a))';              %Write sample indices
        idx(idx <= 0) = []; idx(idx >= M_h) = [];         %Trim at boundaries

        k_idx = c_Ts * idx; %Sample distance  

        if k_idx > k_direct
            h(idx) = h(idx) + (S(q + 1, N) - S(q, N)) * lanczos_kernel(idx - i, a) ./ (k_idx.^((N - 1)/2));  
        end
        
    end

    %Direct method
    if options.T_direct > 0
        h_direct = RIR_ISM_direct(options.T_direct, s, r, l, gamma_pos, gamma_neg, ...
            'kernel_sample_width', options.kernel_sample_width);
        
        idx = 1:min(M_h, numel(h_direct));
        h(idx) = h(idx) + h_direct(idx);
    end
    
end

%Plotting
if options.enable_disp && coder.target('MATLAB')
    %Plot volume function
    %figure; loglog(sqrt((0:Q) / lambda^2) / c * 1000, S, 'linewidth', 1.5); grid on; axis tight; xlabel('Time (ms)'); ylabel('S (Volume)'); legend('location', 'northwest')

    %Plot RIR
    h_fig = plot_ISM_RIR(h, 'Fs', options.Fs, 'name', name, 'clim', options.clim, 'RT60_dB_hi', options.RT60_dB_hi, 'RT60_dB_lo', options.RT60_dB_lo, ...
        'disp_EDC_fig', options.disp_EDC_fig, 'fig_size', options.fig_size, 'font_size', options.font_size, 'legend_location', options.legend_location);
else
    h_fig = [];
end
