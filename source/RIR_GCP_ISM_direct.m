%Compute Gauss circle problem image-source model RIR from N-dimensional lattice

%Method: Direct evaluation of GCP-ISM volumes (slow methods), no LUT, forward construction only

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input

%T:             Scalar, time (sec)

%s:             [1 x N] Source offset from origin (meters)
%r:             [1 x N] Receiver offset from origin (Meters)
%l:             [1 x N] Orthotope dimensions (meters)

%gamma_pos:     [1 x N] Reflection coefficient for wall on +axis
%gamma_neg:     [1 x N] Reflection coefficient for wall on -axis

%options:       struct
%options.mode:  Evaluation method,  'direct'        Compute all image-sources within evaluation distances
%                                   'recursive'     Recursive evaluation of S
%options.Fs:            Sample rate
%options.mode:          String, counting method {'direct', 'recursive'}
%                               'direct'            Generate all lattice points within distance ball
%                               'recursive':        Compute from sub-solutions of lower dimension
%options.enable_disp:   Logical, if true, display RIR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%h:        [ceil(T * Fs) x 1] RIR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate and compare direct GCP_ISM with ISM

% gamma_pos = [0.93, 0.8, 0.9];
% gamma_neg = [0.72, 0.78, 0.8];
% T = 0.2;
% s = [1, 2, 1];
% r = [2, 1, 1];
% l = [5, 6, 3];

%h_ref = RIR_ISM_direct(T, s, r, l, gamma_pos, gamma_neg, 'enable_disp', true, 'ceil_sample_dist', true);
%h = RIR_GCP_ISM_direct(T, s, r, l, gamma_pos, gamma_neg, 'enable_disp', true);
%norm(h_ref - h)
%mag2db(norm(h_ref - h) / norm(h_ref))
% figure; plot((0:numel(h_ref)-1)/48000 * 1000, h_ref - h, 'o-', 'linewidth', 1.25); 
% grid on; axis tight;  xlabel('ms'); ylabel('error');

function [h] = RIR_GCP_ISM_direct(T, s, r, l, gamma_pos, gamma_neg, options)

arguments
    T (1,1) double {mustBeNonnegative} = 0.2;
    s (1,:) double = [1 2 1];
    r (1,:) double = [2 1 1];
    l (1,:) double = [5 6 3];
    gamma_pos (1,:) double = [0.93, 0.8, 0.9];
    gamma_neg (1,:) double = [0.72, 0.78, 0.8];

    %Options struct
    options.Fs (1,1) double {mustBeNonnegative} = 48000;
    options.mode (1,:) char {mustBeMember(options.mode, {'direct', 'recursive'})}  = 'recursive';
    
    %Display options
    options.enable_disp (1,1) logical = false;
    
end

N = numel(s);

c = 343; %Speed of sound
M_h = ceil(T * options.Fs);  %Number of taps
h = zeros(M_h, 1);

%Compute distances
k = (0:M_h) / options.Fs * c; %Meters
q = k.^2; %Squared meters

%Select mode
if strcmp(options.mode, 'direct') 

    S_func = @(q) GCP_ISM_direct(q, s, r, l, gamma_pos, gamma_neg);
    name = 'Direct GCP-ISM';

elseif strcmp(options.mode, 'recursive')

    S_func = @(q) GCP_ISM_recur(q, s, r, l, gamma_pos, gamma_neg);
    name = 'Recursive GCP-ISM';

else
    error('Unknown mode');
end

%Iterate over distances
S_prev = 0;
for n = 2:M_h+1
    S_curr = S_func(q(n));
    h(n-1) = (S_curr - S_prev) / (k(n)^((N - 1)/2));
    S_prev = S_curr;
end

%Plotting
if options.enable_disp
    plot_ISM_RIR(h, 'Fs', options.Fs, 'name', name);
end
