%Generate demo samples

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input    

%T:                 Time (seconds)
%lambda:            Coordinate scaling
%P:                 Number of frequency centers from DC-Nyquist in wall-reflection transfer function 
%Fs:                Sample rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output

%h_stereo_short_cell:       [2 x 6]  cell matrix, stereo [h_left, h_right; ...] RIRs for N dimensions (rows) 1 to 6
%h_stereo_medium_cell:      [2 x 6]  cell matrix, stereo [h_left, h_right; ...] RIRs for N dimensions (rows) 1 to 6
%h_stereo_long_cell:        [2 x 6]  cell matrix, stereo [h_left, h_right; ...] RIRs for N dimensions (rows) 1 to 6

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code generation
%codegen('gen_demo_RIRs', '-o', 'gen_demo_RIRs_mex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage: Generate sample RIRs for demo

%T = 0.5;
%lambda = 1;
%%P = 1;         %Frequency-independent
%P = 128;         %Frequency-dependent
%Fs = 48000;

%tic;
%[h_stereo_short_cell, h_stereo_medium_cell, h_stereo_long_cell] = gen_demo_RIRs(T, lambda, P, Fs); 
%duration_mat = toc

%tic;
%[h_stereo_short_cell_mex, h_stereo_medium_cell_mex, h_stereo_long_cell_mex] = gen_demo_RIRs_mex(T, lambda, P, Fs);
%duration_mex = toc

%mkdir('out_wav');
%for n = 1:6
%   y_n = [h_stereo_long_cell_mex{1, n}, h_stereo_long_cell_mex{2, n}];
%   audiowrite(fullfile('out_wav', ['stereo_long_N_', num2str(n), '.wav']), y_n / max(abs(y_n(:))), Fs);
%end
%for n = 1:6
%   y_n = [h_stereo_medium_cell_mex{1, n}, h_stereo_medium_cell_mex{2, n}];
%   audiowrite(fullfile('out_wav', ['stereo_medium_N_', num2str(n), '.wav']), y_n / max(abs(y_n(:))), Fs);
%end
%for n = 1:6
%   y_n = [h_stereo_short_cell_mex{1, n}, h_stereo_short_cell_mex{2, n}];
%   audiowrite(fullfile('out_wav', ['stereo_short_N_', num2str(n), '.wav']), y_n / max(abs(y_n(:))), Fs);
%end

%err = norm(h_stereo_long_cell_mex{1,6} - h_stereo_long_cell{1,6})


function [h_stereo_short_cell, h_stereo_medium_cell, h_stereo_long_cell, Fs] = gen_demo_RIRs(T, lambda, P, Fs)

arguments
    T (1,1) double {mustBeNonnegative} = 0.5;
    lambda (1,1) double {mustBePositive, mustBeInteger} = 1;
    P (1,1) double {mustBePositive} = 1;
    Fs (1,1) double {mustBeNonnegative} = 48000;
end

h_stereo_short_cell = gen_vary_N_stereo_RIR(T, lambda, P, Fs);
h_stereo_medium_cell = gen_vary_N_stereo_RIR(2 * T, lambda, P, Fs);
h_stereo_long_cell = gen_vary_N_stereo_RIR(8 * T, lambda, P, Fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:             
%T:                 Time (seconds)
%lambda:            Coordinate scaling
%P:                 Number of frequency centers from DC-Nyquist in wall-reflection transfer function 
%Fs:                Sample rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%h_GCP_ISM_cell:    [2 x 6]  cell matrix, stereo [h_left, h_right; ...] RIRs for N dimensions (rows) 1 to 6

function h_GCP_ISM_cell = gen_vary_N_stereo_RIR(T, lambda, P, Fs)

arguments
    T (1,1) double {mustBeNonnegative} = 0.5;
    lambda (1,1) double {mustBePositive, mustBeInteger} = 1;
    P (1,1) double {mustBePositive} = 1;
    Fs (1,1) double {mustBeNonnegative} = 48000;
end

% Define source location per dimension (meters)
%s_full         = [-2,   1,     1,  	     -1,   1,	2];
%s_full          = [1,   lambda,     -1,       3,   -1,  -2] / lambda;
s_full         = [1,    0,      1,  0,      1,       -1];
%s_full         = [4,    3,      -4,  3,     4,       -3]; %Far point

% Define receiver location per dimension (meters), lambda-scale to space them close 
%r_full_left    = [2,    -1,      -1,       3,   -1,  -2] / lambda;
%r_full_right   = [0,    -1,      -1,       3,   -1,  -2] / lambda;

%For apart
% r_full_left    = [2,    -1,     -1,       3,    -1,  -1];
% r_full_right   = [-1,   1,      1,       -3,    1,  2];

r_full_left    = [2,    1,      0,  -2,      -1 ,      1];
r_full_right   = [1,    0,      1,  -1,      -2,       0];

% Define room size per dimension (meters)
%l_full         = [5,    4,      3,  	7,     4,	5];
l_full         = [5,    4,      3,  5,      6,       4];
%l_full         = [5,    5,      5,  5,      5,       5];
%l_full         = [9,    9,      9,  9,      9,       9];

% Number of dimensions
N = numel(s_full);

% Define wall reflection coefficients per dimension (magnitude of intensity reflected)
% gamma_pos_full = [-0.93, 0.8,  -0.9,	0.93,   -0.83,	0.9];
% gamma_neg_full = [0.72, -0.78, 0.93,	-0.67,   0.89,	-0.84];

% Specify number of frequency points between DC and Nyquist
w = linspace(0, pi, P + 1)';

% Design two-tap minimum phase filter with target response at DC and Nyquist
h_refl_2 = filter_two_tap_FIR(0, -6, false);
H_refl_2 = freqz(h_refl_2, 1, [w; pi]); H_refl_2 = H_refl_2(1:end-1);


% Compute reflection coefficients to achieve simulation distance RT60 time
T60 = T * ones(1, N);
gamma_pos_mask = (-1).^(1:N);
gamma_neg_mask = (-1).^((1:N) + 1);

xi = 4;
gamma = compute_RT60_gain_from_dim(l_full(1:N), T60, xi);
gamma_pos = gamma .* gamma_pos_mask;
gamma_neg = gamma .* gamma_neg_mask;

% Specify frequency responses per wall
gamma_pos_freq = H_refl_2 * gamma_pos;
gamma_neg_freq = H_refl_2 * gamma_neg;

% gamma_pos_freq = gamma_pos;
% gamma_neg_freq = gamma_neg;


% Number of room dimensions
N = numel(s_full);	

% Image-source coordinates are sampled from a uniform distribution between [jitter_coord_bnd(1), jitter_coord_bnd(2)]
jitter_coord_bnd = [-1e-1, 1e-1]; % Within +- 10 cm
%jitter_coord_bnd = [-2e-1,  2e-1]; % Within +- 20 cm

% Generate RIRs for increasing number of room dimensions
h_GCP_ISM_cell = cell(2, N);

for n = 1:N
	ndims = 1:n;  %Simulate for subset of dimensions  
    
    h_GCP_ISM_cell{1, n} = RIR_GCP_ISM_LUT_freq(T, s_full(ndims), r_full_left(ndims), l_full(ndims), gamma_pos_freq(:, ndims), gamma_neg_freq(:, ndims), ...
    		'mode', 'ifft', 'lambda', lambda, 'jitter_coord_bnd', jitter_coord_bnd, 'Fs', Fs, 'enable_disp',  true);

    h_GCP_ISM_cell{2, n} = RIR_GCP_ISM_LUT_freq(T, s_full(ndims), r_full_right(ndims), l_full(ndims), gamma_pos_freq(:, ndims), gamma_neg_freq(:, ndims), ...
    		'mode', 'ifft', 'lambda', lambda, 'jitter_coord_bnd', jitter_coord_bnd,  'Fs', Fs, 'enable_disp', true);
end