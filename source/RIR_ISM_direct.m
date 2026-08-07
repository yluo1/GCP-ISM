%Compute image-source model RIR from N-dimensional lattice coordinates.

%options.ceil_sample_dist == true  is equivalent to forward method in GCP-ISM
%options.ceil_sample_dist == false is equivalent to inverse method in GCP-ISM

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

%options.Fs:                    Sample rate
%options.ceil_sample_dist:      Logical, if true round up image-source delay to integer sample
%options.kernel_sample_width:   Half-window size of Lanczos kernel, must be non-negative integer
%options.jitter_coord_bnd:      [1 x 2]    Jitter the image source coordinates by 
%                               unifrnd(min(jitter_coord_bnd), max(jitter_coord_bnd))
%                               (Default = [0, 0] is disabled)
%options.jitter_srand:          Random seed for jitter

%options.enable_disp:           Logical, if true, display RIR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%h:        [ceil(T * Fs) x 1] RIR
%h_fig:    Handle to figure, [] if emable_disp is false

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:

% gamma_pos = [0.93, 0.8, 0.9];
% gamma_neg = [0.72, 0.78, 0.8];
% T = 0.2;
% s = [1, 2, 1];
% r = [2, 1, 1];
% l = [5, 6, 3];

%h = RIR_ISM_direct(T, s, r, l, gamma_pos, gamma_neg, 'enable_disp', true);
%h_rounded = RIR_ISM_direct(T, s, r, l, gamma_pos, gamma_neg, 'ceil_sample_dist', true, 'enable_disp', true);

function [h, h_fig] = RIR_ISM_direct(T, s, r, l, gamma_pos, gamma_neg, options)

arguments
    T (1,:) double = 0.2;
    s (1,:) double = [1 2 1];
    r (1,:) double = [2 1 1];
    l (1,:) double = [5 6 3];
    gamma_pos (1,:) double = [0.93, 0.8, 0.9];
    gamma_neg (1,:) double = [0.72, 0.78, 0.8];

    %Options struct
    options.Fs (1,1) double {mustBeNonnegative} = 48000;
    options.ceil_sample_dist (1,1) logical = false;
   
    options.kernel_sample_width (1,1) double {mustBeNonnegative} = 10;

    options.jitter_coord_bnd (1,2) double = [0, 0];
    options.jitter_srand (1,1) double {mustBeInteger, mustBeNonnegative} = 6452;

    %Display options
    options.enable_disp (1,1) logical = false;
    
end

N = numel(s);

c = 343; %Speed of sound
k = T * c;
q = k^2;
Ts = 1 / options.Fs;
c_Ts = c * Ts;

M_h = ceil(T * options.Fs);  %Number of taps
h = zeros(M_h, 1);


%Compute bounds
ub = zeros(1, N);
lb = zeros(1, N);
x_list = cell(1, N);
for n = 1:N
    while u_vector(ub(n), s(n), r(n), l(n))^2 < q
        ub(n) = ub(n) + 1;
    end
    while u_vector(lb(n), s(n), r(n), l(n))^2 < q
        lb(n) = lb(n) - 1;
    end
    x_list{n} = lb(n):ub(n);
end

%Compute lattice coordinates
V_set = table2array(combinations(x_list{:})); %[M x N]

U_set = u_vector(V_set, s, r, l);   %[M x N]
if ~isequal(options.jitter_coord_bnd, [0, 0]) %Jitter the image-source coordinates
    rng(options.jitter_srand);

    rand_off_LUT = unifrnd(min(options.jitter_coord_bnd), max(options.jitter_coord_bnd), size(U_set));
    rand_off_LUT(ismember(V_set, zeros(1, N), 'rows'), :) = zeros(1, N); %Exclude original source

    U_set = U_set + rand_off_LUT;   %Add random offset to image source coordinates
end
R_sq = sum(U_set.^2, 2);            %[M x 1]
D_dist = sqrt(R_sq);

mask_valid = (R_sq <= q); %[M x 1]

%Compute reflections per valid image-source
R_set = double(mask_valid);
for n = 1:N
    R_set(mask_valid) = R_set(mask_valid) .* ... 
        gamma_pos(n).^abs(floor((V_set(mask_valid, n) + 1)/2) )  .* ...
        gamma_neg(n).^abs(floor((1 - V_set(mask_valid, n))/2) );
end

idx_is = find(mask_valid); %Image-source indices

%Sample fractional indices of valid image-sources
idx_h_frac = D_dist(mask_valid) / c * options.Fs;
idx_h = ceil(idx_h_frac); %Rounded up

%Sum of image-source contributions
if options.ceil_sample_dist
    d_atten = 1 ./ ((idx_h ./ options.Fs * c).^((N - 1)/2)); %Attenuate by rounded distance 
    d_atten(~isfinite(d_atten)) = 0; %Set singularities to 0

    for n = 1:numel(idx_h)
        h(idx_h(n)) = h(idx_h(n)) + R_set(idx_is(n)) .* d_atten(n);
    end

else

    d_atten = 1 ./ D_dist(idx_is).^((N - 1)/2); %Attenuate by true image-source distance
    d_atten(~isfinite(d_atten)) = 0; %Set singularities to 0
 
    a = options.kernel_sample_width; %kernel width

    for n = 1:numel(idx_h_frac)

        % k = sqrt(q / lambda^2);
        % i = k / c * options.Fs; %Fractional index
        % 
        % idx = (ceil(i - a) : floor(i + a))';
        % idx(idx <= 0) = []; idx(idx >= M_h) = [];         %Trim at boundaries
        % 
        % k_idx = c_Ts * idx; %Sample distance        
        % h(idx) = h(idx) + (S(q + 1, N) - S(q, N)) * lanczos_kernel(idx - i, a) ./ (k_idx.^((N - 1)/2));  
        
        idx = (ceil(idx_h_frac(n) - a) : floor(idx_h_frac(n) + a))'; %Write sample index
        idx(idx <= 0) = []; idx(idx >= M_h) = [];         %Trim at boundaries
        k_idx = c_Ts * idx; %Sample distance        

        h(idx) = h(idx) + R_set(idx_is(n))  * lanczos_kernel(idx - idx_h_frac(n), a) ./ (k_idx.^((N - 1)/2));

    end
end

%Name
if options.ceil_sample_dist
    name = 'Integer Delay ISM';
else
    name = 'Fractional Delay ISM';
end

%Plotting
if options.enable_disp
    h_fig = plot_ISM_RIR(h, 'Fs', options.Fs, 'name', name);
else
    h_fig = [];
end
