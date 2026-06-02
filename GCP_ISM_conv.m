%Compute weighted image-source contributions in R^N within distance sqrt(q)
%via weghted Gauss Circle Problem centered on receiver coordinate r,
%source coordinate s, and orthotope with dimensions l, reflection coefficients
%gamma_pos for wall along +axis, gamma_neg for wall along -axis.

%Convolution method: Compute S(q, N) from sub-solutions of weighted S(m, N-1)
%expressed in convolution form. Requires integer Q, s, r, l

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input

%Q:             Scalar, maximum squared distance (meters^2)

%s:             [1 x N] Source offset from origin (meters)
%r:             [1 x N] Receiver offset from origin (Meters)
%l:             [1 x N] Orthotope dimensions (meters)

%gamma_pos:     [1 x N] Reflection coefficient for wall on +axis
%gamma_neg:     [1 x N] Reflection coefficient for wall on -axis

%mode:      String, counting method,    'padded_FFT' (Default)        %O(NQ) memory, O(N Q log(Q)) compute
%                                       'trunc_FFT'                   %O(2Q) memory, O(N Q log(Q)) compute
%                                       'conv'                        %O(NQ) memory, O(N Q^2)      compute

%jitter_coord_bnd:   [1 x 2]    Jitter the image source coordinates by 
%                               unifrnd(min(jitter_coord_bnd), max(jitter_coord_bnd))
%                               (Default = [0, 0] is disabled)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%S:     [(Q+1) x N] LUT of image-source volumes S(p, N) for squared distances (meters) <= p

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:

% gamma_pos = [0.5, 0.3, 0.9];
% gamma_neg = [0.4, 0.2, 0.8];
% Q = 1000;
% s = [1, 2, 1];
% r = [2, 1, 1];
% l = [5, 6, 3];
% 
% tic; S_DP = GCP_ISM_DP(Q, s, r, l, gamma_pos, gamma_neg); toc
% tic; S_conv = GCP_ISM_conv(Q, s, r, l, gamma_pos, gamma_neg, 'conv'); toc
% tic; S_padded_FFT = GCP_ISM_conv(Q, s, r, l, gamma_pos, gamma_neg, 'padded_FFT'); toc
% tic; S_trunc_FFT = GCP_ISM_conv(Q, s, r, l, gamma_pos, gamma_neg, 'trunc_FFT'); toc

% err = norm(S_DP - S_conv)
% err = norm(S_DP - S_padded_FFT)
% err = norm(S_DP - S_trunc_FFT)


function S = GCP_ISM_conv(Q, s, r, l, gamma_pos, gamma_neg, mode, jitter_coord_bnd)

if nargin < 7 || isempty(mode)
    mode = 'padded_FFT';
end

if nargin < 8 || isempty(jitter_coord_bnd)
    jitter_coord_bnd = [0, 0];
end


N = numel(s);
S = zeros([(Q+1), N]);

sqrt_Q = sqrt(Q);

%Compute LUT basecase
if isequal(jitter_coord_bnd, [0, 0]) %No jitter
    S(:, 1) = GCP_ISM_basecase(Q, s(1), r(1), l(1), gamma_pos(1), gamma_neg(1), 'geometric');

else %Add jitter
    S(:, 1) = GCP_ISM_basecase(Q, s(1), r(1), l(1), gamma_pos(1), gamma_neg(1), 'geometric', jitter_coord_bnd);    
end

if strcmp(mode, 'padded_FFT')

    S_spec = complex(zeros(N*Q, N));
    S_spec(:, 1) = fft(S(:, 1), N * Q);

elseif strcmp(mode, 'trunc_FFT')

    S_spec = fft(S(:, 1), 2 * Q);

end

%Iterate over higher dimensions
for n = 2:N     

    [bq] = g_upper((r(n) + sqrt_Q) / l(n), s(n) / l(n));
    [aq] = g_lower((r(n) - sqrt_Q) / l(n), s(n) / l(n));

    if bq >= aq
        m_up = 0:bq;
        m_lo = aq:-1;
    
        f_up = zeros([Q+1, 1]);
        f_lo = zeros([Q+1, 1]);

        if isequal(jitter_coord_bnd, [0, 0]) %No jitter

            tau_up = u_vector(m_up, s(n), r(n), l(n)).^2;
            tau_lo = u_vector(m_lo, s(n), r(n), l(n)).^2;

            f_up(tau_up + 1) = gamma_pos(n).^abs(floor((m_up+1)/2)) .* gamma_neg(n).^abs(floor((1-m_up)/2));
            f_lo(tau_lo + 1) = gamma_pos(n).^abs(floor((m_lo+1)/2)) .* gamma_neg(n).^abs(floor((1-m_lo)/2));

        else %Add jitter

            rand_up_LUT = unifrnd(min(jitter_coord_bnd), max(jitter_coord_bnd), size(m_up));
            rand_up_LUT(1) = 0; %Exclude orignal source            
            rand_lo_LUT = unifrnd(min(jitter_coord_bnd), max(jitter_coord_bnd), size(m_lo));
            
            tau_up = round((u_vector(m_up, s(n), r(n), l(n)) + rand_up_LUT  ).^2);
            tau_lo = round((u_vector(m_lo, s(n), r(n), l(n)) + rand_lo_LUT  ).^2);

            %Bounded
            f_up(min(tau_up + 1, Q + 1)) = gamma_pos(n).^abs(floor((m_up+1)/2)) .* gamma_neg(n).^abs(floor((1-m_up)/2));
            f_lo(min(tau_lo + 1, Q + 1)) = gamma_pos(n).^abs(floor((m_lo+1)/2)) .* gamma_neg(n).^abs(floor((1-m_lo)/2));

        end

            
        if strcmp(mode, 'conv')

            S_n = conv((f_up + f_lo), S(:, n-1));
            S(:, n) = S_n(1:(Q+1));

        elseif strcmp(mode, 'padded_FFT')

            F = fft(f_up + f_lo, N * Q);
            S_spec(:, n) = F .* S_spec(:, n-1);
            S_n = ifft(S_spec(:, n));
            S(:, n) = S_n(1:(Q+1));

        elseif strcmp(mode, 'trunc_FFT')

            F = fft(f_up + f_lo, 2 * Q);
            S_n = ifft(F .* S_spec);
            S(:, n) = S_n(1:(Q+1));
            S_spec = fft(S(:, n), 2 * Q);

        else
            error('Unknown mode');
        end

    end

end