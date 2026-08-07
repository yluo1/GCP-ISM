%Compute weighted image-source contributions in R^N within distance sqrt(q)
%via weghted Gauss Circle Problem centered on receiver coordinate r,
%source coordinate s, and orthotope with dimensions l, reflection coefficients
%gamma_pos for wall along +axis, gamma_neg for wall along -axis.

%Base-case for N = 1 dimensions

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input

%Q:             Scalar, maximum squared distance (meters^2)
%s:             Scalar, source offset from origin (meters)
%r:             Scalar, receiver offset from origin (meters)
%l:             Scalar, orthotope dimensions (meters)

%gamma_pos:     Scalar, reflection coefficient for wall on +axis
%gamma_neg:     Scalar, reflection coefficient for wall on -axis

%mode:          String, counting method,    'direct'
%                                           'geometric' (Default)

%jitter_coord_bnd:   [1 x 2]    Jitter the image source coordinates by 
%                               unifrnd(min(jitter_coord_bnd), max(jitter_coord_bnd))
%                               (Default = [0, 0] is disabled)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%S:     [(Q+1) x 1]  LUT of image-source volumes S(q+1, 1) for squared distances (meters) <= q

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:

% gamma_pos = 0.5;
% gamma_neg = 0.31;
% Q = 10000;
% s = 0.49;
% r = 0.1;
% l = 1;
% 
% GCP_ISM_basecase_direct      = GCP_ISM_basecase(Q, s, r, l, gamma_pos, gamma_neg, 'direct');
% GCP_ISM_basecase_geometric   = GCP_ISM_basecase(Q, s, r, l, gamma_pos, gamma_neg, 'geometric');
% 
% err = norm(GCP_ISM_basecase_direct - GCP_ISM_basecase_geometric)

function S = GCP_ISM_basecase(Q, s, r, l, gamma_pos, gamma_neg, mode, jitter_coord_bnd)

if nargin < 7 || isempty(mode)
    mode = 'geometric';
end

if nargin < 8 || isempty(jitter_coord_bnd)
    jitter_coord_bnd = [0, 0];
end

if coder.target('MATLAB')
    S = zeros(Q+1, 1);
else
    S = complex(zeros(Q+1, 1));
end

if isequal(jitter_coord_bnd, [0, 0])
    use_jitter = false;
else
    use_jitter = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if use_jitter

    sqrt_Q = sqrt(Q);
    [bq] = g_upper((r + sqrt_Q) / l, s / l);
    [aq] = g_lower((r - sqrt_Q) / l, s / l);

    rand_off_LUT = unifrnd(min(jitter_coord_bnd), max(jitter_coord_bnd), abs(bq-aq) + 1, 1);
    rand_off_LUT(1) = 0; %Exclude orignal source  

else

    if ~coder.target('MATLAB')
        sqrt_Q = sqrt(Q);
        [bq] = g_upper((r + sqrt_Q) / l, s / l);
        [aq] = g_lower((r - sqrt_Q) / l, s / l);
        rand_off_LUT = zeros(abs(bq-aq) + 1, 1);
    end
    
end

if strcmp(mode, 'direct')

    if use_jitter
        %With jitter
        for q = 0:Q
            S(q + 1) = GCP_ISM_1D_direct_count_rand(q, s, r, l, gamma_pos, gamma_neg, rand_off_LUT); 
        end
    else
        %Original, no jitter
        for q = 0:Q
            S(q + 1) = GCP_ISM_1D_direct_count(q, s, r, l, gamma_pos, gamma_neg); 
        end
    end

elseif strcmp(mode, 'geometric')

    if use_jitter
        %With jitter
        for q = 0:Q
            S(q + 1) = GCP_ISM_1D_geom_count_rand(q, s, r, l, gamma_pos, gamma_neg, rand_off_LUT);
        end
    else
        %Original, no jitter
        for q = 0:Q      
            S(q + 1) = GCP_ISM_1D_geom_count(q, s, r, l, gamma_pos, gamma_neg);
        end
    end

else

    error('Unknown mode');

end
