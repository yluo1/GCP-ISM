%Get default RIR parameters

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs

%T:             Scalar, time (sec)

%s:             [1 x N] Source offset from origin (meters)
%r:             [1 x N] Receiver offset from origin (Meters)
%l:             [1 x N] Orthotope dimensions (meters)

%gamma_pos:     [1 x N] Reflection coefficient for wall on +axis
%gamma_neg:     [1 x N] Reflection coefficient for wall on -axis

function [s_full, r_full, l_full, gamma_pos_full, gamma_neg_full, T] = get_default_RIR_params

 T = 0.3;

%Original
s_full         = [1,    0,      1,  0,      1,       3];
r_full         = [2,    1,      1,  3,      2,       2];
l_full         = [5,    4,      3,  7,      6,       8];
gamma_pos_full = [0.93, 0.8,  0.9,  0.93,   0.77,   0.82];
gamma_neg_full = [0.72, 0.78, 0.93, 0.67,   0.52,   0.7];

%Longer and more reflective first 3 dimensions
% T = 0.5;
% gamma_pos_full(1:3) = round(sqrt(gamma_pos_full(1:3)), 2);
% gamma_neg_full(1:3) = round(sqrt(gamma_neg_full(1:3)), 2);

%Alt
% s_full         = [1,    0,      1,  0,      1,       3];
% r_full         = [2,    1,      0,  3,      2,       2];
% l_full         = [5,    4,      3,  7,      6,       8];
% gamma_pos_full = [0.93, 0.8,  0.9,  0.93,   0.77,   0.82];
% gamma_neg_full = [0.72, 0.88, 0.8,  0.67,   0.52,   0.96];
% 


%Alt with high reflection coefficients
% s_full         = [1,    0,      1,  0,      1,       3];
% r_full         = [2,    1,      1,  3,      2,       2];
% l_full         = [5,    4,      3,  7,      6,       8];
% gamma_pos_full = [0.93, 0.9,  0.92,  0.93,   0.97,   0.92];
% gamma_neg_full = [0.92, 0.98, 0.96,  0.97,   0.92,   0.96];
