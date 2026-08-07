%Check for valid ISM inputs

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%s:             [1 x N] Source offset from origin (meters)
%r:             [1 x N] Receiver offset from origin (Meters)
%l:             [1 x N] Orthotope dimensions (meters)
%lambda:        Scalar, integer scaling

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%pass:          Logical, true if all sources are inside orthotope, orthotope dimensions are positive, false otherwise
%pass_LUT:      Logical, true if all coordinates are integers, false otherwise

function [pass, pass_LUT] = check_ISM_inputs(s, r, l, lambda)

arguments
    s (1,:) double = [1 2 1];
    r (1,:) double = [2 1 1];
    l (1,:) double = [5 6 3];
    lambda (1,1) double {mustBeInteger, mustBePositive} = 1;
end

s = lambda * s;
r = lambda * r;
l = lambda * l;

pass = all([(abs(s) < l/2),  (abs(r) < l/2), (l > 0)]);

pass_LUT = all([(ceil(s) - s == 0), (ceil(r) - r == 0), (ceil(l) - l == 0)]) && pass;