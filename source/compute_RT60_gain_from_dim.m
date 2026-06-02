%Compute reflection coefficent magnitude for obtaining target T60
%given room dimensions, and expected power of stochastic processes
%(e.g. randomized phase-inversion)

%Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%l:                 [1 x N] Room dimensions in meters
%T60:               [1 x N] Decay time for 60 dB attenuation per dimension
%xi:                Expected power of stochastic processes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%gamma:     [1 x N] reflection coefficient magnitude

function gamma = compute_RT60_gain_from_dim(l, T60, xi)

if nargin < 1 || isempty(xi)
    xi = 1;
end

c = 343;
gamma = db2mag(-60 * l ./ (c * T60 .* sqrt(xi)));