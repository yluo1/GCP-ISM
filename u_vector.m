%Compute vector from receiver to image-source at lattice coordinate

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%v:         [N x M] or [M x N] Lattice coordinate for N dimensions
%s:         [N x 1] or [1 x N] Source coordinate
%r:         [N x 1] or [1 x N] Receiver coordinate
%l:         [N x 1] or [1 x N] Orthotope dimensions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%u:         [N x M] or [1 x N] Vector from receiver to image-source

function u = u_vector(v, s, r, l)

%u = l .* v + ((-1).^v) .* s - r;
u = bsxfun(@times, v, l) + bsxfun(@minus, bsxfun(@times, ((-1).^v), s), r);