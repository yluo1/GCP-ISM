%Find min int m where m + ((-1)^m) * y >= x 

%Author: Yuancheng Luo, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Input
%x:         
%y:
%enable_check:      Logical, if true, check against iterative method

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Output
%m:         Index
%f:         Left-hand side f = m + ((-1)^m) * y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:

% enable_check = true;
% [m] = g_lower(10, 0.9, enable_check)
% [m] = g_lower(10, 0.999, enable_check)
% [m] = g_lower(10, -0.999, enable_check)
% [m] = g_lower(10, -0.001, enable_check)
% [m] = g_lower(10, 0.001, enable_check)
% [m] = g_lower(11, 0.001, enable_check)
% [m] = g_lower(11, -0.001, enable_check)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage:
%q = 523;
%r = 1.2;
%l = 3;
%s = -0.3
%[m_max, f_max] = g_upper((sqrt(q) + r)/l, s/l)
%[m_min, f_min] = g_lower((r - sqrt(q))/l, -s/l)

function [m, f] = g_lower(x, y, enable_check)

if nargin < 3 || isempty(enable_check)
    enable_check = false;
end

m = min(2 * ceil((x - y)/2), 2 * ceil((x + y - 1) / 2 ) + 1  );

if  enable_check && coder.target('MATLAB')

    %Iterative search
    m_ref = floor(x - abs(y));
    while m_ref < x -  ((-1)^m_ref) * y
        m_ref = m_ref + 1;
    end

    if m_ref ~= m
        error('g_lower failed');
    end
    
end

f = m + ((-1)^m) * y;