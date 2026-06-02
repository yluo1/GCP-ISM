%Find max int m where m + ((-1)^m) * y <= x 

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

%enable_check = true;
%[m] = g_upper(10, 0.9, enable_check)
%[m] = g_upper(10, 0.999,  enable_check)
%[m] = g_upper(10, -0.999,  enable_check)
%[m] = g_upper(10, -0.001,  enable_check)
%[m] = g_upper(10, 0.001,  enable_check)
%[m] = g_upper(11, 0.001,  enable_check)
%[m] = g_upper(11, -0.001,  enable_check)

function [m, f] = g_upper(x, y, enable_check)

if nargin < 3 || isempty(enable_check)
    enable_check = false;
end

m = max(2 * floor((x - y)/2), 2 * floor((x + y - 1) / 2 ) + 1  );

if  enable_check && coder.target('MATLAB')

    %Iterative search
    m_ref = ceil(x + abs(y));
    while m_ref > x -  ((-1)^m_ref) * y
        m_ref = m_ref - 1;
    end

    if m_ref ~= m 
        error('g_upper failed');
    end
    
end

f = m + ((-1)^m) * y;