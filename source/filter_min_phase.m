%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Minimum phase reconstruction with oversampling of linear spaced magnitude spectra

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%X_mag:         [N x 1] or [1 x N] DC to Fs bins
%os:            Oversample factor
%enable_disp:   Logical, if true display frequency response at 48 kHz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%yhat:          [os * N x 1] time-domain filter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage

%X_mag_oneside = db2mag([-3, -2, -1, -2, -2.5, -3, -6, -8, -12] / 100)';
%yhat = filter_min_phase([X_mag_oneside; conj(flipud(X_mag_oneside(2:end-1)))], 1, true);

function yhat = filter_min_phase(X_mag, os, enable_disp)

if nargin < 2 || isempty(os)
    os = 1;
end

if nargin < 3 || isempty(enable_disp)
    enable_disp = false;
end

%Triangle kernel (Linear interpolation)
X_mag = X_mag(:);
n = numel(X_mag);
X_mag2 = conv(upsample(X_mag, os), triang(1 + 2 * (os-1)));
X_mag = X_mag2(os:((1 + n)*os - 1));
nq = floor(os * n / 2) + 1;
X_mag(nq : end) = flipud(X_mag(2:nq));

n = numel(X_mag);
xhat = real(ifft(log(X_mag)));
odd = fix(rem(n,2));
wn = [1; 2 * ones([(n+odd)/2-1,1]) ; ones([1-rem(n,2),1]); zeros([(n+odd)/2-1,1])];
yhat = real(ifft(exp(fft(wn.*xhat))));

%Plotting
if enable_disp
    fvtool(yhat, 1, 'fs', 48000);
end