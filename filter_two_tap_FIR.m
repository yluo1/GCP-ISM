%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Two-tap min-phase FIR with DC and Nyquist gain control

%Author: Yuancheng Luo, 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%dB_DC:         dB gain at DC
%dB_NQ:         dB gain at Nyquist
%enable_disp:   Logical, if true, plot frequency response assuming 48 kHz Fs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%h:             [1 x 2] FIR coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sample usage
%h = filter_two_tap_FIR(0, -6);
%fvtool(h, 1, 'fs', 48000);

function [h] = filter_two_tap_FIR(dB_DC, dB_NQ, enable_disp)

if nargin < 3 || isempty(enable_disp)
    enable_disp = false;
end

mag_DC = db2mag(dB_DC);
mag_NQ = db2mag(dB_NQ);

mag_sum  = (mag_DC + mag_NQ) / 2;
mag_diff = (mag_DC - mag_NQ) / 2;

if abs(mag_sum) > abs(mag_diff)
    h = [mag_sum, mag_diff];
else
    h = [mag_diff, mag_sum];
end

%Plotting
if enable_disp
    fvtool(h, 1, 'fs', 48000);
end