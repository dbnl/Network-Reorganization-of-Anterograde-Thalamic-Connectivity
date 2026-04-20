function [tradfilt, thephase, amplitude] = Butterworth_Hilbert_LR(signal, Fs, Fpass)
% signal_filtered = BandpassFilter(signal, Fs, Fpass)
%
% Takes 'signal' and bandpasses it to Fpass frequencies
%
% Arguments
%
% signal - arbitrary signal to be filtered
% Fs - sample frequency
% Fpass - 1x2 vector of F_low and F_high indicating passband

% Lara edited, Andrew Bogard conceived.

Wn_theta = [Fpass(1,1)/(Fs/2) Fpass(1,2)/(Fs/2)]; % normalized by the nyquist frequency

[btheta,atheta] = butter(1,Wn_theta); 

tradfilt = filtfilt(btheta,atheta,signal);
h_signal=hilbert(tradfilt);
thephase=atan2(tradfilt, imag(h_signal));
amplitude=abs(h_signal);
end