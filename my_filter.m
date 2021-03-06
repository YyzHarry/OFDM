function Hd = my_filter
%MY_FILTER Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 8.3 and the Signal Processing Toolbox 6.21.
% Generated on: 28-May-2017 23:24:43

% Equiripple Lowpass filter designed using the FIRPM function.

% All frequency values are in MHz.
Fs = 24;  % Sampling Frequency

N     = 80;     % Order
Fpass = 3.5;    % Passband Frequency
Fstop = 4.5;    % Stopband Frequency
Wpass = 0.005;  % Passband Weight
Wstop = 0.001;  % Stopband Weight
dens  = 20;     % Density Factor

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, [0 Fpass Fstop Fs/2]/(Fs/2), [1 1 0 0], [Wpass Wstop], ...
           {dens});
Hd = dfilt.dffir(b);

% [EOF]
