clear; clc;

path = '09_onbox.txt';
[driverRadius,RE,VT,Vgen,frequency,magnitude,phase] = txtParser(path);
figure; 


% interp_freq = 1:10
% vq2 = interp1(frequency, maginitude,v,xq,'spline');

% semilogx(frequency, magnitude);
loglog(frequency, magnitude);