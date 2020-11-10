clear; clc;

path = '09_offbox.txt';
[driverRadius,RE,VT,Vgen,frequency,magnitude,phase] = txtParser(path);
figure; semilogx(frequency, magnitude);