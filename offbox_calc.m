clear; clc;

%% ============================
% import offbox data
% ============================
path = '09_offbox.txt';
[driverRadius,RE,VT,Vgen,frequency,magnitude,phase] = txtParser(path);

%% ============================
% offbox calc
% ============================

% fs and RES
fs = 17.4392;
RES = 282.3116 - RE;

R1 = sqrt((RE+RES) * RE); % R1 ~= 40.5

% find points on the plot, linear interpolation
% f1 and f2
x1 = 12.53955;
y1 = 38.89935;
x2 = 12.87280;
y2 = 42.21905;
f1 = (x1-x2) * (R1-y1) / (y2-y1) + x1;

x1 = 23.0684;
y1 = 41.7896;
x2 = 23.667;
y2 = 38.3151;
f2 = (x1-x2) * (R1-y1) / (y2-y1) + x1;


% small signal parameters
QMS = fs / (f2-f1) * sqrt((RE+RES)/RE);
QES = QMS * RE/RES;
QTS = QMS*QES/(QMS+QES);

%% ============================
% draw offbox data
% ============================
figure; 
loglog(frequency, magnitude);

hold on;

%% ============================
% draw offbox simulation
% ============================

% manually tweaked parameters
LE = 0.010;
Le = 0.13;
n = 0.651;

Zmot = ZmotMod(frequency, fs, RE, QMS, QES);
L = LMod(frequency, LE, Le, n);
loglog(frequency, abs(RE+Zmot+L), 'r');


