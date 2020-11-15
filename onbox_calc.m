clear; clc;

%% ============================
% import offbox data
% ============================
path = '09_onbox.txt';
[driverRadius,RE,VT,Vgen,frequency,magnitude,phase] = txtParser(path);

%% ============================
% offbox calc
% ============================

% fs and RES
fct = 60.4207;
RES = 82.209 - RE;

R1 = sqrt((RE+RES) * RE); % R1 == 21.8548

% find points on the plot, linear interpolation
% f1 and f2
x1 = 49.33643;
y1 = 21.31834;
x2 = 50.58398;
y2 = 23.86287;
f1ct = (x1-x1) * (R1-y1) / (y2-y1) + x1;


x1 = 68.58936;
y1 = 024.13844;
x2 = 70.33569;
y2 = 020.47655;
f2ct = (x1-x1) * (R1-y1) / (y2-y1) + x1;

% f1ct f2ct by poly
f1ct = 49.609922103605180; 
f2ct = 69.611434609230460;

% small signal parameters
QMCT = fct / (f2ct-f1ct) * sqrt((RE+RES)/RE);
QECT = QMCT * RE/RES;
QTCT = QMCT*QECT/(QMCT+QECT);

%% ===========================
% calculate VAS
% ===========================

% fs and QES are from offbox
fs = 17.4392;
% QES = 0.2440;
QES = 0.2419;

VAS = VT * (fct/fs * QECT/QES - 1);

% %% ============================
% % draw offbox data
% % ============================
% figure; 
% loglog(frequency, magnitude);
% 
% hold on;
% 
% %% ============================
% % draw offbox simulation
% % ============================
% 
% % manually tweaked parameters
% LE = 0.010;
% Le = 0.13;
% n = 0.651;
% 
% Zmot = ZmotMod(frequency, fct, RE, QMCT, QECT);
% L = LMod(frequency, LE, Le, n);
% loglog(frequency, abs(RE+Zmot+L), 'r');


