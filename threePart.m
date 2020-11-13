clear; clc;

%% ============================
% import offbox data
% ============================
path = '09_offbox.txt';
[driverRadius,RE,VT,Vgen,frequency,magnitude,phase] = txtParser(path);

%% ============================
% small signal parameters and some others
% ============================
fs = 17.4392;
QES = 0.2440;
QMS = 11.6136;
QTS = 0.2390;
VAS = 0.3326;

RES = 276.5016;

% inductors
Le = 0.1300;
LE = 0.0100;
n = 0.6510;

rho_0 = 407/345;
c = 345;

kM = sqrt(1 + 0.2699 * (2*pi*fs)^2*VAS / (c^2*driverRadius));

fs = fs / kM;
QES = QES * kM;
QMS = QMS * kM;
QTS = QTS * kM;

%% ============================
% import offbox data
% ============================
SD = pi * driverRadius^2;

CAS = VAS / (rho_0 * c^2);
CMS = CAS / SD^2;
MAS = 1 / ((2*pi*fs)^2 * CAS);
MMS = MAS * SD^2;

MA1 = 8*rho_0 / (3*pi^2*driverRadius);
MAD = MAS - 2*MA1;
MMD = MAD * SD^2;

RAT = 1/QTS * sqrt(MAS/CAS);
RAE = 1/QES * sqrt(MAS/CAS);
RAS = RAT - RAE;
RMS = RAS * SD^2;

Bl2 = RE / QES * sqrt(MMS/CMS);

MMD
RMS
CMS
Bl2

