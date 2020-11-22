clear; clc;

%% ============================
% import offbox data
% ============================
path = '09_offbox.txt';
[driverRadius,RE,VT,Vgen,frequency,magnitude,phase] = txtParser(path);

%% ============================
% small signal parameters and some others
% ============================
% fs = 17.4392;
% QES = 0.2440;
% QMS = 11.6136;
% QTS = 0.2390;
% VAS = 0.3326;

% from poly
% QES = 0.2419;
% QMS = 11.5115;
% QTS = 0.2369;
% VAS = 0.3222;
QES = 0.241885806934716;
QMS = 11.511500009225150;
QTS = 0.236907773837562;
VAS = 0.322172643567582;
fs = 17.439210000000000;

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

RA1 = 0.4410*rho_0*c / (pi*driverRadius^2);
RA2 = rho_0*c / (pi*driverRadius^2);
CA1 = 5.94*driverRadius^3 / (rho_0 * c^2);

RAT = 1/QTS * sqrt(MAS/CAS);
RAE = 1/QES * sqrt(MAS/CAS);
RAS = RAT - RAE;
RMS = RAS * SD^2;

Bl = sqrt(RE / QES * sqrt(MMS/CMS));

%%
figure; 

% measured ZVC
loglog(frequency, magnitude);
% semilogx(frequency, phase);

hold on;

jomega = 1j*2*pi*frequency;
Z_LE = jomega * LE;
Z_Le = (jomega).^n * Le;

ZE = RE + Z_LE.*Z_Le ./ (Z_LE+Z_Le);

ZM = 1./(jomega * CMS) + jomega*MMD + RMS;

RA1pCA1 = RA1*(1./(jomega*CA1)) ./ (RA1 + (1./(jomega*CA1)));

ZAF = (RA1pCA1+RA2) .* jomega*MA1 ./ (RA1pCA1+RA2 + jomega*MA1);
ZAB = ZAF;

ZVC = ZE + ((Bl)^2/SD^2) ./ (ZM/SD^2 + ZAF + ZAB);

% modelled ZVC
loglog(frequency, abs(ZVC), 'r');
%loglog(frequency, abs((Bl)^2./ZM + ZE), 'g');
%semilogx(frequency, angle(phase)/pi*180, 'r');

%% Infinite-baffle

% simulated ZVC with inductance
A = readmatrix("dp-part-2-ib-ZVC.csv");
ms_f = A(:, 1);
ms_ZVC = A(:, 2);

hold on;
loglog(ms_f, abs(ms_ZVC), 'r');

% % modelled pD without inductance
% Gs = (jomega/(2*pi*fs)).^2 ./ ((jomega/(2*pi*fs)).^2 + 1/QTS * jomega/(2*pi*fs) + 1);
% pD = rho_0/(2*pi) * Bl/(SD*RE*MAS) * Gs;
% 
% hold on;
% loglog(frequency, abs(pD), 'r');

% simulated pD without inductance
A = readmatrix("dp-part-2-ib-UD-noL.csv");
ms_f = A(:, 1);
ms_UD = A(:, 2);
ms_jomega = 2*pi*ms_f;
ms_p = rho_0 / (2*pi) * ms_jomega .* ms_UD;

hold on;
loglog(ms_f, abs(ms_p), 'b');


% % modelled pD with inductance
% omega_u1 = RE*MAS/(MAD*LE);
% f_u1 = omega_u1 / (2*pi);
% Ts = 1 ./ (1 + jomega/omega_u1);
% pD = rho_0/(2*pi) * Bl/(SD*RE*MAS) * Gs .* Ts;
% 
% hold on;
% loglog(frequency, abs(pD), 'r');

% simulated pD with inductance
A = readmatrix("dp-part-2-ib-UD.csv");
ms_f = A(:, 1);
ms_UD = A(:, 2);
ms_jomega = 2*pi*ms_f;
ms_p = rho_0 / (2*pi) * ms_jomega .* ms_UD;

hold on;
loglog(ms_f, abs(ms_p), 'b');

%% Closed-box

% simulated ZVC with inductance
A = readmatrix("dp-part-2-cb-ZVC.csv");
ms_f = A(:, 1);
ms_ZVC = A(:, 2);

hold on;
loglog(ms_f, abs(ms_ZVC), 'r');

% simulated pD without inductance
A = readmatrix("dp-part-2-cb-UD-noL.csv");
ms_f = A(:, 1);
ms_UD = A(:, 2);
ms_jomega = 2*pi*ms_f;
ms_p = rho_0 / (2*pi) * ms_jomega .* ms_UD;

hold on;
loglog(ms_f, abs(ms_p), 'b');

% simulated pD with inductance
A = readmatrix("dp-part-2-cb-UD.csv");
ms_f = A(:, 1);
ms_UD = A(:, 2);
ms_jomega = 2*pi*ms_f;
ms_p = rho_0 / (2*pi) * ms_jomega .* ms_UD;

hold on;
loglog(ms_f, abs(ms_p), 'b');

