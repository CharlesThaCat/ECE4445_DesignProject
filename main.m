clear; clc;

%% off box test
path = '09_offbox.txt';
[driverRadius,RE,VT,Vgen,frequency,magnitude,phase] = txtParser(path);
% figure; loglog(frequency, magnitude);

[M,ind] = max(magnitude, [], 'linear');
RES = M - RE;
fs = frequency(ind);
R1 = sqrt(RE*(RE+RES));
p = polyfit(frequency(1:ind), magnitude(1:ind), 14); % 14th or 15th order both good
x1 = linspace(10, fs, 20000);
y1 = polyval(p,x1);
% figure; loglog(frequency, magnitude); hold on;
% loglog(x1,y1); hold off;  xlabel('frequency (Hz)'); ylabel('Z_{VC} magnitude (ohms)');
% legend('measured data', 'fitted polynomial', 'Location', 'southeast')
[~,ix1] = min(abs(y1-R1));
f1 = x1(ix1);
p = polyfit(frequency(ind:45), magnitude(ind:45), 15); % 14th or 15th order both good
x1 = linspace(fs, 27, 20000);
y1 = polyval(p,x1);
% figure; loglog(frequency, magnitude); hold on;
% loglog(x1,y1); hold off;  xlabel('frequency (Hz)'); ylabel('Z_{VC} magnitude (ohms)');
% legend('measured data', 'fitted polynomial', 'Location', 'southeast')
[~,ix1] = min(abs(y1-R1));
f2 = x1(ix1);
QMS = (fs / (f2-f1))*sqrt((RE+RES) / RE);
QES = QMS*RE/RES;
QTS = RE*QMS/(RE+RES);
% The last paragraph of textbook 11.7.3 says use R1 to determine both f1 and f2, then use f1 and f2
% to calculate and compare to fs as a check. The percentage difference between these should be smaller
% than 5%, which is satisfied here.
clear M ind p x1 y1 ix1

%% on box test
path = '09_onbox.txt';
[driverRadius,RE,VT,Vgen,frequencyCT,magnitudeCT,phaseCT] = txtParser(path);
% figure; loglog(frequencyCT, magnitudeCT);
[M,ind] = max(magnitudeCT(1:93), [], 'linear');
RES = M - RE;
fCT = frequencyCT(ind);
R1 = sqrt(RE*(RE+RES));
p = polyfit(frequencyCT(1:ind), magnitudeCT(1:ind), 18); % 18th or 19th order both good
x1 = linspace(10, fCT, 20000);
y1 = polyval(p,x1);
% figure; loglog(frequencyCT, magnitudeCT); hold on;
% loglog(x1,y1); hold off; grid on;
[~,ix1] = min(abs(y1-R1));
f1CT = x1(ix1);
p = polyfit(frequencyCT(ind:92), magnitudeCT(ind:92), 16); % 16th or 17th order both good
x1 = linspace(fCT, 100, 20000);
y1 = polyval(p,x1);
% figure; loglog(frequencyCT, magnitudeCT); hold on;
% loglog(x1,y1); hold off; grid on;
[~,ix1] = min(abs(y1-R1));
f2CT = x1(ix1);
QMCT = (fCT / (f2CT-f1CT))*sqrt((RE+RES) / RE);
QECT = QMCT*RE/RES;
VAS = VT*(-1 + ( (fCT*QECT)/(fs*QES) ));
clear M ind p x1 y1 ix1

%% n, Le, LE
ind1 = 250; % choose different points to get closest approx. (219)
ind2 = 301; % fix at highest frequency
% ind1 = 1:1:301;
% ind2 = 301;
w1 = 2.*pi.*frequency(ind1);
w2 = 2.*pi.*frequency(ind2);
Z1 = magnitude(ind1).*exp(1j.*deg2rad(phase(ind1))); Y1 = 1./Z1;
Z2 = magnitude(ind2).*exp(1j.*deg2rad(phase(ind2))); Y2 = 1./Z2;
n = log10( real(Y1)./real(Y2) ) ./ log10( w2./w1 );
Le = cos(n.*pi./2) ./ (real(Y1).*(w1.^n));
LE = -1 ./ (w1 .* (imag(Y1) + (sin(n.*pi./2) ./ (Le.*(w1.^n))) ));
clear ind1 ind2 w1 w2 Z1 Z2 Y1 Y2

%% Small tweaking
% LE = 0.012;
% Le = 0.11;
% n = 0.665;

Le = 0.1300;
LE = 0.0100;
n = 0.6510;

%% replot ZVC (Zmot, ZEL)
ZL1 = 2*pi.*frequency*1j*LE; ZL2 = ((2*pi.*frequency*1j).^n).*Le;
ZEL = ZL1.*ZL2./(ZL1+ZL2);
RES = RE * QMS / QES;
Zmot = RES * (1/QMS)*(frequency*1j/fs) ./ ((frequency*1j/fs).*(frequency*1j/fs) + (1/QMS)*(frequency*1j/fs) + 1);

% Zmot = Bl_square ./ (RMS + 1j.*w.*MMS + (1./(1j.*w.*CMS)));
ZVC = RE + Zmot + ZEL;
figure; loglog(frequency, magnitude); hold on;
loglog(frequency, abs(ZVC)); hold off; 
xlabel('frequency (Hz)'); ylabel('Z_{VC} magnitude (ohms)');
legend('From measured data', 'From data modeled parameters', 'Location', 'southeast');
figure; semilogx(frequency, deg2rad(phase)); hold on;
semilogx(frequency, angle(ZVC)); hold off; title('Z_{VC} phase');
legend('From original data', 'From data modeled parameters', 'Location', 'southeast');

figure; loglog(frequency, magnitude); hold on;
loglog(frequency, abs(ZL1)); hold on;
loglog(frequency, abs(ZL2)); 
loglog(frequency, abs(ZEL)); hold off; 
xlabel('frequency (Hz)'); title('magnitude of impedances of inductances');
legend('Z_{VC} from measured data', 'Z_{L_E} from data modeled parameters', 'Z_{L_e} from data modeled parameters', 'Parallel combination of Z_{L_E} and Z_{L_e}', 'Location', 'southeast');
% figure; loglog(frequency, deg2rad(phase)); hold on;
% loglog(frequency, angle(ZEL)); hold off; title('Z_{EL} phase');
% legend('From original data', 'From data modeled parameters', 'Location', 'southeast');

%% Conversion to infinite-baffle parameters
% Seems like it's better not doing this step. We can compare the results of ZVC to illustrate.
kM = sqrt( 1 + 0.2699*( (VAS*(2*pi*fs)^2) / (driverRadius*345*345) ) );
fs = fs/kM;
QMS = QMS*kM; QES = QES*kM; QTS = QTS*kM;

%% Three part model parameters (not verified yet)
CMS = VAS/(1.18*(345^2)*(pi*(driverRadius^2))^2);
CAS = VAS/(1.18*(345^2));
MAS = 1/(4*pi*pi*fs*fs*CAS);
MMS = 1/(4*pi*pi*fs*fs*CMS);
RMS = (sqrt(MMS/CMS)) / QMS;
RAS = (sqrt(MAS/CAS)) / QMS;
RAE = -RAS + ((sqrt(MAS/CAS))/QTS);
Bl_square = (RE*sqrt(MMS/CMS)) / QES;
MMD = ((pi*(driverRadius^2))^2) * (MAS-2*( (8*1.18) / (3*pi*pi*driverRadius) ));
