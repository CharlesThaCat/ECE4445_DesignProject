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
% loglog(x1,y1); hold off; grid on;
[~,ix1] = min(abs(y1-R1));
f1 = x1(ix1);
p = polyfit(frequency(ind:45), magnitude(ind:45), 15); % 14th or 15th order both good
x1 = linspace(fs, 27, 20000);
y1 = polyval(p,x1);
% figure; loglog(frequency, magnitude); hold on;
% loglog(x1,y1); hold off; grid on;
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