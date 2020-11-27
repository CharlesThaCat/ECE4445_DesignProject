QTC = 1/sqrt(2);
QMC = 3.5; % CAN CHANGE THIS!!!!!!!!!!!!
QEC = (QMC*QTC)/(QMC-QTC);
alpha = ((QEC/QES)^2) - 1;
VAB = VAS/alpha;
VAB = VAB;
alpha_1 = VAS/VAB;
fc = fs*sqrt(1+alpha);
VB = VAB/1.2; % CAN CHANGE THE DENOMINATOR!!!!!!!!!!!!
k = nthroot(VB/(0.6*1.0*1.6), 3); % ratio to time box dimensions

CAB = VAB/(1.18*345*345);
CAT = (CAS*CAB)/(CAS+CAB);
% rho = 1.5;  % CAN CHANGE THIS!!!!!!!!!!!!
rho = (8*1.18)/(3*0.65*pi);
MAB = (0.65*rho)/(pi*driverRadius);
MAD = MMD/((pi*driverRadius^2)^2);
MA1 = (8*1.18)/(3*driverRadius*(pi^2)); % NEED TO BE VERIFIED
MAC = MAB + MA1 + MAD;
MAS = MAD + 2*MA1;
RAC = sqrt(MAC/CAT)/QMC;
RAB = RAC - RAS;
RATC = RAE + RAC;

%% Infinite baffle plots
rho_0 = 1.18;
jomega = 2*pi.*frequency*1j;
SD = pi*driverRadius^2;
Bl = sqrt(Bl_square);
% figure;
% A = readmatrix("dp-part-2-ib-ZVC.csv");
% ms_f = A(:, 1);
% ms_ZVC = A(:, 2);
% loglog(ms_f, abs(ms_ZVC), 'r');
% 
% figure;
% A = readmatrix("dp-part-2-ib-UD-noL.csv");
% ms_f = A(:, 1);
% ms_UD = A(:, 2);
% ms_jomega = 2*pi*ms_f;
% ms_p = rho_0 / (2*pi) * ms_jomega .* ms_UD;
% ms_spl = 20.*log10(ms_p./(2*10^(-5)));
% loglog(ms_f, abs(ms_spl), 'g'); hold on;
% 
% A = readmatrix("dp-part-2-ib-UD.csv");
% ms_f = A(:, 1);
% ms_UD = A(:, 2);
% ms_jomega = 2*pi*ms_f;
% ms_p = rho_0 / (2*pi) * ms_jomega .* ms_UD;
% ms_spl = 20.*log10(ms_p./(2*10^(-5)));
% loglog(ms_f, abs(ms_spl), 'b'); hold on; title('infinite baffle');
% 
% % legend('ZVC', 'p no L', 'p');

%% Closed box plots
% figure;
% A = readmatrix("dp-part-2-cb-ZVC.csv");
% ms_f = A(:, 1);
% ms_ZVC = A(:, 2);
% loglog(ms_f, abs(ms_ZVC), 'r');

%% SPL
Gs = (jomega/(2*pi*fc)).^2 ./ ((jomega/(2*pi*fc)).^2 + 1/QTC * jomega/(2*pi*fc) + 1);
cb_SPL_mod_nol = SPL(abs(rho_0/(2*pi) * Bl/(SD*RE*MAC) * Gs)/sqrt(2));

omega_u1 = RE*MAC/(MAD*LE);
f_u1 = omega_u1 / (2*pi);
Ts = 1 ./ (1 + jomega/omega_u1);
cb_SPL_mod = SPL(abs(rho_0/(2*pi) * Bl/(SD*RE*MAC) * Gs .* Ts)/sqrt(2)); % convert to prms first

figure;
A = readmatrix("dp-part-2-cb-UD-noL.csv");
ms_f = A(:, 1);
ms_UD = A(:, 2);
ms_omega = 2*pi*ms_f;
% ms_p = rho_0 / (2*pi) * ms_jomega .* ms_UD;
% ms_spl = 20.*log10(ms_p./(2*10^(-5)));
cb_SPL_sim_nol = SPL(rho_0 / (2*pi) * ms_omega .* ms_UD / sqrt(2)); % convert to prms first
cb_UD_sim_nol = ms_UD;
cb_xD_sim_nol = cb_UD_sim_nol ./ ms_omega;
semilogx(frequency, cb_SPL_mod_nol, 'r--'); hold on;
semilogx(ms_f, cb_SPL_sim_nol, 'b'); hold on;

A = readmatrix("dp-part-2-cb-UD.csv");
ms_f = A(:, 1);
ms_UD = A(:, 2);
ms_omega = 2*pi*ms_f;
% ms_p = rho_0 / (2*pi) * ms_jomega .* ms_UD;
% ms_spl = 20.*log10(ms_p./(2*10^(-5)));
cb_SPL_sim = SPL(rho_0 / (2*pi) * ms_omega .* ms_UD / sqrt(2)); % convert to prms first
cb_UD_sim = ms_UD;
cb_xD_sim = cb_UD_sim ./ ms_omega;
semilogx(frequency, cb_SPL_mod, 'm--');
semilogx(ms_f, cb_SPL_sim, 'c'); ylabel('on-axis SPL (dB)'); xlabel('frequency (Hz)'); 
title('closed-box system'); ylim([0 90])
legend('Modeled SPL without L_E','Simulated SPL without L_E', 'Modeled SPL with L_E', 'Simulated SPL with L_E','Location', 'southeast');

A = readmatrix("dp-part-2-ib-UD.csv");
ms_f = A(:, 1);
ms_UD = A(:, 2);
ms_omega = 2*pi*ms_f;
ib_SPL_sim = SPL(rho_0 / (2*pi) * ms_omega .* ms_UD / sqrt(2)); % convert to prms first
figure;
semilogx(ms_f, ib_SPL_sim, 'r'); hold on;
semilogx(ms_f, cb_SPL_sim, 'b'); ylabel('on-axis SPL (dB)'); xlabel('frequency (Hz)'); 
legend('Simulated infinite baffle SPL','Simulated closed box SPL','Location', 'southeast');
%% different sizes
figure;
A = readmatrix("dp-part-2-cb-UD.csv");
ms_f = A(:, 1);
ms_UD = A(:, 2);
ms_omega = 2*pi*ms_f;
cb_SPL_sim = SPL(rho_0 / (2*pi) * ms_omega .* ms_UD / sqrt(2));
semilogx(ms_f, cb_SPL_sim, 'r'); hold on;

A = readmatrix("dp-part-2-cb-UD-quarter vol.csv");
ms_f = A(:, 1);
ms_UD = A(:, 2);
ms_omega = 2*pi*ms_f;
cb_SPL_sim = SPL(rho_0 / (2*pi) * ms_omega .* ms_UD / sqrt(2));
semilogx(ms_f, cb_SPL_sim, 'g'); title('closed box'); hold on;

A = readmatrix("dp-part-2-cb-UD-four times vol.csv");
ms_f = A(:, 1);
ms_UD = A(:, 2);
ms_omega = 2*pi*ms_f;
cb_SPL_sim = SPL(rho_0 / (2*pi) * ms_omega .* ms_UD / sqrt(2));
semilogx(ms_f, cb_SPL_sim, 'b'); ylabel('on-axis SPL (dB)'); xlabel('frequency (Hz)'); 
title('closed-box system');
legend('Simulated SPL with original box size', 'Simulated SPL with one-quarter box size', 'Simulated SPL with four times box size','Location', 'southeast');

%% UD
cb_jomega_sim = 1j*2*pi*ms_f; % use the same freq points as the simulated curve does
cb_UD_mod = abs(SD/Bl * RAE/RATC * (1/QTC * cb_jomega_sim / (2*pi*fc)) ./ ((cb_jomega_sim/(2*pi*fc)).^2 + 1/QTC*cb_jomega_sim/(2*pi*fc) + 1));
figure;
semilogx(ms_f, cb_UD_mod, 'b');
% semilogx(ms_f, cb_UD_sim_nol, 'r');
hold on;
semilogx(ms_f, cb_UD_sim, 'r');
% max point of modeled curve: (fc, 0.002581)
cb_UD_max_mod = 0.002581;
xline(58.8844); % modeled resonance freq
yline(cb_UD_max_mod/sqrt(2)); % modeled -3db line
% max point of simulated curve: (72.4436, 0.0026279)
cb_UD_max_sim = 0.0026279;
cb_fs_sim = 72.4436;
xline(cb_fs_sim, '--')
yline(cb_UD_max_sim/sqrt(2), '--'); 
% plot cutoffs for the modeled curve
% % lower cutoff: (3.90714947735, 0.00224)
% plot(3.90714947735, 0.00224, 'b*');
% % higher cutoff: (72.4436, 0.00224)
% plot(72.4436, 0.00224, 'b*');
legend('Modeled U_D','Simulated U_D');
xlabel('frequency (Hz)');
ylabel('U_D (m^3/s)'); title('closed-box system');


%% xD
figure;
% semilogx(ms_f, cb_xD_sim_nol, 'r');
% hold on;
semilogx(ms_f, cb_xD_sim, 'b');
% legend('Simulated x_D without L_E', 'Simulated x_D with L_E');
xlabel('frequency (Hz)');
ylabel('x_D (m)'); title('closed-box system');