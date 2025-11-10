clear
close all
clc
% parameter
a = 6371 * 1e3;
D0 = 30 * 1e3;
D = a - D0;
H = 10 * 1e3;
M = 5.965e24;
G = 6.672e-11;
nmax = 360;
[n, m] = creat_nm(nmax);
sigmalog10 = @(x) log10( sqrt(degVar([n, m, real(x), imag(x)], nmax)) );
%% three density functions
rho_poly = [100, 1.0e-3, 1.9e-6, -1.0e-10, 4.0e-14];
rho_exp = [100, 1.79 * 1e-4];
funRhoPoly = @(xr, i, j) rho_poly(1) + rho_poly(2) * xr + rho_poly(3) * xr .^ 2 + rho_poly(4) * xr .^ 3 + rho_poly(5) * xr .^ 4;
funRhoExp = @(xr, i, j) rho_exp(1) * exp(rho_exp(2) * xr);
Qnm_ref_poly = sph_shell_Q00_poly(H, D, rho_poly, M);
Qnm_ref_exp = sph_shell_Q00_exp(H, D, rho_exp, M);
%% plot
HRange = [-D0, D0];
x = linspace(HRange(1), HRange(2), 100);
fontsize = 9;
linewidth = 1.2;
fig = figure;
tiledlayout(1, 1, "Padding", "compact", 'TileSpacing', 'tight', Parent=fig);
nexttile
plot(x * 1e-3, funRhoPoly(x), x * 1e-3, funRhoExp(x), "LineWidth", linewidth);
xlim([0, H] * 1e-3);
xlabel('Topography $H$ (km)', 'Interpreter', 'latex');
ylabel('Radial density (kg m$^{-3}$)', 'Interpreter', 'latex');
grid on
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);
leg = legend('Polynomial', 'Exponential', 'fontsize', fontsize, 'Location', 'northwest');
set(gcf,'Units','centimeters','Position', [1.3 3 8, 6]);
% exportgraphics(gcf, '../manuscript/CG/figures/shell_density_profiles.pdf');
%% dis tesseroid
dlat = 0.5; dlon = dlat;
lon = -180 + dlon / 2 : dlat : 180;
lat = 90 - dlat / 2 : -dlat : -90;
[Lon, Lat] = meshgrid(lon, lat);
tol = 1e-8;
tolH = 1e-4;
%% find best Nq and gamma
gamma = find_gamma(H, D, nmax, tol);
Nq_exp = find_Nq(H, funRhoExp, gamma, D, tol);
Nq_poly = find_Nq(H, funRhoPoly, gamma, D, tol);
%% forward
Qnm_cal_exp = H2Qnm(lon, lat, H, funRhoExp, D, nmax, "a", a, "M", M, "Nq", Nq_exp, "szm", [dlon, dlat], "gamma", gamma);
Qnm_cal_poly = H2Qnm(lon, lat, H, funRhoPoly, D, nmax, "a", a, "M", M, "Nq", Nq_poly, "szm", [dlon, dlat], "gamma", gamma);

dQnm_exp = Qnm_cal_exp;
dQnm_exp(1) = dQnm_exp(1) - Qnm_ref_exp;
dQnm_poly = Qnm_cal_poly;
dQnm_poly(1) = dQnm_poly(1) - Qnm_ref_poly;
%% plot
fig = figure;
td = tiledlayout(2, 1, "TileSpacing", "tight", "Padding", "tight", Parent=fig);
nexttile
plot(0 : nmax, sigmalog10(dQnm_poly), 'LineWidth', 0.5);
xlim([0, 360])
ylim([-26, -14])
grid on
title('(a)', 'fontsize', fontsize);
set(gca, 'XTick', 0:60:360, 'YTick', -26:4:-14);
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);

nexttile
plot(0 : nmax, sigmalog10(dQnm_exp), 'LineWidth', 0.5);
xlim([0, 360])
ylim([-26, -14])
grid on
title('(b)', 'fontsize', fontsize);
set(gca, 'XTick', 0:60:360, 'YTick', -26:4:-14);
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);
set(gcf,'Units','centimeters','Position',[1.3, 3, 10, 8]);

xlabel(td, 'Spherical harmonic degree', 'fontsize', fontsize);
ylabel(td, 'Degree standard deviation in log10', 'fontsize', fontsize);
% exportgraphics(gcf, '../manuscript/CG/figures/shell_degVar_dQnm.pdf');
%% inverse
r_nm = inf(size(Qnm_cal_poly));
r_nm(1) = 0;

[HinvExp, iterCorrectionExp, iterRMSExp] = Qnm2H(lon, lat, Qnm_cal_exp, funRhoExp, D, nmax, "r_nm", r_nm, "mu", 1, "a", a, "M", M, "gamma", gamma, "HRange", HRange, "szm", [dlon, dlat], "Nq", Nq_exp, "tolH", tolH, Hcheck=H);
[HinvPoly, iterCorrectionPoly, iterRMSPoly] = Qnm2H(lon, lat, Qnm_cal_poly, funRhoPoly, D, nmax, "r_nm", r_nm, "mu", 1, "a", a, "M", M, "gamma", gamma, "HRange", HRange, "szm", [dlon, dlat], "Nq", Nq_poly, "tolH", tolH, Hcheck=H);
%% plot
fontsize = 9;
fig = figure;
tiledlayout(1, 2, "TileSpacing", "tight", "Padding", "tight", Parent=fig);
nexttile
plot(1:3, iterCorrectionPoly / 1e3, linewidth=1.2);
ylabel('Correction (km)', 'fontsize', fontsize);
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true, 'XTick', 1:3, 'YTick', (0:2:8) * 1e-3);
yyaxis right
plot(1:3, iterRMSPoly / 1e3, linewidth=1.2);
ylim([1.2, 1.5]*1e-4);
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true, 'XTick', 1:3, 'YTick', (1.2:0.1:1.5) * 1e-4);
xlabel('Iterative step', 'fontsize', fontsize);
title('(a)', 'fontsize', fontsize);

nexttile
plot(1:3, iterCorrectionExp / 1e3, linewidth=1.2);
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true, 'XTick', 1:3, 'YTick', (0:2:10) * 1e-3);
axb = gca;
axb.YAxis.Exponent = -3;
yyaxis right
plot(1:3, iterRMSExp / 1e3, linewidth=1.2);
ylim([1.3, 1.6]*1e-5);
ylabel('Misfit (km)', 'fontsize', fontsize);
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true, 'XTick', 1:3, 'YTick', (1.3:0.1:1.6) * 1e-5);
xlabel('Iterative step', 'fontsize', fontsize);
title('(b)', 'fontsize', fontsize);

set(gcf,'Units','centimeters','Position',[1.3, 3, 12, 8]);
% exportgraphics(gcf, '../manuscript/CG/figures/shell_dH.pdf', 'ContentType', 'vector');
