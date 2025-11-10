clear
clc
close all
%% parameters
a = 6371 * 1e3;
M = 5.965e24;
G = 6.672e-11;
D0 = 0;
D = a - D0;
%% true model
load ../data/Earth2014_truncate_360degree.mat
dlat = 0.5; dlon = dlat;
nmax = 360;
[n, m] = creat_nm(nmax);
sigmalog10 = @(x) log10( sqrt(degVar([n, m, real(x), imag(x)], nmax)) );
lon = -180 + dlon / 2 : dlat : 180;
lat = 90 - dlat / 2 : -dlat : -90;
[Lon, Lat] = meshgrid(lon, lat);
H = shs_grid(CS, lon, lat, nmax, 0);
ind = H > 0;
H(ind) = 0;
rho = 630 - 5.6 * 1e-2 * (-H) + 2.8 * 1e-8 * (-H) .^ 2;
MIN = (5.6 * 1e-2 - sqrt( (5.6 * 1e-2) ^ 2 - 4 * 630 * 2.8 * 1e-8 )) / -(2 * 2.8 * 1e-8);
%% plot
cmp1 = load('cpt\batlowW.mat'); cmp1 = cmp1.batlowW;
cmp2 = load('cpt\roma.mat'); cmp2 = cmp2.roma;
cmp1(1, :) = [0.6, 0.6, 0.6];
cmp2(1, :) = [0.6, 0.6, 0.6];
Hdraw = H;
Hdraw(ind) = nan;
rhodraw = rho;
rhodraw(ind) = nan;
fontsize = 9;
fig = figure;
tiledlayout(2, 1, "TileSpacing", "tight", "Padding", "compact", Parent=fig);
nexttile
bm = axesm("MapProjection", "robinson", "MeridianLabel", "on", "ParallelLabel", "on");
setm(bm, FLineWidth=1, FontSize=fontsize, MLabelLocation=90, Frame='on', Grid='on', LabelFormat='compass', PLabelLocation=45, MLabelParallel=-10, MLineLocation=90, PLineLocation=45);
geoshow(Lat, Lon, Hdraw / 1e3, 'DisplayType','texturemap')
set(gca, 'FontSize', fontsize, 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
tightmap
title('(a)', 'fontsize', fontsize);
cba = colorbar;
title(cba, '(km)', 'fontsize', fontsize);
cba.TickDirection = "out";
cba.Units = "centimeters";
cba.Position(1) = 12;
cba.Position(4) = 3.5;
cba.TickLength = 0.02;
clim([-10, 0]);
colormap(cmp1);
axd = nexttile;
bm = axesm("MapProjection", "robinson", "MeridianLabel", "on", "ParallelLabel", "on");
setm(bm, FLineWidth=1, FontSize=fontsize, MLabelLocation=90, Frame='on', Grid='on', LabelFormat='compass', PLabelLocation=45, MLabelParallel=-10, MLineLocation=90, PLineLocation=45);
geoshow(Lat, Lon, rhodraw, 'DisplayType','texturemap')
set(gca, 'FontSize', fontsize, 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
tightmap
title('(b)', 'fontsize', fontsize);
cbb = colorbar;
title(cbb, '(kg m^{-3})', 'fontsize', fontsize);
cbb.TickDirection = "out";
cbb.Units = "centimeters";
cbb.Position(1) = 12;
cbb.Position(4) = 3.5;
cbb.TickLength = 0.02;
clim([0, 600]);
cmb = colormap(axd, cmp2);
set(gcf,'Units','centimeters','Position',[1.3, 3, 14, 10]);
% exportgraphics(gcf, '../manuscript/CG/figures/seafloor_true_model.pdf', Resolution=350);
%% Forward
funRho = @(r, i, j) 630 - 5.6 * 1e-2 * r + 2.8 * 1e-8 * r .^ 2;
tol = 1e-8;
tolH = 1e-4;
gamma = find_gamma(H, D, nmax, tol);
Nq = find_Nq(H, funRho, gamma, D, tol);

Qnm = H2Qnm(lon, lat, H, funRho, D, nmax, "a", a, "M", M, "Nq", Nq, "szm", [dlon, dlat], "gamma", gamma);
rng(1)
noise = randn(size(Qnm)) * 0.05 * mean(abs(Qnm));
QnmAddNoise = Qnm + noise;
QnmAddNoise(1) = Qnm(1);
%% plot
fontsize = 9;
linewidth = 1.2;
fig = figure;
tiledlayout(1, 1, "Padding", "compact", 'TileSpacing', 'tight', Parent=fig);
nexttile
plot(sigmalog10(Qnm), linewidth=linewidth)
hold on
plot(sigmalog10(noise), linewidth=linewidth)
plot(sigmalog10(QnmAddNoise), linewidth=linewidth)
xlim([0, 360]);
xlabel('Spherical harmonic degree');
ylabel('Degree standard deviation in log10');
grid on
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true, 'XTick', 0:60:360);
leg = legend('Synthesized signal', 'Noise', 'Signal plus noise', 'fontsize', fontsize, 'Location', 'best');
set(gcf,'Units','centimeters','Position', [1.3 3 10, 8]);
% exportgraphics(gcf, '../manuscript/CG/figures/seafloor_Qnm.pdf');
%% inversion
[Hinv, iterCorrection] = Qnm2H(lon, lat, QnmAddNoise, funRho, D, nmax, "sigma_nm", 0.05 * mean(abs(Qnm)), "a", a, "M", M, "gamma", gamma, "HRange", [MIN, D0], "szm", [dlon, dlat], "Nq", Nq, "tolH", tolH, "nc", 180);
STD1 = sim_sta(H(~ind), Hinv(~ind));
%% Comparison
cmp2 = load('cpt\vik.mat'); cmp2 = cmp2.vik;
cmp2(1, :) = [0.6, 0.6, 0.6];

fig = figure;
tiledlayout(1, 1, "Padding", "compact", 'TileSpacing', 'tight', Parent=fig);
nexttile
plot(1:7, iterCorrection / 1e3, linewidth, 1.5);
ylim([0, 0.15]);
xlabel('Iterative step', 'fontsize', fontsize);
ylabel('Correction (km)', 'fontsize', fontsize);
grid on
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true, 'XTick', 1:7);
title('(a)', 'fontsize', fontsize);
set(gcf,'Units','centimeters','Position',[1.3, 3, 10, 6]);
exportgraphics(gcf, '../manuscript/CG/figures/seafloor_inv_model_a.pdf');

Hdraw = Hinv;
Hdraw(ind) = nan;

fig = figure;
td = tiledlayout(2, 1, "TileSpacing", "tight", "Padding", "compact", Parent=fig);
nexttile
bm = axesm("MapProjection", "robinson", "MeridianLabel", "on", "ParallelLabel", "on");
setm(bm, FLineWidth=1, FontSize=fontsize, MLabelLocation=90, Frame='on', Grid='on', LabelFormat='compass', PLabelLocation=45, MLabelParallel=-10, MLineLocation=90, PLineLocation=45);
geoshow(Lat, Lon, Hdraw / 1e3, 'DisplayType','texturemap')
set(gca, 'FontSize', fontsize, 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
tightmap
title('(b)', 'fontsize', fontsize);
cba = colorbar;
title(cba, '(km)', 'fontsize', fontsize);
cba.TickDirection = "out";
cba.Units = "centimeters";
cba.Position(1) = 12;
cba.Position(4) = 3.5;
cba.TickLength = 0.02;
clim([-10, 0]);
colormap(cmp1);

axd = nexttile;
bm = axesm("MapProjection", "robinson", "MeridianLabel", "on", "ParallelLabel", "on");
setm(bm, FLineWidth=1, FontSize=fontsize, MLabelLocation=90, Frame='on', Grid='on', LabelFormat='compass', PLabelLocation=45, MLabelParallel=-10, MLineLocation=90, PLineLocation=45);
geoshow(Lat, Lon, (H - Hdraw) / 1e3, 'DisplayType','texturemap')
set(gca, 'FontSize', fontsize, 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
tightmap
title('(c)', 'fontsize', fontsize);
cbb = colorbar;
title(cbb, '(km)', 'fontsize', fontsize);
cbb.TickDirection = "out";
cbb.Units = "centimeters";
cbb.Position(1) = 12;
cbb.Position(4) = 3.5;
cbb.TickLength = 0.02;
colormap(axd, cmp2);
clim([-3, 3]);
set(gcf,'Units','centimeters','Position',[1.3, 3, 14, 10]);
% exportgraphics(gcf, '../manuscript/CG/figures/seafloor_inv_model_bc.pdf', Resolution=350);