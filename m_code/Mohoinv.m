clear
clc
close all
%% parameters
nmax = 300;
a = 6371 * 1e3;
M = 5.965e24;
G = 6.672e-11;
D0 = 30 * 1e3;
area_ex = [60, 150, -5, 65];
area = [65, 145, 0, 60];
ex = 5;
dlon = 1; dlat = dlon;
hdclp = @(x) gridclip(x, dlon, dlat, ex);
[n, m] = creat_nm(nmax);
sigmalog10 = @(x) log10( degVar([n, m, real(x), imag(x)], nmax) );
%% true model
load ../data/ECM1_crust_60E_150E_5S_65N.mat
Lon = rho.Lon; Lat = rho.Lat;
[row, col] = size(Lon);
D = rho.crs.dp3;
lon = Lon(1, :);
lat = Lat(:, 1);

MIN = -70 * 1e3;
rhom = 3300;

rho0 = zeros(row, col) + 2500;
rho1 = rho.crs.rho1 * 1e3;
rho2 = rho.crs.rho2 * 1e3;
rho3 = rho.crs.rho3 * 1e3;
rho4 = zeros(row, col) + 3300;

dp0 = zeros(row, col) + D0;
dp1 = -rho.crs.dp1 * 1e3 + D0;
dp2 = -rho.crs.dp2 * 1e3 + D0;
dp3 = -rho.crs.dp3 * 1e3 + D0;
dp4 = zeros(row, col) + MIN;

dp = [dp0(:), dp1(:), dp2(:), dp3(:), dp4(:)];
drho = rhom - [rho0(:), rho1(:), rho2(:), rho3(:), rho4(:)];

rho0 = zeros(row * col, 4);
rho1 = rho0;
% piecewise linear fitting
for j = 1 : 4
    for i = 1 : row * col
        s = polyfit(dp(i, [j, j + 1]), drho(i, [j, j + 1]), 1);
        rho0(i, j) = s(2);
        rho1(i, j) = s(1);
    end
end

rho0 = reshape(rho0, row, col, []);
rho1 = reshape(rho1, row, col, []);
dp = reshape(dp, row, col, []);
funRho = @(x, i, j) funRho_piecelinear(x, rho0, rho1, dp, i, j);
shs = @(x) shs_grid([n, m, real(x), imag(x)], lon, lat, nmax, 0, "a", a, "GM", G * M, "height", 0, "fldTyp", "vz");
H = D0 - D * 1e3;
D = a - D0;
%% plot
figure;
fontsize = 9;
cmp = load("cpt\batlowW.mat"); cmp = cmp.batlowW;
bm = axesm("MapProjection", "eqdconic", "MapLonLimit", area(1:2), "MapLatLimit", area(3:4), "MeridianLabel", "on", "ParallelLabel", "on");
setm(bm, FLineWidth=1, FontSize=fontsize, MLabelLocation=65:20:145, MLineLocation=65:20:145, Frame='on', Grid='on', LabelFormat='compass', MLabelParallel='south', PLabelLocation=0:20:60, PLineLocation=0:20:60);
geoshow(Lat, Lon, (D0 - H) / 1e3, 'DisplayType', 'texturemap')
set(gca, 'FontSize', fontsize);
cba = colorbar;
title(cba, '(km)', "Fontsize", fontsize);
cba.Units = "centimeters";
cbapos = cba.Position;
cbapos(4) = 5;
cbapos(3) = 0.4;
cbapos(1) = 10.3;
cbapos(2) = 2.5;
cba.Position = cbapos;
cba.TickDirection = "out";
cba.TickLength = 0.02;
clim([0, 80]);
colormap(cmp);
tightmap
set(gca,'Units','centimeters','Position', [1 0.2 9, 9], 'Color', 'none', 'XColor','none','YColor','none');
set(gcf,'Units','centimeters','Position', [1.3 3 12, 9]);
% exportgraphics(gcf, '../manuscript/CG/figures/Moho_depth.pdf', Resolution=350);
%% plot
fig = figure;
tiledlayout(4, 2, "Padding", "loose", 'TileSpacing', 'compact', Parent=fig);
cmp_density = load("cpt\roma.mat"); cmp_density = cmp_density.roma;
for i = 1 : 4
    nexttile
    bm = axesm("MapProjection", "eqdconic", "MapLonLimit", area(1:2), "MapLatLimit", area(3:4), "ParallelLabel", "on");
    setm(bm, FLineWidth=1, FontSize=fontsize, Frame='on', Grid='on', MLineLocation=65:20:145, MLabelParallel='south', PLabelLocation=0:20:60, PLineLocation=0:20:60);
    if i == 4
        setm(bm, MLabelLocation=65:20:145, LabelFormat='compass', MLabelParallel='south', MeridianLabel= "on");
    end
    if i == 3
        cb = colorbar;
        title(cb, '(kg m^{-3})');
        cb.Units = "centimeters";
        cbpos = cb.Position;
        cb.Position = [cbpos(1), cbpos(2)+1, 0.3, 4];
        cb.Ticks = -800:400:800;
        cb.TickDirection = "out";
        cb.TickLength = 0.02;
    end
    geoshow(Lat, Lon, rho0(:, :, i), 'DisplayType','texturemap')
    set(gca, 'Color','none','XColor','none','YColor','none');
    clim([-800, 800])
    colormap(cmp_density)
    tightmap
    title(['(',char(97 + i - 1),')'], 'FontSize', fontsize);

    nexttile
    bm = axesm("MapProjection", "eqdconic", "MapLonLimit", area(1:2), "MapLatLimit", area(3:4), "ParallelLabel", "on");
    setm(bm, FLineWidth=1, FontSize=fontsize, Frame='on', Grid='on', MLabelParallel='south', MLineLocation=65:20:145, PLabelLocation=0:20:60, PLineLocation=0:20:60);
    if i == 4
        setm(bm, MLabelLocation=65:20:145, LabelFormat='compass', MLabelParallel='south', MeridianLabel= "on");
    end
    if i == 3
        cb = colorbar;
        title(cb, '(kg m^{-4})');
        cb.Units = "centimeters";
        cbpos = cb.Position;
        cb.Position = [cbpos(1) - 0.2, cbpos(2)+1, 0.3, 4];
        cb.TickDirection = "out";
        cb.TickLength = 0.02;
    end
    clim([0, 0.04]);
    geoshow(Lat, Lon, rho1(:, :, i), 'DisplayType','texturemap')
    set(gca, 'FontSize', fontsize, 'Color','none','XColor','none','YColor','none');
    tightmap
    title(['(',char(101 + i - 1),')'], 'FontSize', fontsize);
end
set(gcf,'Units','centimeters','Position', [1.3 3 14, 14]);
% exportgraphics(gcf, '../manuscript/CG/figures/Moho_linear_density.pdf', Resolution=350);
%% Forward
tol = 1e-8;
tolH = 1e-4;
gamma = find_gamma(H, D, nmax, tol);
Nq = find_Nq(H, funRho, gamma, D, tol);
Fwd = @(xH, xfunRho, xD) H2Qnm(lon, lat, xH, xfunRho, xD, nmax, "a", a, "M", M, "Nq", Nq, "szm", [dlon, dlat], "gamma", gamma);
Qnm = Fwd(H, funRho, D);
% noise
sigmaNoise = rd_gfc('../data/GOCO06s.gfc');
sigmaNoise = sigmaNoise.CS(:, 5:6);
rng(1)
QnmNoise = randn(size(Qnm)) .* sigmaNoise;
QnmNoise = complex(QnmNoise(:, 1), QnmNoise(:, 2));
QnmAddNoise = Qnm + QnmNoise;

sigma_nm = 10 .^ sigmalog10(QnmNoise);
sigma_nm = cs2cnm(repmat(sigma_nm, 1, nmax + 1), nmax);
sigma_nm = sigma_nm(:, 3);
%% plot
fontsize = 9;
linewidth = 1.2;
fig = figure;
tiledlayout(1, 1, "Padding", "compact", 'TileSpacing', 'tight', Parent=fig);
nexttile
plot(sigmalog10(Qnm), linewidth=linewidth)
hold on
plot(sigmalog10(QnmNoise), linewidth=linewidth)
plot(sigmalog10(QnmAddNoise), linewidth=linewidth)
xlim([0, 300]);
xlabel('Spherical harmonic degree');
ylabel('Degree standard deviation in log10');
grid on
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true, 'XTick', 0:60:360);
leg = legend('Synthesized signal', 'Noise', 'Signal plus noise', 'fontsize', fontsize, 'Location', 'Northeast');
set(gcf,'Units','centimeters','Position', [1.3 3 10, 8]);
% exportgraphics(gcf, '../manuscript/CG/figures/Moho_Qnm.pdf');
%% cp
cpd = hdclp(D0 + H);
indcp = randperm(numel(cpd), 20);
cpd = cpd(indcp);
cplon = hdclp(Lon); cplon = cplon(indcp);
cplat = hdclp(Lat); cplat = cplat(indcp);
%% estimation D
vz = hdclp(shs(QnmAddNoise));
vzcp = vz(indcp) * -1e5;
s = polyfit(vzcp, cpd, 1);
D_ = a - s(2);
%% plot
fig = figure;
tiledlayout(2, 1, "Padding", "tight", 'TileSpacing', 'loose', Parent=fig);
nexttile
bm = axesm("MapProjection", "eqdconic", "MapLonLimit", area(1:2), "MapLatLimit", area(3:4), "MeridianLabel", "on", "ParallelLabel", "on");
setm(bm, FLineWidth=1, FontSize=fontsize, MLabelLocation=65:20:145, MLineLocation=65:20:145, Frame='on', Grid='on', LabelFormat='compass', MLabelParallel='south', PLabelLocation=0:20:60, PLineLocation=0:20:60);
scatterm(cplat, cplon, 15, cpd / 1e3, "filled");
cba = colorbar;
title(cba, '(km)', "Fontsize", fontsize);
cba.Units = "centimeters";
cbapos = cba.Position;
cba.Position = [cbapos(1)-1.9, cbapos(2)+3, cbapos(3)-0.1, cbapos(4)-2];
cba.TickDirection = "out";
cba.TickLength = 0.02;
clim([0, 80]);
colormap(cmp);
tightmap
set(gca, 'Color','none','XColor','none','YColor','none');
title('(a)', FontSize=fontsize)
nexttile
plot(vzcp, s(1) * 1e-3 * vzcp + s(2) * 1e-3, linewidth=linewidth, Color=[0.8500 0.3250 0.0980]);
hold on 
scatter(vzcp, cpd / 1e3, 15, [0 0.4470 0.7410], 'filled');
box on
grid on
title('(b)', FontSize=fontsize)
xlabel('Gravitational acceleration (mGal)');
ylabel('Moho depth (km)');
ylim([0, 60]);
text(-500, 40, ['$y=', num2str(s(1) * 1e-3), 'x+', 'a-', num2str((a - s(2)) * 1e-3),'$'], 'Interpreter', 'latex');
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);
set(gcf,'Units','centimeters','Position', [1.3, 3, 10, 14]);
% exportgraphics(gcf, '../manuscript/CG/figures/Moho_estimate_D.pdf', Resolution=350);
%% inversion
HRange = [MIN + D - D_, a - D_];
Inv = @(xQnm) Qnm2H(lon, lat, xQnm, funRho, D_, nmax, "sigma_nm", sigma_nm, "a", a, "M", M, "gamma", gamma, "HRange", HRange, "szm", [dlon, dlat], "Nq", Nq, "tolH", tolH, "extend", ex, "nc", 200);
% comparison
[Hinv, iterInfo] = Inv(QnmAddNoise);
STD1 = sim_sta(hdclp(H), hdclp(Hinv));
%% plot
cmp2 = load("cpt\vik.mat"); cmp2 = cmp2.vik;
fig = figure;
tiledlayout(2, 2, "Padding", "loose", 'TileSpacing', 'compact', Parent=fig);
nexttile([1, 2])
plot(iterInfo / 1e3, 'LineWidth', linewidth);
ylim([0, 3])
yticks(0:1:3)
title('(a)', 'FontSize', fontsize);
ylabel('Correction (km)');
xlabel('Iteration step');
grid on
set(gca, 'FontSize', fontsize, 'GridLineStyle', '--', 'XMinorTick', true, 'YMinorTick', true);

ax1 = nexttile;
bm = axesm("MapProjection", "eqdconic", "MapLonLimit", area(1:2), "MapLatLimit", area(3:4), "ParallelLabel", "on");
setm(bm, FLineWidth=1, FontSize=fontsize, Frame='on', Grid='on', MLineLocation=65:20:145, MLabelLocation=65:20:145, MLabelParallel='south', MeridianLabel="on", LabelFormat='compass', PLabelLocation=0:20:60, PLineLocation=0:20:60);
cb = colorbar;
cb.Units = "centimeters";
cb.Label.String = '(km)';
cb.TickDirection = "out";
cb.Location = "southoutside";
cbpos = cb.Position;
cb.Position = [cbpos(1) + 0.5, cbpos(2) - 0.6, 4, 0.3];
cb.TickLength = 0.02;
geoshow(Lat, Lon, (a - (D_ + Hinv)) / 1e3, 'DisplayType','texturemap')
set(gca, 'Color','none','XColor','none','YColor','none');
clim([0, 80])
colormap(ax1, cmp)
tightmap
title('(b)', 'FontSize', fontsize);
ax2 = nexttile;
bm = axesm("MapProjection", "eqdconic", "MapLonLimit", area(1:2), "MapLatLimit", area(3:4), "ParallelLabel", "on");
setm(bm, FLineWidth=1, FontSize=fontsize, Frame='on', Grid='on', MLineLocation=65:20:145, MLabelLocation=65:20:145, MLabelParallel='south', MeridianLabel="on", LabelFormat='compass', ParallelLabel='off');
cb = colorbar;
cb.Units = "centimeters";
cb.Label.String = '(km)';
cb.TickDirection = "out";
cb.Location = "southoutside";
cbpos = cb.Position;
cb.Position = [cbpos(1) + 0.65, cbpos(2) - 0.6, 4, 0.3];
cb.TickLength = 0.02;
geoshow(Lat, Lon, (H - Hinv) / 1e3, 'DisplayType','texturemap')
set(gca, 'Color','none','XColor','none','YColor','none');
clim([-2, 2])
colormap(ax2, cmp2)
tightmap
title('(c)', 'FontSize', fontsize);

set(gcf,'Units','centimeters','Position', [1.3, 3, 15, 14]);
% exportgraphics(gcf, '../manuscript/CG/figures/Moho_inv.pdf', Resolution=350);

