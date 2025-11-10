function [Hk, iterCorrection, iterRMS] = Qnm2H(lon, lat, Qnm, funRho, D, nmax, opt)
%{
    Invert potential (multipoles) coefficients to interface relief H.
Input:
    lon: the longitude vector of the center of tesseroids [deg]
    lat: the latitude vector of the center of tesseroids [deg]
    Qnm: [ (nmax+1)(nmax+2)/2 × 1 complex] multipoles due to interface
    relief
    funRho: function handle, funRho(r,i,~) -> density at (r + D) for row i
           (vectorized over columns). r is in meters.
    D: reference radius of interface [m]
    nmax: maximum computational harmonic degree
    opt
    ::
    szm: [Latitude dimension, longitude dimension] of tesseroids [deg]
    sigma_nm:   noise std per (n,m) (scalar or same size as Qnm)   [default 1]
    r_nm:       Regularized spectrum per (n,m) (scalar or same size)  [default 1]
    M: scaling constant of mass [m^3 kg^-1 s^-2]
    a: scaling constant of radius [m]
    maxit:      max iterations 
    tolH:       stopping tol on H (max abs change)
    extend:     grid clip extension cells for edge handling
    gamma: largest radial power used
    HRange:     [Hmin Hmax] bracket for H solver (zeta2H)
    mu:         regularization parameter
    nc:         cutoff degree for weighted filter
    Nq: the number of quadrature nodes, if is empty, the automatic determination.
    Hcheck:     reference H for RMS tracking (scalar or grid)       [default 0]
Output:
    Hk              [nthe×nlam] inverted relief grid
    iterCorrection  [k×1] max-abs correction per iteration
    iterRMS         [k×1] RMS(Hk - Hcheck) per iteration
%}

arguments
    lon (1, :)
    lat (:, 1)
    Qnm (:, :)
    funRho  % function handle
    D (1, 1)
    nmax (1, 1)
    opt.szm = [] % (1, 2)
    opt.sigma_nm (:, :) = 1
    opt.r_nm (:, :) = 1
    opt.M (1, 1) = 5.965e24
    opt.a (1, 1) = 6371 * 1e3
    opt.maxit = 200
    opt.tolH = 1e-2
    opt.extend = 0;
    opt.gamma = nmax + 2
    opt.HRange = []
    opt.mu = 0
    opt.nc = []
    opt.Nq = 2
    opt.Hcheck = 0
end

if isempty(opt.HRange)
    opt.HRange(2) = opt.a - D;
    opt.HRange(1) = -2 * opt.HRange(2);
end
if isempty(opt.szm)
    dlat = abs(lat(2) - lat(1));
    dlon = abs(lon(2) - lon(1));
else
    dlat = opt.szm(1);
    dlon = opt.szm(2);
end
if isscalar(opt.r_nm)
    opt.r_nm = repmat(opt.r_nm, size(Qnm));
end
if isscalar(opt.sigma_nm)
    opt.sigma_nm = repmat(opt.sigma_nm, size(Qnm));
end

hdclp = @(x) gridclip(x, dlon, dlat, opt.extend);

[xq, Wq] = glnw(opt.Nq, -1, 1);
[n, m] = creat_nm(nmax);
shs = @(x) shs_grid([n, m, real(x), imag(x)], lon, lat, nmax, 0);

m = (0 : nmax)';
ll = (0 : opt.gamma)';

% Analytic longitude integral (per order m)
lon = m .* (lon - dlon / 2);
costhe2 = cosd( min(180, 90 - lat + dlat / 2) );
costhe1 = cosd( max(0, 90 - lat - dlat / 2) );
Emj = (cosd(m * dlon) + 1i * sind(m * dlon) - 1) ./ (1i * m) .* (cosd(lon) + 1i * sind(lon));
Emj(1, :) = deg2rad(dlon);
Emj = Emj.';

beta = 4 * pi * D ^ 2 ./ (opt.M * (2 * n + 1)) .* (D / opt.a) .^ n;
if isempty(opt.nc)
    omega = (1 + opt.mu .* opt.r_nm .* opt.sigma_nm .^ 2 ./ beta .^ 2) .^ (-1);
elseif isinf(opt.nc)
    omega = 1;
else
    ind = n == opt.nc;
    r_nc = rms(opt.r_nm(ind));
    sigma_nc = rms(opt.sigma_nm(ind));
    mu = (4 * pi * D ^ 2 / (opt.M * (2 * opt.nc + 1)) * (D / opt.a) ^ opt.nc) ^ 2 ./ r_nc ./ sigma_nc .^ 2;
    omega = (1 + mu .* opt.r_nm .* opt.sigma_nm .^ 2 ./ beta .^ 2) .^ (-1);
end

zeta1_0nm = omega .* Qnm ./ beta;

FindH = @(x) zeta2H(x, funRho, opt.HRange, tol=opt.tolH, Nq=opt.Nq);
Hk = FindH(shs(zeta1_0nm));
[nthe, nlam] = size(Hk);
n = m;
Cn2l = Cnk(n + 2, ll);
Hk0 = Hk;

iterCorrection = zeros(opt.maxit, 1);
iterRMS = zeros(opt.maxit, 1);

for k = 1 : opt.maxit
    rnm = zeros((nmax + 2) * (nmax + 1) / 2, 1, 'like', 1 + 1i);
    Rnm = rnm;
    for i = 1 : nthe
        zetal = zeros(opt.gamma + 1, nlam);
        zeta = zeros(nmax + 1, nlam);
        Hki = Hk(i, :) / 2;
        rq = Hki .* (xq + 1);
        rhoq = Wq .* funRho(rq, i, []);
        rq = rq / D;
        for q = 1 : opt.Nq
            zetal = zetal + rhoq(q, :) .* rq(q, :) .^ ll;
        end
        zetal = zetal .* Hki;
        for l = 1 : opt.gamma
            zeta = zeta + zetal(l + 1, :) .* Cn2l(:, l + 1);
        end
        inds = 1;
        zeta = zeta * Emj;
        for nn = 0 : nmax
            inde = inds + nn;
            rnm(inds : inde) = zeta(nn + 1, 1 : nn + 1);
            inds = inde + 1;
        end
        rnm = rnm .* IPnm(nmax, costhe1(i), costhe2(i));
        Rnm = Rnm + rnm;
    end
    Rnm = Rnm / (4 * pi);

    zetak_0nm = zeta1_0nm - omega .* Rnm;
    Hk = FindH(shs(zetak_0nm));
    dH = Hk - Hk0;
    Hk0 = Hk;
    iterCorrection(k) = max(abs(hdclp(dH)), [], "all");
    iterRMS(k) = rms(hdclp(Hk - opt.Hcheck), "all");
    disp(iterCorrection(k))
    if iterCorrection(k) < opt.tolH
        break
    end
end

iterCorrection(k + 1 : end) = [];
iterRMS(k + 1 : end) = [];

end