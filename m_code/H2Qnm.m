function Qnm = H2Qnm(lon, lat, H, funRho, D, nmax, opt)
%{
    Potential spherical-harmonic coefficients (multipoles) of tesseroids grid with
    general radial density using Gauss–Legendre quadrature (fixed Nq).
Input:
    lon: the longitude vector of the center of tesseroids [deg]
    lat: the latitude vector of the center of tesseroids [deg]
    H: tesseroid thickness (interface relief) per cell [m]
    funRho: function handle, funRho(r,i,~) -> density at (r + D) for row i
           (vectorized over columns). r is in meters.
    D: reference radius of interface [m]
    nmax: maximum computational harmonic degree
    opt
    ::
    M: scaling constant of mass [m^3 kg^-1 s^-2]
    a: scaling constant of radius [m]
    szm: [Latitude dimension, longitude dimension] of tesseroids [deg]
    gamma: largest radial power used
    Nq: the number of quadrature nodes, if is empty, the automatic determination.
Return:
    Qnm: [ (nmax+1)(nmax+2)/2 × 1 complex] multipoles
%}

arguments
    lon (1, :)
    lat (:, 1)
    H (:, :)
    funRho % function handle
    D (1, 1)
    nmax (1, 1)
    opt.M (1, 1)
    opt.a (1, 1)
    opt.szm = [] % (1, 2)
    opt.gamma = nmax + 2
    opt.Nq = 2;
end

nthe = length(lat);
nlam = length(lon);

H = H / 2;
if isscalar(H)
    H = repmat(H, nthe, nlam);
end

n = (0 : nmax)';
m = n;
ll = (0 : opt.gamma)';

% Binomial coefficients C(n+2, l)
Cn2l = Cnk(n + 2, ll);

[xq, Wq] = glnw(opt.Nq, -1, 1);

if isempty(opt.szm)
    dlat = abs(lat(2) - lat(1));
    dlon = abs(lon(2) - lon(1));
else
    dlat = opt.szm(1);
    dlon = opt.szm(2);
end

% Analytic longitude integral (per order m)
lon = m .* (lon - dlon / 2);
costhe2 = cosd( min(180, 90 - lat + dlat / 2) );
costhe1 = cosd( max(0, 90 - lat - dlat / 2) );
Emj = (cosd(m * dlon) + 1i * sind(m * dlon) - 1) ./ (1i * m) .* (cosd(lon) + 1i * sind(lon));
Emj(1, :) = deg2rad(dlon);
Emj = Emj.';

qnm = zeros((nmax + 2) * (nmax + 1) / 2, 1, 'like', 1 + 1i);
Qnm = qnm;

for i = 1 : nthe
    zetal = zeros(opt.gamma + 1, nlam);
    zeta = zeros(nmax + 1, nlam);
    Hi = H(i, :);
    rq = Hi .* (xq + 1);
    rhoq = Wq .* funRho(rq, i, []);
    rq = rq / D;
    for q = 1 : opt.Nq
        zetal = zetal + rhoq(q, :) .* rq(q, :) .^ ll;
    end
    zetal = zetal .* Hi;
    for l = 0 : opt.gamma
        zeta = zeta + zetal(l + 1, :) .* Cn2l(:, l + 1);
    end
    inds = 1;
    zeta = zeta * Emj;
    for nn = 0 : nmax
        inde = inds + nn;
        qnm(inds : inde) = zeta(nn + 1, 1 : nn + 1);
        inds = inde + 1;
    end
    qnm = qnm .* IPnm(nmax, costhe1(i), costhe2(i));
    Qnm = Qnm + qnm;
end

[n, ~] = creat_nm(nmax);
Qnm = Qnm ./ (2 * n + 1) .* (D / opt.a) .^ n * (D ^ 2 / opt.M);

end
