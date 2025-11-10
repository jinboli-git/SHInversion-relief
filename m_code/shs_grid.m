function fld = shs_grid(SHCs, lon, lat, nmax, nmin, parm)
%{
    Calculate the spherical harmonic synthesis of any rectangular grid
Input:
   SHCs: Spherical harmonic coefficient, [n, m, Cnm, Snm]
   lon: the longitude vector of the computation grid, [deg]
   lat: the latitude vector of the computation grid, [deg]
   nmax: maximum computational degree/order
   nmin: minimum computational degree/order
   parm
    ::
   fldTyp:  the type of calculation
   'shs', only spherical harmonic synthesis [none]
   'v', potiental [m^2 s^-2] 
   'vz', first-order vertical gradient [m s^-2]
   'vzz', second-order vertical gradient [s^-2]
    a: radius of the reference sphere [m]
    GM: G*M constant
    height: computation height [m]
Ref:
    Sneeuw N (1994) Global spherical harmonic analysis by least-squares and numerical quadrature methods in historical perspective. Geophys J Int 118:707–716. https://doi.org/10.1111/j.1365-246X.1994.tb03995.x
%}

arguments
    SHCs (:, :)
    lon (1, :)
    lat (:, 1)
    nmax = []
    nmin = 0;
    parm.fldTyp char {mustBeMember(parm.fldTyp, {'shs', 'v', 'vz', 'vzz'})} = 'shs'
    parm.a (1, 1) = 6371 * 1e3
    parm.GM (1, 1) = 5.965e24 * 6.672e-11
    parm.height (1, 1) = 0
end

nthe = length(lat);
r = parm.a + parm.height;

n = SHCs(:, 1);
nmaxSHCs = n(end);

if isempty(nmax)
    nmax = nmaxSHCs;
end
if nmin > nmax
    nmin = nmax;
end
if nmax > nmaxSHCs
    nmax = nmaxSHCs;
end

SHCs = SHCs(n >= nmin & n <= nmax, :);
n = SHCs(:, 1);
SHCs = SHCs(:, 3 : 4);

costhe = cosd(90 - lat);
Pnm = @(x) Pnm_Bel(nmax, costhe(x));
u = kapa(n, parm.fldTyp, parm.GM, parm.a, r);

SHCs = SHCs .* u;
SHCs = complex(SHCs(:, 1), SHCs(:, 2));

AB = zeros(nthe, nmax + 1, 'like', 1 + 1i);
mlam = lon .* (0 : nmax)';

% 1st step
for ii = 1 : nthe
    pnm = Pnm(ii) .* SHCs;
    tp = zeros(1, nmax + 1, 'like', 1 + 1i); tp(end) = pnm(end);
    for mfor = 0 : nmax - 1
        ind = mfor : nmax; ind = ind .* (ind + 1) / 2 + mfor + 1;
        tp(mfor + 1) = sum(pnm(ind));
    end
    AB(ii, :) = tp;
end

% 2nd step
fld = real(AB) * cosd(mlam) + imag(AB) * sind(mlam);

end
%% subroutine
function u = kapa(n, typ, GM, a, r)

arguments
    n (:, 1)
    typ
    GM (1, 1)
    a (1, 1)
    r (1, 1)
end

ir = 1 / r;
ir2 = ir * ir;
ir3 = ir2 * ir;

logq = log1p((a - r) / r);
qn = exp(n * logq);

switch typ
    case 'shs'
        u = 1;
    case {'v'}
        u = GM * ir;
    case 'vz'
        u = -(n + 1) * (GM * ir2);
    case {'vzz'}
        u = (n + 1) .* (n + 2) * (GM * ir3);
end

u = u .* qn;

end