function pnm = Pnm_Bel(nmax, t)
%{
    Calculate the fully normalized Pnm using the column-wise recurrence formula
Input:
    nmax: maximum computational degree/order
    t: t = cos(x)
Return:
    pnm: (nmax + 1, length(x))
Reference:
    Xing, Z., Li, S., Tian, M., Fan, D., & Zhang, C. (2020). Numerical experiments on column-wise recurrence formula to compute fully normalized associated Legendre functions of ultra-high degree and order. Journal of Geodesy, 94(1), 2. https://doi.org/10.1007/s00190-019-01331-0
%}

arguments
    nmax (1, 1)
    t (1, :)
end

lent = length(t);
pnm = zeros((nmax + 1) * (nmax + 2) / 2, lent);
pnm(1, :) = 1; u = sqrt(1 - t .^ 2); pnm(2, :) = sqrt(3) * t; pnm(3, :) = sqrt(3) * u;

if nmax == 1
    return
elseif nmax == 0
    pnm = pnm(1, :);
else
    pnm = Pnm(nmax, t, u, pnm);
end

end
%% subroutine
function pnm = Pnm(nmax, t, u, pnm)
for n = 2 : nmax
    an = sqrt( (2 * n + 1) / (2 * n - 1) );
    bn = sqrt( 2 * (n - 1) * (2 * n + 1) / (2 * n - 1) / n );
    k = n * (n + 1) / 2 + 1;
    k1 = n * (n - 1) / 2 + 1;
    k2 = k1 + 1;
    pnm(k, :) = an * t .* pnm(k1, :) - bn * u / 2 .* pnm(k2, :);
    for m = 1 : n
        k = n * (n + 1) / 2 + m + 1;
        k1 = n * (n - 1) / 2 + m + 1;
        k2 = k1 + 1;
        k3 = k1 - 1;
        cnm = an / n * sqrt(n * n - m * m);
        dnm = an / n / 2 * sqrt((n - m) * (n - m - 1));
        enm = an / n / 2 * sqrt((n + m) * (n + m - 1));
        if(m == 1)
            enm = enm * sqrt(2);
        end
        pnm(k, :) = cnm * t .* pnm(k1, :) - dnm * u .* pnm(k2, :) + enm * u .* pnm(k3, :);
    end
end
end