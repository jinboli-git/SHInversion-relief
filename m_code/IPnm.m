function ipnm = IPnm(nmax, t1, t2)
%{
    The integral of the associated Legendre function \int^t2_t1 Pnm(t) dt is recursively calculated using the cross-order method
Input:
    nmax maximum computational degree/order
    t t = cos(x)
Return:
    P ((n * n + 3 * n + 2) / 2, 1) the integral of Pnm in cnm-format sort by order
Reference:
    邢志斌. (2019). GOCE卫星重力梯度数据恢复地球重力场理论与方法研究. 博士论文
%}

arguments
    nmax (1, 1)
    t1 (1, :)
    t2 (1, :)
end

len = length(t1);
ipnm = zeros((nmax + 1) * (nmax + 2) / 2, len);
u1 = sqrt(1 - t1 .^ 2); u2 = sqrt(1 - t2 .^ 2);

ipnm(1, :) = t2 - t1;

if nmax == 0
    ipnm = -ipnm;
    return
end

ipnm(2, :) = sqrt(3) / 2 * (t2 .^ 2 - t1 .^ 2);
ipnm(3, :) = sqrt(3) / 2 * (t2 .* u2 - t1 .* u1 - acos(t2) + acos(t1));

if nmax == 1
    ipnm = -ipnm;
    return
end

u1 = u1 .* u1; u2 = u2 .* u2;
pn01 = u1 .* Pnm(nmax, t1, 0)'; pn11 = u1 .* Pnm(nmax, t1, 1)';
pn02 = u2 .* Pnm(nmax, t2, 0)'; pn12 = u2 .* Pnm(nmax, t2, 1)';

indi = 3;
for n = 2 : nmax
    % m = 0
    anm = sqrt( (2 * n - 1) * (2 * n + 1) / n / n );
    anm1 = sqrt( (2 * n - 1) * (2 * n - 3) / (n - 1) / (n - 1) );
    ipnm(indi + 1, :) = (n - 2) / (n + 1) * anm / anm1 * ipnm(indi - n - n + 2, :) - anm / (n + 1) * (pn02(n, :) - pn01(n, :));
    % m = 1
    anm = sqrt( (2 * n - 1) * (2 * n + 1) / (n - 1) / (n + 1) );
    anm1 = sqrt( (2 * n - 1) * (2 * n - 3) / (n - 2) / n );
    ipnm(indi + 2, :) = (n - 2) / (n + 1) * anm / anm1 * ipnm(indi - n - n + 3, :) - anm / (n + 1) * (pn12(n, :) - pn11(n, :));
    for m = 2 : n
        if(n - m == 0)
            anm = 0;
        else
            anm = sqrt((n + n + 1) * (n - m) * (n - m - 1) / ((n + n - 3) * (n + m) * (n + m - 1)));
        end
        if(m == 2)
            bnm = sqrt(2 * (n + n + 1) * (n + m - 2) * (n + m - 3) / ((n + n - 3) * (n + m) * (n + m - 1)));
            cnm = sqrt(2 * (n - m + 1) * (n - m + 2) / ((n + m) * (n + m - 1)));
        else
            bnm = sqrt((n + n + 1) * (n + m - 2) * (n + m - 3) / ((n + n - 3) * (n + m) * (n + m - 1)));
            cnm = sqrt((n - m + 1) * (n - m + 2) / ((n + m) * (n + m - 1)));
        end
        ipnm(indi + m + 1, :) = anm * ipnm(indi - n - n + m + 2, :) + bnm * ipnm(indi - n - n + m, :) - cnm * ipnm(indi + m - 1, :);
    end
    indi = indi + n + 1;
end

ipnm = -ipnm;

end

%% subroutine
function pnm = Pnm(nmax, t, m)

arguments
    nmax (1, 1)
    t (1, :)
    m (1, 1)
end

lent = length(t);
u = sqrt(1 - t .^ 2);
root3 = sqrt(3);
dmax = nmax + 1;

t = t'; u = u';
pnm(lent, dmax) = 0;
mind = m + 1;

if m == 0
    pnn = 1;
else
    pnn = root3 * u;
    for i = 2 : m
        pnn = sqrt(1 + 0.5 / i) * u .* pnn;
    end
end
pnm(:, mind) = pnn;
for n = mind : nmax
    nind = n + 1;
    anm = sqrt((4 * n ^ 2 - 1) / (n ^ 2 - m ^ 2));
    if n == mind
        pnm(:, nind) = anm * t .* pnm(:, n);
    else
        bnm = sqrt( (n + n + 1) * (n + m - 1) * (n - m - 1) / (n + n - 3) / (n ^ 2 - m ^ 2) );
        pnm(:, nind) = anm * t .* pnm(:, n) - bnm * pnm(:, n - 1);
    end
end

end