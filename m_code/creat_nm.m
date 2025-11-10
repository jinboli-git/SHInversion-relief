function [n, m] = creat_nm(nmax, typ)
%{
    Calculate the degree and order of the nmax spherical harmonic coefficients
Input:
    nmax:  maximum computational degree/order
    typ: the returned sorting order: Sort by order (cnm_ord), Sort by degree (cnm_degree)
Return:
    n: (:, 1) degree
    m: (:, 1) order
%}

arguments
    nmax (1, 1)
    typ {mustBeMember(typ, {'cnm_deg', 'cnm_ord'})} = 'cnm_deg'
end

n = zeros((nmax + 1) * (nmax + 2) / 2, 1);
m = n;

rows = nmax + 1;
first = 1;
last = nmax + 1;

for i = 1:rows
    n(first:last) = ((i - 1) : (rows - 1))';
    m(first:last, 1) = (i - 1) * ones(length(first:last), 1);
    first = last + 1;
    last = last + 1 + nmax - i;
end

if strcmp(typ, 'cnm_deg')
    out = [n, m];
    out = sortrows(out, 1);
    n = out(:, 1);
    m = out(:, 2);
end

end