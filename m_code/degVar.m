function rs = degVar(SHCs, nmax)
%{
    Used for calculating the degree variance
Input:
    SHCs: [n, m, Cnm, Snm]
    nmax: maximum computational degree/order
Return:
    rs: (:, 1) degree variance
%}

arguments
    SHCs (:, :)
    nmax = [];
end

nmaxsig = SHCs(end, 1);
SHCs = cnm2sc(SHCs, nmaxsig);

if isempty(nmax)
    nmax = nmaxsig;
end

rs = sum(SHCs .^ 2, 2);
rs = rs(1 : nmax + 1);

end

