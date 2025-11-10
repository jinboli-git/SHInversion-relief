function c = Cnk(n, k)
%{
    Calculate the binomial coefficients C_n^k without overflow   
%}
arguments
    n (:, 1)
    k (1, :)
end

nk1 = n - k + 1;
nk1(nk1 < 1) = 0;
% Apply log-gamma trick only for valid elements
c = exp(gammaln(n + 1) - gammaln(k + 1) - gammaln(nk1));
% Round to remove floating-point artifacts
c = round(c);

end

