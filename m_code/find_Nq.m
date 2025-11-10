function Nq = find_Nq(H, funRho, gamma, D, tol)
% find_Nq  Choose an adequate Gauss–Legendre quadrature order Nq
%           for ∫_D^(D+H) ρ(r) (r-D)^γ dr, by row-wise adaptive refinement.
% USAGE
%   Nq = find_Nq(H, funRho, gamma, D)
%   Nq = find_Nq(H, funRho, gamma, D, tol)
%
% INPUTS
%   H      [R×C]   relief per grid cell (same units as D).
%   funRho handle  funRho(r,i,~) → density evaluated at (r+D) for row i.
%                  r is 1×C (or broadcastable); i is scalar row index.
%   gamma  scalar  Exponent applied to (r/D)^gamma inside the integral.
%   D      scalar  Reference radius for nondimensionalization.
%   tol    scalar  Relative change tolerance between consecutive orders (default 1e-4).
%
% OUTPUT
%   Nq     scalar  Selected quadrature order (max over rows), with a floor of 2.

arguments
    H (:, :)
    funRho 
    gamma (1, 1)
    D (1, 1)
    tol (1, 1) = 1e-4
end

H2 = H / 2;
numts = 50;
[row, col] = size(H2);
Nq = zeros(row, 1);
for i = 1 : row
    vv = zeros(numts, col);
    H2i = H2(i, :);
    for nq = 1 : numts
        t = 0;
        [xq, Wq] = glnw(nq, -1, 1);
        for q = 1 : nq
            rq = H2i .* (xq(q) + 1);
            t = t + Wq(q) * funRho(rq, i, []) .* (rq / D) .^ gamma;
        end
        vv(nq + 1, :) = t;
        rerr = max(abs((t - vv(nq, :)) ./ t));
        if isnan(rerr) || rerr <= tol
            break;
        end
    end
    Nq(i) = nq;
end

Nq = max(Nq);
if Nq < 2
    Nq = 2;
end

end