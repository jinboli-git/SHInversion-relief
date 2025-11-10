function H = zeta2H(zeta, funRho, range, opt)
%{
    Solve for relief H from target zeta by 1D root finding (per element).
Input:
   zeta: [M×N double] target zeta values for each grid cell.
   funRho: handle, funRho(r, i, j) -> density at (r - D) for cell (i,j).
%            r is radius can be scalar or array, must return same size as r.
   range: bracket [Hmin, Hmax] for the root of f(H)=0.
   opt
    ::
   tol:  tolerance of H
   Nq: Gaussian quadrature order in H2zeta.
Output:
    H: [M×N double] solution per cell.
%}
arguments
    zeta (:, :)
    funRho
    range (1, 2)
    opt.tol (1, 1) = 1e-2;
    opt.Nq = 2;
end

% Preallocate output
[numRow, numCol] = size(zeta);
H = zeros(numRow, numCol);
f = @(x, i, j) zeta(i, j) - H2zeta(x, @(x) funRho(x, i, j), opt.Nq);
opt = optimset('TolX', opt.tol, 'MaxIter', 100);
% Solve cell-by-cell
for i = 1 : numRow
    for j = 1 : numCol
        fij = @(x) f(x, i, j);
        fij1 = fij( range(1) );
        fij2 = fij( range(2) );
        if fij1 * fij2 >= 0
            H(i, j) = nan;
            if abs(fij2) < abs(fij1)
                H(i, j) = range(2);
            else
                H(i, j) = range(1);
            end
        else
            H(i, j) = fzero(fij, range, opt);
        end
    end
end

end
%% subroutine
function zeta = H2zeta(H, funRho, Nq)

arguments
    H (:, :)
    funRho 
    Nq (1, 1) = 2
end

H = H / 2;
[xq, W] = glnw(Nq, -1, 1);
zeta = 0;

for q = 1 : Nq
    rq = H .* (xq(q) + 1);
    zeta = zeta + W(q) * funRho(rq);
end

zeta = zeta * H;

end
