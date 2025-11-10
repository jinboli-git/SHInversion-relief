function rho = funRho_piecelinear(r, rho0, rho1, rlay, i, j)
% funRho_piecelinear  Evaluate a piecewise-linear radial density profile.
%
% USAGE
%   rho = funRho_piecelinear(r, rho0, rho1, rlay)
%   rho = funRho_piecelinear(r, rho0, rho1, rlay, i, j)
%
% INPUTS
%   r     : query radii; scalar/array. If i,j empty, expanded to [row×col].
%   rho0  : [row×col×lay] intercepts of linear density in each layer
%   rho1  : [row×col×lay] slopes of linear density in each layer
%   rlay  : [row×col×(lay+1)] layer boundaries (descending or ascending);
%           layer k is defined by rlay(:,:,k+1) <= r <= rlay(:,:,k)
%   i, j  : optional row/column index to select a subset before evaluation
%
% OUTPUT
%   rho   : density at r, same size as r (or [row×col] if r scalar and i,j empty)
arguments
    r (:, :)
    rho0 (:, :, :)
    rho1 (:, :, :)
    rlay (:, :, :)
    i = []
    j = []
end

[row, col, lay] = size(rho0);

if isempty(i) && isempty(j)
    r = zeros(row, col) + r;
elseif ~isempty(i) && ~isempty(j)
    rho0 = rho0(i, j, :);
    rho1 = rho1(i, j, :);
    rlay = rlay(i, j, :);
elseif ~isempty(i) && isempty(j)
    rho0 = rho0(i, :, :);
    rho1 = rho1(i, :, :);
    rlay = rlay(i, :, :);
else
    rho0 = rho0(:, j, :);
    rho1 = rho1(:, j, :);
    rlay = rlay(:, j, :);
end
rho = zeros(size(r));
t = ones(size(r));
for k = 1 : lay
    ind = rlay(:, :, k + 1) <= r & r <= rlay(:, :, k);
    rho0k = rho0(:, :, k) .* t; rho1k = rho1(:, :, k) .* t;
    rho(ind) = rho0k(ind) + rho1k(ind) .* r(ind);
end

end

