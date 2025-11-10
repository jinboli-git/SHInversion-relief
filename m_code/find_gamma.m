function best_gamma = find_gamma(H, D, N, tol)
% find_gamma  Pick the smallest integer γ ∈ [0, N+2] such that an a priori
%             bound ε_γ < tol
%
% USAGE
%   best_gamma = find_gamma(H, D, N)
%   best_gamma = find_gamma(H, D, N, tol)
%
% INPUTS
%   H   [M×N]    Interface relief (can be positive or negative).
%   D   scalar   Reference radius (same unit as H).
%   N   scalar   Maximum harmonic degree
%   tol scalar   Tolerance for the bound ε_γ (default 1e-4).
%
% OUTPUT
%   best_gamma   The first γ (0..N+2) with ε_γ < tol.
%
arguments
    H (:, :)
    D (1, 1)
    N (1, 1)
    tol (1, 1) = 1e-4
end

ind1 = find(H >= 0); ind2 = find(H < 0);
if isempty(ind1)
    Hmax1 = 0;
else
    Hmax1 = max(H(ind1));
end
if isempty(ind2)
    Hmax2 = 0;
else
    Hmax2 = max(-H(ind2));
end
t = (D / (D - Hmax2)) ^ (N + 2);
for gamma = 0 : N + 2
    bio = Cnk(N + 2, gamma + 1);
    supU1 = bio * (Hmax1 / (D + Hmax1)) ^ (gamma + 1);
    supU2 = bio * t * (Hmax2 / D) ^ (gamma + 1);
    epsilon_gamma = max(supU1, supU2);
    if epsilon_gamma < tol
        best_gamma = gamma;
        break
    end
end
