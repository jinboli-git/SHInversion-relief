function Q00 = sph_shell_Q00_poly(H, D, rho, M)
% sph_shell_Q00_poly  Compute the isotropic (degree-0) potential of a spherical shell.
%
% USAGE
%   Q00 = sph_shell_Q00_poly(H, D, rho, M)
%
% INPUTS
%   H   (scalar)   Shell thickness [same unit as D].
%   D   (scalar)   Reference radius (e.g., shell inner radius).
%   rho (1×P+1)    Radial density polynomial coefficients rho_s (s=0..P),
%                  i.e.,  Δρ(r) = Σ_{s=0..P} rho(s+1) * (r - D)^s  (or equivalent).
%   M   (scalar)   Mass scaling constant.
%
% OUTPUT
%   Q00 (scalar)   Degree-0 potential coefficient (C00-like scale): 4π D^2 / M * Σ_l ...
%
% NOTES
%   - High-precision arithmetic (VPA, 100 digits) is used to stabilize
%     small-difference powers and divisions by D^l; final value is cast to double.
%   - Cnk(n,k) is assumed to return the binomial coefficient "n choose k".
%
% FORMULA (implemented)
%   Q00 = (4π/M) * D^2 * Σ_{l=0}^{2} [ C(2,l)/D^l * Σ_{s=0}^{P} rho_{s} * H^{l+s+1}/(l+s+1) ]
%   where P = numel(rho)-1.
%
arguments
    H (1, 1)
    D (1, 1)
    rho (1, :)
    M (1, 1)
end

digits(100);
H = vpa(H);
D = vpa(D);

P = length(rho) - 1;
tl = 0;
for l = 0 : 2
    ts = 0;
    for s = 0 : P
        ts = ts + rho(s + 1) / (l + s + 1) * H ^ (l + s + 1);
    end
    tl = tl + ts * Cnk(2, l) / D ^ l;
end
Q00 = 4 * pi / M * D ^ 2 * tl;
Q00 = double(Q00);
end
