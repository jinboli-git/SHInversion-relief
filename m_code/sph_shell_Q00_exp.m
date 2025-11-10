function Q00 = sph_shell_Q00_exp(H, D, rho, M)
% sph_shell_Q00_exp  Degree-0 (isotropic) potential of a spherical shell
%                    with exponential radial density.
%
% USAGE
%   Q00 = sph_shell_Q00_exp(H, D, rho, M)
%
% INPUTS
%   H   (scalar)        Shell thickness (same unit as D).
%   D   (scalar)        Reference radius (e.g., shell inner radius).
%   rho (1×2 vector)    Parameters of exponential density:
%                         rho0 = rho(1), kappa = rho(2),
%                         Δρ(r') = rho0 * exp(kappa * r'),  r'∈[0,H].
%   M   (scalar)        Mass scaling constant.
%
% OUTPUT
%   Q00 (scalar)        Degree-0 potential coefficient.
%
% FORMULA
%   Using binomial terms l = 0..2 and the integral
%      ∫_0^H r'^l exp(κ r') dr' = κ^(-l-1) * Γ(l+1) * [gammainc(-κH,l+1,'upper') - 1],
%   where MATLAB's gammainc(·,a,'upper') is the **normalized** upper incomplete
%   gamma Γ(a,x)/Γ(a). High precision (VPA) is used for numerical stability,
%   then cast back to double.
% NOTE
%   Cnk(n,k) is assumed to return the binomial coefficient "n choose k".
arguments
    H (1, 1)
    D (1, 1)
    rho (1, 2)
    M (1, 1)
end

digits(100);
rho0 = vpa(rho(1));
kappa = rho(2);
D = vpa(D);

tl = 0;
for l = 0 : 2
    ts = rho0  * Cnk(2, l) * (-D) ^ (-l) * gamma(l + 1) * kappa ^ (-l - 1) * (gammainc(-kappa * H, l + 1, "upper") - 1);
    tl = tl + ts;
end
Q00 = 4 * pi / M * D ^ 2 * tl;
Q00 = double(Q00);
end
