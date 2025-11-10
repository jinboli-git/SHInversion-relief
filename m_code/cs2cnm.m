function vec = cs2cnm(mat, lmax)
% cs 2 clm
% to Colombo ordering vectors [l m Clm Slm]
% IN:
%    mat ..... CS, SC or degree vector
%    lmax .... Maximum degree of the Spherical Harmonic development
%
% OUT: 
%    vec ..... [l m Clm Slm] vector with Colombo ordering

vec = (1:lmax+1)'*ones(1,lmax+1);
ind = find(tril(vec)~=0);
l   = vec(ind)-1;
vec = vec';
m   = vec(ind)-1;
clm = mat(ind);
mat(end,:) = [];
mat = [zeros(1,lmax+1);mat]';
slm = mat(ind);
vec = sortrows([l m clm slm], 1);
end