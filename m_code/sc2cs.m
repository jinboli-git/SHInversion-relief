function [cs, rows, cols] = sc2cs(sc)
% SC2CS(FIELD) converts the rectangular (L+1)x(2L+1) matrix FIELD, containing
% spherical harmonics coefficients in /S|C\ storage format into a 
% square (L+1)x(L+1) matrix in |C\S| format.
%
% IN:
%    field .... the rectangular (L+1)x(2L+1) matrix FIELD, containing
%               spherical harmonics coefficients in /S|C\ storage format
%
% OUT: 
%    cs ....... square (L+1)x(L+1) matrix in |C\S| format
%    rows ..... the number of rows for matrix FIELD
%    cols ..... the number of columns for matrix FIELD
% -------------------------------------------------------------------------
[rows, cols] = size(sc);
if cols == rows
    cs = sc;
else
    lmax = rows-1;
    c    = sc(:,lmax+1:2*lmax+1);
    s    = [zeros(lmax+1,1) sc(:,1:lmax)];
    cs   = tril(c) + triu(rot90(s),1);
end

end