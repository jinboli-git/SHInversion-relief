function sc = cnm2sc(cnm, max_nm)
%{ 
     converts coefficients in cnm-format to /S|C\-format
Input:
    cnm: SHCs in cnm-format
    max_nm:  maximum computational degree/order
Return:
    sc: coefficients in /S|C\ format
%}
arguments
    cnm (:, :)
    max_nm (1, 1)
end

sc = zeros(max_nm + 1, 2 * (max_nm + 1) - 1);
idx_s = sub2ind(size(sc), cnm(:, 1) + 1, max_nm + 1 - cnm(:, 2));
idx_c = sub2ind(size(sc), cnm(:, 1) + 1, max_nm + 1 + cnm(:, 2));
sc(idx_s) = cnm(:, 4);
sc(idx_c) = cnm(:, 3);

end