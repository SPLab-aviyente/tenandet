function T = m2t(M, sz, n)
% T = m2t(M, sz, n)
% Matrix-to-Tensor reshaping operator.
mode_row = setdiff(1:length(sz),n);
T = reshape(M,[sz(n), sz(mode_row)]);
T = ipermute(T, [n, mode_row]);
end