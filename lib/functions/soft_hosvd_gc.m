function [L, nuc_norm] = soft_hosvd_gc(Y, Lam1, Sig, U, psi, beta_1)
% [L, nuc_norm] = soft_hosvd_gc(Y, Lam1, Sig, U, psi, beta_1)
% Function that returns the update to nuclear norm minimization of graph
% core using ADMM.
N = ndims(Y);
sz = cellfun(@size, U, num2cell(2*ones(1,length(U))));
L = cell(1,N);
nuc_norm = L;
for i = 1:N
    tmp = beta_1*tmprod(Y+Lam1{i}, U, 1:4, 'T');
    [tempL, nuc_norm{i}] = soft_mode_wtd(tmp, psi(i).*(Sig{i}.^-1), i );
    nuc_norm{i} = nuc_norm{i}/(beta_1);
    tempL = m2t(tempL, sz, i)/(beta_1);
    L{i} = tmprod(tempL, U, 1:4);
end
end

function [X, nuc_norm] = soft_mode_wtd(T, tau, n)
[U, S, V] = svd(t2m(T, n), 'econ');
s = diag(S)-tau;
smask = s>0;
S = diag(s(smask));
nuc_norm = sum(s(smask));
X = U(:, smask)*S*V(:, smask)';
end