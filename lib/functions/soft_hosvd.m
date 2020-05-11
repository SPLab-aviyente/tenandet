function [X] = soft_hosvd(Y, Lam1, tau)
% [X] = soft_hosvd(Y, Lam, tau)
% Function that returns tensor whose singular values are soft thresholded
% using parameter tau at each mode.
N = ndims(Y);
sz = size(Y);
X = cell(1,N);
for i = 1:N
    X{i} = soft_moden(Y+Lam1{i}, tau, i );
    X{i} = m2t(X{i}, sz, i);
end
end

function X = soft_moden(T, tau, n)
[U, S, V] = svd(t2m(T, n), 'econ');
s = diag(S)-tau;
smask = s>0;
S = diag(s(smask));

X = U(:, smask)*S*V(:, smask)';
end