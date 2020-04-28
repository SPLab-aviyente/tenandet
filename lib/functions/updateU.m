function [U, objVal] = updateU(F, Y, U, modes)
%%
%
n  = length(modes);
N  = ndims(Y);
theta = 10^3;
% indModes = setdiff(1:N,modes);

szs = cell(1,2*n);
for i=1:2*n
    szs{i} = size(U{i});
end

FsU = mergeTensors(F, U(n+1:end), modes(end:-1:1));
Ur = mergeFactors(U(n+1:end));
Ur = Runfold(Ur); Ur = Ur*Ur';
[U, objVal]  = updateFormer(Y, FsU, Ur, U, modes, theta, F, 1);

YU  = mergeTensors(Y, U(1:n), modes);
YU  = permute(YU, [1:N-n, N-n+2, N-n+1]);
[U, val] = updateLatter(F, YU, U, modes, theta, 1);
objVal = [objVal, val];

end

function [U, objAll] = updateFormer(Y, FsU, Ur, U, modes, theta, F, store_convergence)
%
%
if nargin < 8
    store_convergence=false;
end
n  = length(modes);
N  = ndims(Y);
szs = cell(1,2*n);
for j=1:2*n
    szs{j} = size(U{j});
end

if store_convergence
    objT = ipermute(squeeze(mergeFactors({squeeze(mergeTensors(Y, U(1:n), ...
        modes)), U{n+1:end}})),[setdiff(1:N, modes), modes(end:-1:1)])- F;
    objAll = obj_func(objT, U, theta);
end
for i=1:n
    % Find the mode shift of the current mode after the projection.
    ni = modes(i)-i+1;                                                          
    if i == 1
        Ytemp = mergeTensors(Y, U(2:n), modes(2:n));
        setInd = setdiff(1:N-n+3, [ni, N-n+2]);
        A = mergeTensors(mergeTensors(Ytemp, Ur, N-n+3, 1),...
            Ytemp, setInd, setInd);
        A = mMat(A, 2);
        b = mergeTensors(Ytemp, FsU, setInd, 1:N-n+1);
        U{1} = pinv(A+theta*eye(size(A)))*b(:);
    elseif i == n
        Ytemp = mergeTensors(Y, U(1:n-1), modes(1:n-1));
        Ytemp = permute(Ytemp, [1:N-n+1, N-n+3, N-n+2]);
        setInd = setdiff(1:N-n+2, [ni,N-n+2]);
        b = mergeTensors(Ytemp, FsU, setInd, 1:N-n);
        b = permute(b, [2,1,3]);
        A = mergeTensors(Ytemp, Ytemp, setInd, setInd);
        A = mMat(permute(A, [2,1,4,3]),2);
%         tic
%         A = pinv(kron(Ur, A)+theta*eye(numel(b)));
%         toc
%         tic
        A     = comp_inv_frmr_n(Ur, A, theta);
%         toc
        U{n}  = A*b(:);
%         U{n}  = reshape(U{n},szs{n});
    else
        Ytemp = mergeTensors(mergeTensors(Y, U(1:i-1), modes(1:i-1)), ...
            U(i+1:n), modes(i+1:n)-i+1);
        Ytemp = permute(Ytemp, [1:N-n+1, N-n+3:N-n+5, N-n+2]);
        setInd = setdiff(1:N-n+4, [ni,N-n+2,N-n+3]);
        A = mergeTensors(mergeTensors(Ytemp, Ur, N-n+4, 1),...
            Ytemp, setInd, setInd);
        A = permute(A, [2,1,3,5,4,6]);
        A = mMat(A, 3);
        b = mergeTensors(Ytemp, FsU, setInd, 1:N-n+1);
        b = permute(b, [2,1,3]);
        U{i} = pinv(A+theta*eye(size(A)))*b(:);
%         U{i} = reshape(U{i},szs{i});
    end
    [Q, R] = qr(reshape(U{i}, prod(szs{i}(1:2)), []), 0);
    U{i} = reshape(Q, szs{i});
    U{i+1} = reshape(R*Runfold(U{i+1}, 1), szs{i+1});
    if store_convergence
        objT = ipermute(squeeze(mergeFactors({squeeze(mergeTensors(Y, U(1:n), ...
            modes)), U{n+1:end}})),[setdiff(1:N, modes), modes(end:-1:1)])- F;
        objAll = [objAll, obj_func(objT, U, theta)];
    end
end
end

function [U, objAll] = updateLatter(F, YU, U, modes, theta, store_convergence)
%
%
if nargin < 6
    store_convergence=false;
end
n  = length(modes);
N  = ndims(F);
Ytemp  = squeeze(mergeTensors(YU, YU, 1:N-n, 1:N-n));
szs = cell(1,2*n);
for j=1:2*n
    szs{j} = size(U{j});
end
objAll = [];
for i=n:-1:1
    ni = modes(n+1-i)-n+i;                                                          % Find the mode shift of the current mode after the projection.
    if i == 1
        Ftemp   = mergeTensors(F, U(n+2:end), modes(n-1:-1:1));
        b       = mergeTensors(YU, Ftemp, 1:N-n, [1:ni-1, ni+1:N-n+1]);
        A       = mergeFactors(U(n+2:end));
        A       = mergeTensors(A, A, 2:n, 2:n);
        A       = comp_inv_lat_1(A, Ytemp, theta, szs{n+1}(2));
%         A       = pinv(kron(A,kron(eye(szs{n+1}(2)),Ytemp))+theta*eye(numel(b)));
    elseif i == n
        Ftemp   = mergeTensors(F, U(n+1:end-1), modes(n:-1:2));
        b       = mergeTensors(YU, Ftemp, 1:N-n+1, [1:ni-1, ni+1:N-n+2])';
        Utemp   = mergeFactors(U(n+1:end-1));
        Utemp   = mergeTensors(YU, Utemp, N-n+1, 1);
        Utemp   = Lunfold(Utemp, N);
        A       = comp_inv_lat_n(Utemp, szs{2*n}(2), theta);
%         A = pinv(kron(eye(szs{2*n}(2)), Utemp'*Utemp)+theta*eye(numel(b)));
%         U{2*n}  = reshape(U{2*n}, szs{2*n});
    else
        Ftemp   = mergeTensors(mergeTensors(F, U(n+(1:i-1)), ...
            modes(n+1-(1:i-1))), U(n+(i+1:n)), modes(n+1-(i+1:n)));
        b       = mergeTensors(YU, Ftemp, 1:N-n+1, [1:ni-1, ni+1:N-n+2]);
        b       = permute(b, [2,1,3]);
        A       = mergeFactors(U(n+1:n+i-1));
        A       = Lunfold(mergeTensors(Ytemp, A, 1, 1), i+1)'*Lunfold(A, i+1);
        Utemp   = Runfold(mergeFactors(U(n+i+1:n*2)));
        A       = comp_inv_lat(Utemp*Utemp', A, szs{n+i}(2), theta);
%         A = pinv(kron(Utemp*Utemp', kron(eye(szs{n+i}(2)), A))+theta*eye(numel(b)));
    end
    U{n+i} = A*b(:);
    U{n+i} = reshape(U{n+i}, szs{n+i});
    [Q, R] = qr(reshape(U{n+i}, szs{n+i}(1), [])', 0);
    U{i+n} = reshape(Q', szs{i+n});
    U{i+n-1} = reshape(Lunfold(U{i+n-1}, 3)*R', szs{i+n-1});
    if store_convergence
        if i == 1
            YU = mergeTensors(YU,R, N-n+1, 1);
        end
        objT    = ipermute(squeeze(mergeFactors({mergeTensors(YU, U{n+1}, ...
            N-n+1, 1), U{n+2:end}})), [setdiff(1:N, modes), modes(end:-1:1)])-F;
        objAll = [objAll, obj_func(objT, U, theta)];
    end
end

end

function In = comp_inv_lat(M1, M2, s1, theta)
% In = comp_inv_lat(M1, M2, s1, theta)
%   Compute inverse of the intermediary matrices in the latter half.
%
[U1, S1, ~] = svd(M1, 'econ');
[U2, S2, ~] = svd(M2, 'econ');

U = kron(U1, kron(eye(s1), U2));
S = kron(S1, kron(eye(s1), S2))+theta*eye(size(U,2));
S = diag(diag(S).^-1);
In = U*S*U';
end

function In = comp_inv_lat_1(M1, M2, theta, s1)
% In = comp_inv_lat_1(M1, M2, theta, s1)
%   Compute inverse of the first matrix in the latter half.
%
[U1, S1, ~] = svd(M1, 'econ');
[U2, S2, ~] = svd(M2, 'econ');

U = kron(U1, kron(eye(s1), U2));
S = kron(S1, kron(eye(s1), S2))+theta*eye(size(U,2));
S = diag(diag(S).^-1);
In = U*S*U';
end

function In = comp_inv_lat_n(M1, s1, theta)
% In = comp_inv_lat_n(M1, s1, theta)
%   Compute inverse of the last matrix in the latter half.
%
[U1, S1, ~] = svd(M1', 'econ');

U = kron(eye(s1), U1);
S = kron(eye(s1), S1.^2)+theta*eye(size(U,2));
S = diag(diag(S).^-1);
In = U*S*U';
end

function In = comp_inv_frmr_n(M1, M2, theta)
%   Compute inverse of the nth matrix in the former half.
%
if M1'*M1 == eye(size(M1,2))
    U1 = M1;
    S1 = eye(size(M1,2));
else
    [U1, S1, ~] = svd(M1,'econ');
end

[U2, S2, ~] = svd(M2,'econ');

U = kron(U1, U2);
S = kron(S1, S2)+theta*eye(size(U,2));
S = diag(diag(S).^-1);
In = U*S*U';
end

function val = obj_func(objT, U, theta)

uval = 0;
for i=1:length(U)
    uval = uval+theta/2*norm(ten2vec(U{i}))^2;
end
val = norm(objT(:))^2+uval;
end