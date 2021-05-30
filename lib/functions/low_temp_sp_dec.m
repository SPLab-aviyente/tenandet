function [L, S, Nt, times] = low_temp_sp_dec(Y, param)
% [S] = low_temp_sp_dec(Y, param)
% Low rank+Temporally Smooth Sparse Decomposition

N = ndims(Y);
sz = size(Y);
mask_Y = ones(sz)>0;
mask_Y(param.ind_m) = ~mask_Y(param.ind_m);
max_iter = param.max_iter;
err_tol = param.err_tol;
alpha = param.alpha;
lambda = param.lambda;
gamma = param.gamma;
psi = param.psi;
beta_1 = param.beta_1;
beta_2 = param.beta_2;
beta_3 = param.beta_3;
beta_4 = param.beta_4;
Lx = cell(1, N);
for i=1:N
    Lx{i} = zeros(size(Y));
end
L = zeros(size(Y));
S = zeros(size(Y));
Nt = zeros(size(Y));
W = zeros(size(Y));
Z = zeros(size(Y));
D = convmtx([1,-1], size(Y,1));
D(:,end) = [];
D(end,1) = -1;
if beta_4~=0 && beta_3~=0
    invD = (beta_4*eye(sz(1))+beta_3*(D'*D))^-1;
else
    invD = zeros(sz(1));
end
Lam{2} = cell(1, N);
for i=1:N
    Lam{2}{i} = zeros(size(Y));
end
Lam{1} = zeros(size(Y));
Lam{3} = zeros(size(Y));
Lam{4} = zeros(size(Y));

times = [];
iter = 1;
% obj_val = compute_obj(Y,Lx,S,N,W,Z,Lam1,Lam2,Lam3,D,lambda,gamma,alpha,beta_1,beta_2,beta_3);
while true
    %% L Update
    tstart = tic;
    temp = zeros(size(Y));
    for i=1:4
        temp = temp + (Lx{i}-Lam{2}{i})*beta_2;
    end
    T1 = Y-S+Lam{1};
    L(mask_Y) = (beta_1*T1(mask_Y)+temp(mask_Y))/(beta_1+4*(beta_2));
    L(~mask_Y) = temp(~mask_Y)/(4*beta_2);
    times(iter,1)=toc(tstart);
    
    %% Lx Update
    tstart = tic;
    Lx = soft_hosvd(L, Lam{2}, psi, 1/beta_2);
    times(iter,2) = toc(tstart);
    
    %% S Update
    tstart = tic;
    temp1 = beta_1*(Y-L-Nt+Lam{1});
    temp1(~mask_Y) = 0;
    temp2 = beta_4*(W+Lam{4});
    Sold = S;
    S = soft_threshold((temp1+temp2), lambda)./(beta_1+beta_4);
    S(~mask_Y) = soft_threshold(temp2(~mask_Y), lambda)./beta_4;
    times(iter,3) = toc(tstart);
    %% N update
    tstart = tic;
    Nt = (beta_1/(beta_1+alpha)).*(Y+Lam{1}-L-S);
    Nt(~mask_Y) = 0;
    times(iter,4) = toc(tstart);
    %% W Update
    tstart = tic;
    W = invD*(beta_4*Runfold(S-Lam{4})+beta_3*D'*Runfold(Z+Lam{3}));
    W = reshape(W, sz);
    times(iter,5) = toc(tstart);
    
    %% Z Update
    tstart = tic;
    Z = soft_threshold(mergeTensors(D, W, 2, 1)-Lam{3}, gamma/beta_3);
    times(iter,6) = toc(tstart);
    %% Dual Updates
    tstart = tic;
    temp = 0;
    for i=1:N
        Lam2_up = L-Lx{i};
        temp = temp + norm(Lam2_up(:))^2;
        Lam{2}{i} = Lam{2}{i}+Lam2_up;
    end
    Lam{1} = Lam{1} + Y-L-S-Nt;
    Lam{1}(~mask_Y) = 0;
    temp = sqrt(temp)/(sqrt(N)*norm(Y(:)));
    Lam{3} = Lam{3} - mergeTensors(D, W, 2, 1) + Z;
    Lam{4} = Lam{4} - S + W;
    times(iter,7) = toc(tstart);
    
    %% Error and objective calculations
%     obj_val(iter+1) = compute_obj(Y,Lx,S,Nt,W,Z,Lam1,Lam2,Lam3,D,lambda,gamma,alpha,beta_1,beta_2,beta_3);
%     err = abs(obj_val(iter)-obj_val(iter+1))/obj_val(iter);
    err = max(norm(S(:)-Sold(:))/norm(Sold(:)), temp);
    iter = iter+1;
    
    if err<=err_tol
        disp('Converged!')
        break;
    end
    if iter>=max_iter
        disp('Max iter')
        break;
    end
end

end

% function [val, term] = compute_obj(Y,L,S,Nt,W,Z,Lam1,Lam2,Lam3,D,lambda,gamma,alpha,beta_1,beta_2,beta_3)
% N = length(L);
% term = 0;
% lag1 = 0;
% for i=1:N
%     temp = beta_1/2*sum((Y-S-L{i}-Nt+Lam1{i}).^2,'all');
%     lag1 = lag1 + temp;
%     term = term + comp_nuclear(L{i}, i);
% end
% lag2 = (beta_2/2)*sum((mergeTensors(D, W, 2, 1)-Lam2-Z).^2,'all');
% lag3 = (beta_3/2)*sum((S-W-Lam3).^2,'all');
% term = term + lag1;
% term(2) = lambda*sum(abs(S),'all')+lag1+lag3;
% term(3) = (alpha/2)*norm(Nt(:))^2+lag1;
% term(4) = lag2+lag3;
% term(5) = gamma*sum(abs(Z),'all')+lag2;
% val = sum(term)-2*lag1-lag3-lag2;
% end