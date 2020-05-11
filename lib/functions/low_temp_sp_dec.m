function [L, S] = low_temp_sp_dec(Y, param)
% [S] = low_temp_sp_dec(Y, param)
% Low rank+Temporally Smooth Sparse Decomposition

N = ndims(Y);
sz = size(Y);
max_iter = param.max_iter;
err_tol = param.err_tol;
lambda = param.lambda;
gamma = param.gamma;
beta_1 = param.beta_1;
beta_2 = param.beta_2;
beta_3 = param.beta_3;
L = cell(1, N);
for i=1:N
    L{i} = zeros(size(Y));
end
S = zeros(size(Y));
W = zeros(size(Y));
Z = zeros(size(Y));
D = convmtx([1,-1], size(Y,1));
D(:,end) = [];
D(end,1) = -1;
Lam1 = cell(1, N);
for i=1:N
    Lam1{i} = zeros(size(Y));
end
Lam2 = zeros(size(Y));
Lam3 = zeros(size(Y));

err = inf;
iter = 1;
obj_val = compute_obj(Y,L,S,W,Z,Lam1,Lam2,Lam3,D,lambda,gamma,beta_1,beta_2,beta_3);
while true
    %% L Update
    L = soft_hosvd(Y-S, Lam1, 1/beta_1);
    
    %% S Update
    temp1 = zeros(size(Y));
    for i=1:N
        temp1 = temp1+beta_1*(Y-L{i}+Lam1{i});
    end
    temp2 = beta_3*(W+Lam3);
    Sold = S;
    S = soft_threshold((temp1+temp2), lambda)./(N*beta_1+beta_3);
    
    %% W Update
    Dtemp = D'*D;
    Dtemp2 = D';
    W = (beta_3*eye(sz(1))+beta_2*Dtemp)^-1 ...
        *(beta_3*Runfold(S-Lam3)+beta_2*Dtemp2*Runfold(Z+Lam2));
    W = reshape(W, sz);
    
    %% Z Update
    Z = soft_threshold(mergeTensors(D, W, 2, 1)-Lam2, gamma/beta_2);
    
    %% Dual Updates
    temp = 0;
    for i=1:N
        Lam1_up = Y-L{i}-S;
        temp = temp + norm(Lam1_up(:))/(sqrt(N)*norm(Y(:)));
        Lam1{i} = Lam1{i}+Lam1_up;
    end
    Lam2 = Lam2-mergeTensors(D, W, 2, 1)+Z;
    Lam3 = Lam3-S+W;
    
    %% Error and objective calculations
%     obj_val(iter+1) = compute_obj(Y,L,S,W,Z,Lam1,Lam2,Lam3,D,lambda,gamma,beta_1,beta_2,beta_3);
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
temp = zeros(size(Y));
for i=1:N
    temp = temp+L{i};
end
L = temp/N;

end

function val = compute_obj(Y,L,S,W,Z,Lam1,Lam2,Lam3,D,lambda,gamma,beta_1,beta_2,beta_3)
N = length(L);
term = 0;
for i=1:N
    term = term + comp_nuclear(L{i}, i)+beta_1/2*sum((Y-S-L{i}-Lam1{i}).^2,'all');
end
term(2) = lambda*sum(abs(S),'all');
term(3) = gamma*sum(abs(Z),'all');
term(4) = beta_2/2*sum((mergeTensors(D, W, 2, 1)-Lam2-Z).^2,'all') +...
    beta_3/2*sum((S-W-Lam3).^2,'all');
val = sum(term);
end