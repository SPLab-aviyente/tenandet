function [L, S, Nt] = logss(Y, param)
% [L, S, Nt] = logss(Y, param)
% Low rank On Graphs plus Smooth-Sparse Decomposition

N = ndims(Y);
sz = size(Y);
max_iter = param.max_iter;
err_tol = param.err_tol;
alpha = param.alpha;
theta = param.theta;
lambda = param.lambda;
gamma = param.gamma;
psi = param.psi;
beta_1 = param.beta_1;
beta_2 = param.beta_2;
beta_3 = param.beta_3;
% beta_4 = param.beta_4;
L = cell(1, N);
for i=1:N
    L{i} = zeros(size(Y));
end
% La = L;
S = zeros(size(Y));
Nt = zeros(size(Y));
W = zeros(size(Y));
Z = zeros(size(Y));
Phi = cell(1, N);
U = Phi;
Sig = Phi;
flags = [0,0,0,1];
for i = 1:N
    mods = setdiff(1:N, i);
    Phi{i} = get_graphL(permute(Y, [mods, i]), 10, flags(i));
    [U{i}, Sig{i}, ~] = svdtrunc2(pinv(Phi{i}), 0.01);
end
D = convmtx([1,-1], size(Y,1));
D(:,end) = [];
D(end,1) = -1;
if beta_3~=0 && beta_2~=0
    invD = (beta_3*eye(sz(1))+beta_2*(D'*D))^-1;
else
    invD = zeros(sz(1));
end
Lam{1} = cell(1, N);
for i=1:N
    Lam{1}{i} = zeros(size(Y));
end
% Lam{4} = Lam{1};
Lam{2} = zeros(size(Y));
Lam{3} = zeros(size(Y));

timeL = [];
timeS = []; 
timeN = [];
timeZ = [];
timeW = [];
timeDual = [];
iter = 1;
while true
    %% L Update
    tstart = tic;
    [L, ~] = soft_hosvd_gc(Y-S-Nt, Lam{1}, Sig, U, psi, beta_1);
    timeL(end+1)=toc(tstart);
%     %% La Update
%     La = graph_reg_update(L,Lam{4},inv_Phi);
    %% S Update
    tstart = tic;
    temp1 = zeros(size(Y));
    for i=1:N
        temp1 = temp1+beta_1*(Y-L{i}-Nt+Lam{1}{i});
    end
    temp2 = beta_3*(W+Lam{3});
    Sold = S;
    S = soft_threshold((temp1+temp2), lambda)./(N*beta_1+beta_3);
    timeS(end+1)=toc(tstart);
    %% N update
    tstart = tic;
    Nt = zeros(size(Y));
    for i=1:N
        Nt = Nt+(beta_1/(N*beta_1+alpha)).*(Y+Lam{1}{i}-L{i}-S);
    end
    timeN(end+1)=toc(tstart);
    %% W Update
    tstart = tic;
    W = invD*(beta_3*Runfold(S-Lam{3})+beta_2*D'*Runfold(Z+Lam{2}));
    W = reshape(W, sz);
    timeW(end+1)=toc(tstart);
    
    %% Z Update
    tstart = tic;
    Z = soft_threshold(mergeTensors(D, W, 2, 1)-Lam{2}, gamma/(beta_2+eps));
    timeZ(end+1)=toc(tstart);
    %% Dual Updates
    tstart = tic;
    temp = 0;
    for i=1:N
        Lam1_up = Y-L{i}-S-Nt;
        temp = temp + norm(Lam1_up(:))^2;
        Lam{1}{i} = Lam{1}{i}+Lam1_up;
%         Lam{4}{i} = Lam{4}{i}-(L{i}-La{i});
    end
    temp = sqrt(temp)/(sqrt(N)*norm(Y(:)));
    Lam{2} = Lam{2}-mergeTensors(D, W, 2, 1)+Z;
    Lam{3} = Lam{3}-S+W;
    timeDual(end+1) = toc(tstart);
    
    %% Error calculations
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