clear
load nyc_bikedata.mat
load regions.mat
arrs = squeeze(sum(reshape(arrs,6,24,365,[]),1));
Y = double(reshape(arrs(:,1:364,:), 24, 7, 52,[]));
ind_removed = [];

rng(1);
param.lambda = 1/8/sqrt(max(size(Y)));
param.gamma = 1/16/sqrt(max(size(Y)));
param.theta = 1/(sum(std(t2m(Y,4),[],2).^2));
param.alpha = 0;
param.psi = [.01,1,10,.001];
param.beta_1 = 1/(5*std(Y(:)));
param.beta_2 = param.beta_1;
param.beta_3 = param.beta_1;
param.beta_4 = param.beta_1;
param.beta_5 = param.beta_1;
param.max_iter = 100;
param.err_tol = 0.01;
[L_gloss,S_gloss,~] = gloss(Y, param);
[L, S_lrt, ~] = low_temp_sp_dec(Y, param);
param.lambda = 1/sqrt(max(size(Y)));
param.gamma = 1/sqrt(max(size(Y)));
[~,S_whorpca, ~] = horpca(Y, param);
param.psi = [1,1,1,1];
[~,S_horpca, ~] = horpca(Y, param);
param.lambda = 1/15/(max(size(Y)));
param.gamma = 1/80/(max(size(Y)));
param.theta = 1/(sum(std(t2m(Y,4),[],2).^2));
param.alpha = 0;
param.psi = [.01,.1,10,.001];
[L_gloss_3,S_gloss_3,~] = gloss_3(Y, param);
a_L = apply_lof(Y, 10);

results_real_nyc(S_gloss, S_gloss_3, S_lrt, S_whorpca, S_horpca, Y, a_L, ind_removed, param, true)

