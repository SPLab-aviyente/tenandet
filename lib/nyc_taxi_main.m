clear
sizes = [144,7,52,10];
load nyc_tensors.mat
load regions.mat
arrs = squeeze(sum(reshape(arrs,6,24,365,[]),1));
Y = double(reshape(arrs(:,1:364,regions), 24, 7, 52,[]));
Y(:,1,53,:) = arrs(:,365,regions);
Y(:,2:7,53,:) = mean(Y(:,2:7,1:52,:),3);
ind_removed = [];

param.lambda = 1/sqrt(max(size(Y)));
param.gamma = 1/sqrt(max(size(Y)));
param.theta = 1/sum(std(t2m(Y,4),[],2).^2);
param.alpha = 0;
param.psi = [.01,1,10,.001];
param.beta_1 = 1/(5*std(Y(:)));
param.beta_2 = param.beta_1;
param.beta_3 = param.beta_1;
param.beta_4 = param.beta_1;
param.max_iter = 100;
param.err_tol = 0.01;
[~, S_wlrt, ~] = low_temp_sp_dec(Y, param);
[~,S_whorpca, ~] = horpca(Y, param);
param.psi = [1,1,1,1];
[~, S_lrt, ~] = low_temp_sp_dec(Y, param);
[~,S_horpca, ~] = horpca(Y, param);
a_L = apply_lof(Y, 10);

list_k = [20,100,500,1000,2000,5000,7000,14000,21000];
det_gloss = zeros(size(list_k));
det_glos = det_gloss;
det_lrt = det_gloss;
det_whorpca = det_gloss;
det_horpca = det_gloss;
det_ee = det_gloss;
det_lof = det_gloss;
for i=1:length(list_k)
    det_gloss(i) = sum(detect_real_events(mahal_dist(S_wlrt), ind_removed, list_k(i)));
    det_glos(i) = sum(detect_real_events(mahal_dist(S_gh), ind_removed, list_k(i)));
    det_lrt(i) = sum(detect_real_events(mahal_dist(S_lrt), ind_removed, list_k(i)));
    det_whorpca(i) = sum(detect_real_events(mahal_dist(S_whorpca), ind_removed, list_k(i)));
    det_horpca(i) = sum(detect_real_events(mahal_dist(S_horpca), ind_removed, list_k(i)));
    det_ee(i) = sum(detect_real_events(mahal_dist(Y), ind_removed, list_k(i)));
    det_lof(i) = sum(detect_real_events(a_L, ind_removed, list_k(i)));
end

figure,
plot(100*list_k/numel(Y), det_gloss,'DisplayName','WLOSS','LineWidth',3)
hold
plot(100*list_k/numel(Y), det_lrt,'DisplayName','LOSS','LineWidth',3)
plot(100*list_k/numel(Y), det_whorpca,'DisplayName','WHORPCA','LineWidth',3)
plot(100*list_k/numel(Y), det_horpca,'DisplayName','HORPCA','LineWidth',3)
plot(100*list_k/numel(Y), det_ee,'DisplayName','EE','LineWidth',3)
plot(100*list_k/numel(Y), det_lof,'DisplayName','LOF','LineWidth',3)
legend
grid
ax = gcf;
ax.CurrentAxes.FontSize=19;
ax.CurrentAxes.FontWeight = 'bold';
% plot_nyc(Y, S, S_gh, S_lrt, S_horpca, a_L, 1, 75)

