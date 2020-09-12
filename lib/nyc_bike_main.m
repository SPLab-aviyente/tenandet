clear
load nyc_bikedata.mat
load regions.mat
arrs = squeeze(sum(reshape(arrs,6,24,365,[]),1));
Y = double(reshape(arrs(:,1:364,:), 24, 7, 52,[]));
Y(:,1,53,:) = arrs(:,365,:);
Y(:,2:7,53,:) = mean(Y(:,2:7,1:52,:),3);
ind_removed = [];
for i=1:4
    std_m(i) = sum(std(t2m(Y,i)));
end
list_k = [20,100,500,1000,2000,5000,7000,14000,21000,30000,50000,70000,90000];
lam_list = 10.^[3];
gam_list = 10.^[3];
len_g = length(gam_list);
len_k = length(list_k);
gloss_k = zeros(length(lam_list),len_g,len_k);
% logss_k = gloss_k;
loss_k = gloss_k;
% whorpca_k = gloss_k;
det_gloss = zeros(length(lam_list),len_g,len_k,20);
det_logss = det_gloss;
det_lrt = det_gloss;
det_whorpca = det_lrt;
for i=1:length(lam_list)
    param = [];
    param.ind_m = ind_removed;
    param.lambda = 1/sqrt(max(size(Y)));
    param.theta = 1/(sum(std(t2m(Y,4),[],2).^2))^2;
    param.alpha = 0;
    param.psi = max(std_m)*std_m.^-1;%[0.1, 1,5,0.01];
    param.beta_1 = 1/(5*std(Y(:)));
    param.beta_2 = param.beta_1;
    param.beta_3 = param.beta_1;
    param.beta_4 = param.beta_1;
    param.beta_5 = param.beta_1;
    param.max_iter = 100;
    param.err_tol = 0.01;
    [L_whorpca,S_whorpca, ~] = horpca(Y, param);
    param.psi = [1,1,1,1];
    [~,S_horpca, ~] = horpca(Y, param);
%     m_S_whorpca = mahal_dist(S_whorpca);
    a_L = apply_lof(Y, 10);
    out_fr = 0.1;
    OCS = one_class_svm(Y, out_fr);
    for j=1:len_g
        param.psi = max(std_m)*std_m.^-1;
        param.lambda = 1/sqrt(max(size(Y)));
        param.gamma = 1/sqrt(max(size(Y)));
        [L_lrt, S_lrt, ~] = low_temp_sp_dec(Y, param);
%         param.psi = min(std_m)*std_m.^-2;
%         [L_logss,S_logss,~] = logss(Y,param);
%         param.psi = [0.1, 1,5,0.01];%min(std_m)*std_m.^-1;
        param.lambda = 1/sqrt(numel(Y));
        param.gamma = 1/sqrt(numel(Y));
        [L_gloss, S_gloss,~] = gloss_3(Y, param);

        results_real_nyc(S_gloss, S_lrt, S_whorpca, S_horpca, Y, a_L, OCS, ind_removed, param, true)
%         m_S_gloss = mahal_dist(S_gloss);
%         o_S_gloss = one_class_svm(S_gloss, out_fr);
%         lof_S_gl = apply_lof(S_gloss, 10);
%         m_S_logss = mahal_dist(S_logss);
%         m_S_lrt = mahal_dist(S_lrt);
%         m_S_horpca = mahal_dist(S_horpca);
%         for k=1:len_k
%             [det_gloss(i,j,k,:), gloss_k(i,j,k)] = (detect_real_events((m_S_gloss), ind_removed, list_k(k)));
%             [det_o_gl(i,j,k,:), gloss_k(i,j,k)] = (detect_real_events((o_S_gloss), ind_removed, list_k(k)));
%             [det_lo_gl(i,j,k,:), gloss_k(i,j,k)] = (detect_real_events((lof_S_gl), ind_removed, list_k(k)));
%             [det_logss(i,j,k,:), logss_k(i,j,k)] = (detect_real_events((m_S_logss), ind_removed, list_k(k)));
%             [det_lrt(i,j,k,:), loss_k(i,j,k)] = (detect_real_events((m_S_lrt), ind_removed, list_k(k)));
%             [det_whorpca(i,j,k,:), whorpca_k(i,j,k)] = (detect_real_events((m_S_whorpca), ind_removed, list_k(k)));
%             [det_horpca(i,j,k,:), horpca_k(i,j,k)] = (detect_real_events((m_S_horpca), ind_removed, list_k(k)));
%             [det_lof(i,j,k,:), horpca_k(i,j,k)] = (detect_real_events((m_S_horpca), ind_removed, list_k(k)));
%         end
    end
end


