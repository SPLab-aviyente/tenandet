clear
nodisp = true;
sizes = [24,7,52,25];
Y = [];
%%Uncomment for synthesizing using NYC data
% load nyc_tensors.mat
% load regions.mat
% arrs = squeeze(sum(reshape(arrs,6,24,365,[]),1));
% Y = double(reshape(arrs(:,1:364,regions), 24, 7, 52,[]));
% sizes=size(Y);
%%_________________________________________________________________________
%%Uncomment for synthesizing using traffic data
% [~, Y, ~] = get_traffic_data;
param.err_tol = 0.01;
anom_list = 700;
n_missing = round([0.01].*(prod(sizes)/24));%,0.05,0.1,0.2:.2:0.8
for ind_outer=1:length(n_missing)
    num_anom = anom_list(1);
    len_anom = 7;
    amp_anom = 1.5;
    [X, Y_gen, Yn, ind_removed, mat_anom] = gendata(sizes, num_anom, len_anom, amp_anom, Y, n_missing(ind_outer), false);
    [mah_Yn, L_ee] = mahal_dist(Yn);
    rmse_ee(ind_outer) = norm(L_ee(:)-Y_gen(:))/sqrt(numel(Y_gen));
    mape_ee(ind_outer) = sum(abs(L_ee-Y_gen),'all')/numel(Y_gen);
    if ind_outer == length(anom_list)
%         [fpr, recall_e] = analyze_envelope(Yn, X, ind_removed);
    end
    [k_list(:,ind_outer), precision_or(:,ind_outer), recall_or(:,ind_outer), fpr_or(:,ind_outer)] = analyze_top_K(mah_Yn, X, ind_removed);
    
    param.psi = [.01,1,3,.001];
    param.lambda = 1/4/sqrt(max(size(Yn)));
    param.gamma = 1/4/sqrt(max(size(Yn)));
    gloss_subscript
    lrtssd_subscript
%     Yn_d = S_lrt;
%     precision_wlrs = precision_lrs;
%     recall_wlrs = recall_lrs;
%     fpr_wlrs = fpr_lrs;
    horpca_subscript
    rmse_whorpca(ind_outer) = rmse_horpca(ind_outer);
    mape_whorpca(ind_outer) = mape_horpca(ind_outer);
    precision_whorpca(:,ind_outer) = precision_horpca(:,ind_outer);
    recall_whorpca(:,ind_outer) = recall_horpca(:,ind_outer);
    fpr_whorpca(:,ind_outer) = fpr_horpca(:,ind_outer);
    param.psi = [1,1,1,1];
    lof_subscript
%     param.gamma = 1/5/(max(size(Y)));
    horpca_subscript
end

% figure,
% plot(fpr, recall_e,'DisplayName','Original')
% hold on;
% plot(fpr_lof,recall_lof_en,'DisplayName','LOF','LineWidth',3)
% plot(fpr_lrs,recall_lrs_en,'DisplayName','LRTSSD','LineWidth',3)
% plot(fpr_horpca_en,recall_horpca_en,'DisplayName','HORPCA','LineWidth',3)
% plot(fpr_ghorpca_en,recall_ghorpca_en,'DisplayName','GHORPCA','LineWidth',3)
% plot(fpr_gloss_en,recall_gloss_en,'DisplayName','GLOSS','LineWidth',3)
% legend, grid
% title('ROC of the envelope')

% k_list = 24*100*n_missing/numel(Yn);
figure,
plot(k_list, precision_or,'DisplayName','EE','LineWidth',3);
hold on;
plot(k_list, precision_lof,'DisplayName','LOF','LineWidth',3);
plot(k_list, precision_horpca,'DisplayName','HORPCA','LineWidth',3);
plot(k_list, precision_whorpca,'DisplayName','WHORPCA','LineWidth',3);
plot(k_list, precision_lrs,'DisplayName','LOSS','LineWidth',3);
% plot(k_list, precision_lrs_db,'DisplayName','LOSS-DB','LineWidth',3);
plot(k_list, precision_gloss,'DisplayName','GLOSS','LineWidth',3)
% plot(k_list, precision_gloss_2,'DisplayName','GLOSS-2','LineWidth',3)
% plot(k_list, precision_dbscan,'DisplayName','DBSCAN','LineWidth',3)
% plot(k_list, precision_wlrs,'DisplayName','WLOSS','LineWidth',3);
legend, grid
title('Precision')
ax = gcf;
ax.CurrentAxes.FontSize=19;
ax.CurrentAxes.FontWeight = 'bold';

figure,
plot(k_list, recall_or,'DisplayName','EE','LineWidth',3);
hold on;
plot(k_list, recall_lof,'DisplayName','LOF','LineWidth',3);
plot(k_list, recall_horpca,'DisplayName','HORPCA','LineWidth',3);
plot(k_list, recall_whorpca,'DisplayName','WHORPCA','LineWidth',3);
plot(k_list, recall_lrs,'DisplayName','LOSS','LineWidth',3);
% plot(k_list, recall_lrs_db,'DisplayName','LOSS-DB','LineWidth',3);
plot(k_list, recall_gloss,'DisplayName','GLOSS','LineWidth',3)
% plot(k_list, recall_gloss_2,'DisplayName','GLOSS-2','LineWidth',3)
% plot(k_list, recall_dbscan,'DisplayName','DBSCAN','LineWidth',3)
% plot(k_list, recall_wlrs,'DisplayName','WLOSS','LineWidth',3);
legend, grid
title('Recall')
ax = gcf;
ax.CurrentAxes.FontSize=19;
ax.CurrentAxes.FontWeight = 'bold';

figure,
plot(fpr_or, recall_or,'DisplayName','EE','LineWidth',3);
hold on;
plot(fpr_lof, recall_lof,'DisplayName','LOF','LineWidth',3);
plot(fpr_horpca, recall_horpca,'DisplayName','HORPCA','LineWidth',3);
plot(fpr_whorpca, recall_whorpca,'DisplayName','WHORPCA','LineWidth',3);
plot(fpr_lrs, recall_lrs,'DisplayName','LOSS','LineWidth',3);
% plot(fpr_lrs_db, recall_lrs_db,'DisplayName','LOSS-DB','LineWidth',3);
plot(fpr_gloss,recall_gloss,'DisplayName','GLOSS','LineWidth',3)
% plot(fpr_gloss_2,recall_gloss_2,'DisplayName','GLOSS-2','LineWidth',3)
% plot(fpr_dbscan, recall_dbscan,'DisplayName','DBSCAN','LineWidth',3)
% plot(fpr_wlrs, recall_wlrs,'DisplayName','WLOSS','LineWidth',3);
legend, grid
title('ROC')
ax = gcf;
ax.CurrentAxes.FontSize=19;
ax.CurrentAxes.FontWeight = 'bold';