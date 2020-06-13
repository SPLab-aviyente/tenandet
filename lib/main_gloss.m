clear
nodisp = true;
sizes = [144,7,52,10];
load nyc_tensors.mat
load regions.mat
arrs = squeeze(sum(reshape(arrs,6,24,365,[]),1));
Y = double(reshape(arrs(:,1:364,regions), 24, 7, 52,[]));
% [~, Y, ~] = get_traffic_data;
anom_list = 50;
n_missing = [100,200,500,1000,5000,10000];
for ind_outer=1:length(n_missing)
    num_anom = anom_list(1);
    len_anom = 7;
    amp_anom = 1.5;
    [X, ~, Yn, ind_removed, mat_anom] = gendata(sizes, num_anom, len_anom, amp_anom, Y, n_missing(ind_outer));
    if ind_outer == length(anom_list)
%         [fpr, recall_e] = analyze_envelope(Yn, X, ind_removed);
    end
    [k_list(:,ind_outer), precision_or(:,ind_outer), recall_or(:,ind_outer), fpr_or(:,ind_outer)] = analyze_top_K(mahal_dist(Yn), X, ind_removed);
    
    param.psi = [.01,5,20,.001];
    param.gamma = 1/sqrt(max(size(Y)));
    gloss_subscript
    lrtssd_subscript
%     precision_wlrs = precision_lrs;
%     recall_wlrs = recall_lrs;
%     fpr_wlrs = fpr_lrs;
    horpca_subscript
    lof_subscript
    precision_whorpca(:,ind_outer) = precision_horpca(:,ind_outer);
    recall_whorpca(:,ind_outer) = recall_horpca(:,ind_outer);
    fpr_whorpca(:,ind_outer) = fpr_horpca(:,ind_outer);
    param.psi = [1,1,1,1];
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

k_list = n_missing;
figure,
plot(k_list, precision_or,'DisplayName','EE','LineWidth',3);
hold on;
plot(k_list, precision_lof,'DisplayName','LOF','LineWidth',3);
plot(k_list, precision_horpca,'DisplayName','HORPCA','LineWidth',3);
plot(k_list, precision_whorpca,'DisplayName','WHORPCA','LineWidth',3);
plot(k_list, precision_lrs,'DisplayName','LOSS','LineWidth',3);
plot(k_list, precision_gloss,'DisplayName','GLOSS','LineWidth',3)
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
plot(k_list, recall_gloss,'DisplayName','GLOSS','LineWidth',3)
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
plot(fpr_gloss,recall_gloss,'DisplayName','GLOSS','LineWidth',3)
% plot(fpr_wlrs, recall_wlrs,'DisplayName','WLOSS','LineWidth',3);
legend, grid
title('ROC')
ax = gcf;
ax.CurrentAxes.FontSize=19;
ax.CurrentAxes.FontWeight = 'bold';