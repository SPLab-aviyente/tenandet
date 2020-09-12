clear
nodisp = true;
sizes = [24,7,52,25];
Y = [];
%%Uncomment for synthesizing using NYC data
load nyc_tensors.mat
% load nyc_bikedata.mat
load regions.mat
arrs = squeeze(sum(reshape(arrs,6,24,365,[]),1));
Y = double(reshape(arrs(:,1:364,regions), 24, 7, 52,[]));
% Y = double(reshape(arrs(:,1:364,1:81), 24, 7, 52,[]));
sizes = size(Y);
%%_________________________________________________________________________
%%Uncomment for synthesizing using traffic data
% [~, Y, ~] = get_traffic_data;
anom_list = 700;
amp_anom = [2.5];
% load(['preliminary_results\gloss_amp', num2str(amp_anom),'.mat'])
n_missing = round([0.2,0.4,0.6].*(prod(sizes)/24));%,0.05,0.1,0.2:.2:0.8
psi_list = 10.^(0);
lam_list_1 = (10.^[-4,0])*1/(max(size(Y)));
lam_list_2 = (10.^[-1,0])*1/(numel(Y));
gam_list_1 = (10.^[0])*1/(max(size(Y)));
gam_list_2 = (10.^[2])*1/(numel(Y));
len_lam = length(lam_list_1);
len_gam = length(gam_list_1);
len_n = length(n_missing);
len_a = length(amp_anom);
precision = cell(1,len_n);
recall = precision;
fpr = recall;

geo_mean_std = geo_mean([1/(sum(std(t2m(Y,4),[],2).^2)),...
    1/(sum(std(t2m(Y,3),[],2).^2)),...
    1/(sum(std(t2m(Y,2),[],2).^2)),...
    1/(sum(std(t2m(Y,1),[],2).^2))]);
for i=1:4
    std_m(i) = sum(std(t2m(Y,i)));
end
parfor ind_n=1:len_n
    for ind_amp = 1:len_a
    num_anom = anom_list(1);
    len_anom = 7;
    [X, Y_gen, Yn, ind_removed, mat_anom] = gendata(sizes, num_anom, len_anom, amp_anom(ind_amp), Y, n_missing(ind_n), false);
    [mah_Yn, L_ee] = mahal_dist(Yn);
    [~, precision{ind_n}(:, 1, ind_amp, 1, 1), recall{ind_n}(:, 1, ind_amp, 1,  1), fpr{ind_n}(:, 1, ind_amp, 1, 1)] = analyze_top_K(mah_Yn, X, ind_removed);
    [~, precision{ind_n}(:, 1, ind_amp, 1, 3), recall{ind_n}(:, 1, ind_amp, 1, 3),fpr{ind_n}(:, 1, ind_amp, 1, 3)] = ocsvm_subscript(Yn, X, param);
    [~, precision{ind_n}(:, 1, ind_amp, 1, 2), recall{ind_n}(:, 1, ind_amp, 1, 2),fpr{ind_n}(:, 1, ind_amp, 1, 2)] = lof_subscript(Yn, X, param);
    param =[];
    param.err_tol = 0.01;
    param.ind_m = ind_removed;
    param.alpha = 0;
    param.psi = [1,1,1,1];
    param.lambda = 1/sqrt(max(size(Yn)));
    param.beta_1 = 1/(5*std(Yn(setdiff(1:numel(Yn), ind_removed))));
    [L_h, S_h, precision{ind_n}(:, 1, ind_amp, 1, 4), recall{ind_n}(:, 1, ind_amp, 1, 4),fpr{ind_n}(:, 1, ind_amp, 1, 4)] = horpca_subscript(Yn, Y_gen, X, param);
        for ind_lam=1:len_lam
            [L_wh, S_wh, precision{ind_n}(:, 1, ind_amp, ind_lam, 5), recall{ind_n}(:, 1, ind_amp, ind_lam, 5),fpr{ind_n}(:, 1, ind_amp, ind_lam, 5)] = horpca_subscript(Yn, Y_gen, X, param);
            param.theta = geo_mean_std;
            param.beta_2 = param.beta_1;
            param.beta_3 = param.beta_1;
            param.beta_4 = param.beta_1;
            param.beta_5 = param.beta_1;
            param.max_iter = 100;
            for ind_gam =1:len_gam
                param.lambda = lam_list_1(ind_lam);
                param.psi = psi_list(1)*min(std_m)*std_m.^-1;
                param.gamma = gam_list_1(ind_gam);
                [L_loss,S_loss,precision{ind_n}(:, ind_gam, ind_amp, ind_lam, 6), recall{ind_n}(:, ind_gam, ind_amp, ind_lam, 6), fpr{ind_n}(:, ind_gam, ind_amp, ind_lam, 6)] = lrtssd_subscript(Yn, Y_gen, X, param);
                
                param.lambda = lam_list_2(ind_lam);
                param.gamma = gam_list_2(ind_gam);
                param.psi = psi_list(1)*min(std_m)*std_m.^-1;
                [L_g, S_g, precision{ind_n}(:, ind_gam, ind_amp, ind_lam, 7:9), recall{ind_n}(:, ind_gam, ind_amp, ind_lam, 7:9), fpr{ind_n}(:, ind_gam, ind_amp, ind_lam, 7:9)] = gloss_3_subscript(Yn, Y_gen, X, param);
    %             [~, ~, precision{ind_outer}(:, ind_inner, k,10), recall{ind_outer}(:, ind_inner, k,10), fpr{ind_outer}(:, ind_inner, k,10)] = logss_subscript(Yn, Y_gen, X, param);
            end
        end
    end
end

auc = zeros(len_lam, len_gam, len_a, len_n, 9);
for l = 1:len_n
    for i=1:len_lam
        for j=1:len_gam
            for k=1:len_a
                figure,
                plot(fpr{l}(:,1,k,1,1), recall{l}(:,1,k,1,1),'Color','#0072BD','DisplayName','EE','LineWidth',7);
                hold on;
                plot(fpr{l}(:,1,k,1,3), recall{l}(:,1,k,1,3),'Color','#EDB120','DisplayName','OCSVM','LineWidth',7)
                plot(fpr{l}(:,1,k,1,2), recall{l}(:,1,k,1,2),'Color','#D95319','DisplayName','LOF','LineWidth',7);
                plot(fpr{l}(:,1,k,1,4), recall{l}(:,1,k,1,4),'Color','#7E2F8E','DisplayName','HORPCA','LineWidth',7);
                plot(fpr{l}(:,1,k,i,5), recall{l}(:,1,k,i,5),'--','Color','#77AC30','DisplayName','WHORPCA','LineWidth',7);
                plot(fpr{l}(:,j,k,i,6), recall{l}(:,j,k,i,6),'--+','Color','#4DBEEE','DisplayName','LOSS','LineWidth',7,'MarkerSize',20);
                plot(fpr{l}(:,j,k,i,7), recall{l}(:,j,k,i,7),'--x','Color','#0072BD','DisplayName','GLOSS-EE','LineWidth',7,'MarkerSize',20)
                plot(fpr{l}(:,j,k,i,8), recall{l}(:,j,k,i,8),'--x','Color','#EDB120','DisplayName','GLOSS-SVM','LineWidth',7,'MarkerSize',20)
                plot(fpr{l}(:,j,k,i,9), recall{l}(:,j,k,i,9),'--x','Color','#D95319','DisplayName','GLOSS-LOF','LineWidth',7,'MarkerSize',20)
                %         plot(fpr{10}(:,i), recall{10}(:,i),'--x','Color','#0072BD','DisplayName','LOGSS','LineWidth',7,'MarkerSize',20)
                % plot(fpr_dbscan, recall_dbscan,'DisplayName','DBSCAN','LineWidth',3)
                % plot(fpr_wlrs, recall_wlrs,'DisplayName','WLOSS','LineWidth',3);
                legend('location','southeast'), grid
                title(['ROC ', num2str(round(100*24*n_missing(l)/prod(sizes))), '% miss, ', num2str(amp_anom(k)),' amp ', ' 0.1 std'])
                ax = gcf;
                ax.CurrentAxes.FontSize = 26;
                ax.CurrentAxes.FontWeight = 'bold';
                auc(i,j,k,l,1) = trapz(fpr{l}(:,1,k,1,1), recall{l}(:,1,k,1,1));
                auc(i,j,k,l,2) = trapz(fpr{l}(:,1,k,1,3), recall{l}(:,1,k,1,3));
                auc(i,j,k,l,3) = trapz(fpr{l}(:,1,k,1,2), recall{l}(:,1,k,1,2));
                auc(i,j,k,l,4) = trapz(fpr{l}(:,1,k,1,4), recall{l}(:,1,k,1,4));
                auc(i,j,k,l,5) = trapz(fpr{l}(:,1,k,i,5), recall{l}(:,1,k,i,5));
                auc(i,j,k,l,6) = trapz(fpr{l}(:,j,k,i,6), recall{l}(:,j,k,i,6));
                auc(i,j,k,l,7) = trapz(fpr{l}(:,j,k,i,7), recall{l}(:,j,k,i,7));
                auc(i,j,k,l,8) = trapz(fpr{l}(:,j,k,i,8), recall{l}(:,j,k,i,8));
                auc(i,j,k,l,9) = trapz(fpr{l}(:,j,k,i,9), recall{l}(:,j,k,i,9));
            end
        end
    end
end

save(['preliminary_results\gloss_prctgs.mat'], 'recall', 'fpr', 'gam_list_1', 'gam_list_2', 'lam_list_1', 'lam_list_2', 'amp_anom', 'n_missing', 'auc')