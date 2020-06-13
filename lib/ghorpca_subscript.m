
param.lambda = 1/sqrt(max(size(Y)));
param.theta = 1/(norm(Y(:)));
param.alpha = 0;
param.psi = [.01,1,10,.001];
param.beta_1 = 1/(5*std(Y(:)));
param.beta_4 = param.beta_1/1;
param.max_iter = 100;
param.err_tol = 0.01;
tic;
[L,S_gh, N, obj_val] = ghorpca(Yn, param);
time_glos = toc
%% Envelope Analysis
if ind_outer == length(anom_list)
%     [fpr_ghorpca_en, recall_ghorpca_en] = analyze_envelope(S_gh, X, ind_removed);
end
%% Top-K Analysis
mahal_S = mahal_dist(S_gh);
[~, precision_ghorpca, recall_ghorpca, fpr_ghorpca] = analyze_top_K(mahal_S, X, ind_removed);
%% Visualize Decomposition
% plot_sensor_new(X, Yn, S, L, 5)
% plot_sensor_new(permute(X,[1,3,2,4]), permute(Yn,[1,3,2,4]), permute(S,[1,3,2,4]), permute(L,[1,3,2,4]), 5)
if ~nodisp
    plot_sensor_new(permute(X,[1,3,4,2]), permute(Yn,[1,3,4,2]), permute(mahal_S,[1,3,4,2]), permute(S_gh,[1,3,4,2]), 3,40)
end

if ind_outer == length(anom_list)
%     plot_sensor_new(permute(X,[1,3,2,4]), permute(Yn,[1,3,2,4]), permute(S,[1,3,2,4]), permute(L,[1,3,2,4]), 5)
end

if ~nodisp
    a_G = zeros(size(Y));
    for i=1:size(Yn,1)
        for j=1:size(Yn,2)
            for k=1:size(Yn,4)
                [~,a_G(i,j,:,k)] = LOF(squeeze(S_gh(i,j,:,k)),1);
            end
        end
    end
    plot_sensor_new(permute(X,[1,3,4,2]), permute(Yn,[1,3,4,2]), permute(a_G,[1,3,4,2]), permute(S_gh,[1,3,4,2]), 3,10)

    figure,
    plot(fpr_ghorpca_en,recall_ghorpca_en,'DisplayName','GHORPCA')
    legend, grid
    title('ROC of the envelope')

    figure,
    plot(k_list, precision_ghorpca,'DisplayName','GHORPCA');
    legend, grid
    title('Precision')

    figure,
    plot(k_list, recall_ghorpca,'DisplayName','GHORPCA');
    legend, grid
    title('Recall')
end