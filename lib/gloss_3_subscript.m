
param.theta = 1/(sum(std(t2m(Yn,4),[],2).^2));
param.alpha = 0;
param.beta_1 = 1/(5*std(Yn(setdiff(1:numel(Yn), ind_removed))));
param.beta_2 = param.beta_1;
param.beta_3 = param.beta_1;
param.beta_4 = param.beta_1;
param.beta_5 = param.beta_1;
param.max_iter = 100;

[L_3,S_3,~] = gloss_3(Yn, param);
rmse_gloss_3(ind_outer) = norm(L_3(:)-Y_gen(:))/sqrt(numel(Y_gen));
mape_gloss_3(ind_outer) = sum(abs(L_3-Y_gen),'all')/numel(Y_gen);
time_gloss_3 = toc

out_fr = 0.7;
S_svm = one_class_svm(S_3, out_fr);

S_lof = apply_lof(S_3, 10);
%% Top-K Analysis
[k_list(:,ind_outer), precision_gloss_3(:,ind_outer), recall_gloss_3(:,ind_outer), fpr_gloss_3(:,ind_outer)] = analyze_top_K(mahal_dist(S_3), X, ind_removed);
[~, precision_gl_svm(:,ind_outer), recall_gl_svm(:,ind_outer),fpr_gl_svm(:,ind_outer)] = analyze_top_K(S_svm, X, ind_removed, true);
[~, precision_gl_lof(:,ind_outer), recall_gl_lof(:,ind_outer),fpr_gl_lof(:,ind_outer)] = analyze_top_K(S_lof, X, ind_removed);
%% Visualize Decomposition66

if ~nodisp
    a_S = apply_lof(S_3, 10);
    plot_sensor_new(permute(X,[1,3,4,2]), permute(Yn,[1,3,4,2]), permute(mahal_S,[1,3,4,2]), permute(S_3,[1,3,4,2]), 3,30)
    
    figure,
    plot(fpr_gloss_en,recall_gloss_en,'DisplayName','GLOSS')
    legend, grid
    title('ROC of the envelope')
    
    figure,
    plot(k_list, precision_gloss_3,'DisplayName','GLOSS');
    legend, grid
    title('Precision')
    
    figure,
    plot(k_list, recall_gloss_3,'DisplayName','GLOSS');
    legend, grid
    title('Recall')
end