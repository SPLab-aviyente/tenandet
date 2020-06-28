
% param.lambda = 1/2/sqrt(max(size(Yn)));
% param.gamma = 1/(max(size(Y)));
param.theta = 1/(sum(std(t2m(Yn,4),[],2).^2));
param.alpha = 0;
% param.lambda = 1/5/sqrt(max(size(Y)));
% param.gamma = 3/sqrt(max(size(Y)));
% param.theta = param.lambda/5000;
% param.lambda = 15/sqrt(max(size(Y)));
% param.gamma = 15/sqrt(max(size(Y)));
% param.theta = param.lambda/1000;
% param.alpha = param.lambda/(std(Y(:))^2);
% param.psi = [.01,1,10,.001];
param.beta_1 = 1/(5*std(Yn(setdiff(1:numel(Yn), ind_removed))));
param.beta_2 = param.beta_1/1;
param.beta_3 = param.beta_1/1;
param.beta_4 = param.beta_1/1;
param.max_iter = 100;
tic;
[L,S,~, obj_val] = gloss(Yn, param);
rmse_gloss(ind_outer) = norm(L(:)-Y_gen(:))/sqrt(numel(Y_gen));
mape_gloss(ind_outer) = sum(abs(L-Y_gen),'all')/numel(Y_gen);
time_gloss = toc
%% Envelope Analysis
if ind_outer == length(anom_list)
%     [fpr_gloss_en, recall_gloss_en] = analyze_envelope(S, X, ind_removed);
end
%% Top-K Analysis
[k_list(:,ind_outer), precision_gloss(:,ind_outer), recall_gloss(:,ind_outer), fpr_gloss(:,ind_outer)] = analyze_top_K(mahal_dist(S), X, ind_removed);
%% Visualize Decomposition
% plot_sensor_new(X, Yn, S, L, 5)
% plot_sensor_new(permute(X,[1,3,2,4]), permute(Yn,[1,3,2,4]), permute(S,[1,3,2,4]), permute(L,[1,3,2,4]), 5)
if ~nodisp
%     plot_sensor_new(permute(X,[1,3,4,2]), permute(Yn,[1,3,4,2]), permute(S,[1,3,4,2]), permute(L,[1,3,4,2]), 3,30)
end
if ind_outer == length(anom_list)
    %     plot_sensor_new(permute(X,[1,3,2,4]), permute(Yn,[1,3,2,4]), permute(S,[1,3,2,4]), permute(L,[1,3,2,4]), 5)
end
% a_S = apply_lof(mahal_dist(S), 10);
% plot_sensor_new(permute(a_S,[1,3,4,2]), permute(Yn,[1,3,4,2]), permute(mahal_dist(S),[1,3,4,2]), permute(S,[1,3,4,2]), 3,30)
if ~nodisp
    a_S = apply_lof(S, 10);
    plot_sensor_new(permute(X,[1,3,4,2]), permute(Yn,[1,3,4,2]), permute(mahal_S,[1,3,4,2]), permute(S,[1,3,4,2]), 3,30)
    
    figure,
    plot(fpr_gloss_en,recall_gloss_en,'DisplayName','GLOSS')
    legend, grid
    title('ROC of the envelope')
    
    figure,
    plot(k_list, precision_gloss,'DisplayName','GLOSS');
    legend, grid
    title('Precision')
    
    figure,
    plot(k_list, recall_gloss,'DisplayName','GLOSS');
    legend, grid
    title('Recall')
end