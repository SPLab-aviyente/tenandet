% param.lambda = 1/4/sqrt(max(size(Yn)));
% param.lambda = 1/2/sqrt(max(size(Y)));
% param.gamma = 1/15/sqrt(max(size(Y)));
param.alpha = 0;
param.beta_1 = 1/(5*std(Yn(setdiff(1:numel(Yn), ind_removed))));
param.beta_2 = param.beta_1;
param.beta_3 = param.beta_1;
param.max_iter = 100;
% param.err_tol = 0.001;

t = tic;
[L, S_lrt, N] = low_temp_sp_dec(Yn, param);
rmse_loss(ind_outer) = norm(L(:)-Y_gen(:))/sqrt(numel(Y_gen));
mape_loss(ind_outer) = sum(abs(L-Y_gen),'all')/numel(Y_gen);
time_lrtssd = toc(t)
%% Envelope analysis
if ind_outer == length(anom_list)
%     [fpr_lrs, recall_lrs_en] = analyze_envelope(S_lrt, X, ind_removed); 
end
%% Top-K Analysis
mahal_S = mahal_dist(S_lrt);
[~, precision_lrs(:,ind_outer), recall_lrs(:,ind_outer), fpr_lrs(:,ind_outer)] = analyze_top_K(mahal_S, X, ind_removed);
%% Visualized Decomposition
if ind_outer == length(anom_list)
%     plot_sensor_new(permute(X,[1,3,4,2]), permute(Yn,[1,3,4,2]), permute(mahal_S,[1,3,4,2]), permute(S_lrt,[1,3,4,2]), 3, 10)
end
