function [L, S, precision, recall, fpr, rmse, mape] =  logss_subscript(Y, Y_gen, X, param)

param.max_iter = 100;
[L,S,~] = logss(Y, param);
rmse = norm(L(:)-Y_gen(:))/sqrt(numel(Y_gen));
mape = sum(abs(L-Y_gen),'all')/numel(Y_gen);

%% Top-K Analysis
[k_list, precision, recall, fpr] = analyze_top_K(mahal_dist(S), X, param.ind_m);
end