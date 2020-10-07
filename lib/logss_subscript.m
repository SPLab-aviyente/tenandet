function [L, S, precision, recall, fpr, time_l, times, rmse, mape] =  logss_subscript(Y, Y_gen, X, param)

tic;
[L, S, times] = logss(Y, param);
time_l = toc;
rmse = norm(L(:)-Y_gen(:))/sqrt(numel(Y_gen));
mape = sum(abs(L-Y_gen),'all')/numel(Y_gen);

%% Top-K Analysis
[~, precision, recall, fpr] = analyze_top_K(mahal_dist(S), X, param.ind_m);
end