function [S, precision, recall, fpr, rmse, mape]= ocsvm_subscript(Y, X, param)

out_fr = 0.1;
S = one_class_svm(Y, out_fr);

% Top-K Analysis
[~, precision, recall, fpr] = analyze_top_K(S, X, param.ind_m, true);
end