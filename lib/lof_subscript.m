function [a_L, precision, recall, fpr] = lof_subscript(Y, X, param)
%
%
a_L = apply_lof(Y, 10);


% Top-K Analysis
[~, precision, recall,fpr] = analyze_top_K(a_L, X, param.ind_m);
end