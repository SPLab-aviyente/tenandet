
out_fr = 0.1;
timer_ocsvm = tic;
S_svm = one_class_svm(Yn, out_fr);
time_ocsvm = toc(timer_ocsvm)

%% Top-K Analysis
[~, precision_svm(:,ind_outer), recall_svm(:,ind_outer),fpr_svm(:,ind_outer)] = analyze_top_K(S_svm, X, ind_removed, true);