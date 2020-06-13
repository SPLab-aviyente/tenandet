function [k_list, precision, recall, fpr] = analyze_top_K(S, X, ind_removed)
% [k_list, precision, recall] = analyze_top_K(S, X, ind_removed)
% Function that provides analysis of precision and recall of top estimates
% of anomaly based on amplitude.
ind_rem = setdiff(1:numel(S), ind_removed);
S(isnan(S)) = 0;
[~, ind] = sort(abs(S(ind_rem)),'descend');
ind = ind_rem(ind);
k_list = [1000];%,5,10,50,100:500:3000];%, 4000:1000:9000, 10^4:10^4:10^5];%, 10^5:10^5:7*10^5];
precision = zeros(length(k_list),1);
recall = precision;
fpr = recall;
for i=1:length(k_list)
    true_pos = sum(X(ind(1:k_list(i))));
    false_pos = k_list(i)-true_pos;
    false_neg = sum(X,'all')-true_pos;
    true_neg = sum(1-X,'all')-false_pos;
    precision(i) = true_pos/(true_pos+false_pos);
    recall(i) = true_pos/(true_pos+false_neg);
    fpr(i) = false_pos/(false_pos+true_neg);
end
end