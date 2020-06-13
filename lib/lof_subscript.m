

tic;
a_L = apply_lof(Yn, 10);
time_LOF = toc

%% Envelope Analysis
if ind_outer == length(anom_list)
%     [fpr_lof, recall_lof_en] = analyze_envelope(a_L, X, ind_removed);
end

%% Top-K Analysis
[~, precision_lof(:,ind_outer), recall_lof(:,ind_outer),fpr_lof(:,ind_outer)] = analyze_top_K(a_L, X, ind_removed);
%% Visualized Decomposition
if ind_outer == length(anom_list)
%     plot_sensor_new(permute(X,[1,3,4,2]), permute(Yn,[1,3,4,2]), permute(a_L,[1,3,4,2]), permute(L,[1,3,4,2]), 3,30)
end

if ~nodisp    
    figure,
    plot(fpr_gloss_en,recall_lof_en,'DisplayName','LOF')
    legend, grid
    title('ROC of the envelope')
    
    figure,
    plot(k_list, precision_lof,'DisplayName','LOF');
    legend, grid
    title('Precision')
    
    figure,
    plot(k_list, recall_lof,'DisplayName','LOF');
    legend, grid
    title('Recall')
end