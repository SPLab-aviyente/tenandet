
param.lambda = 1/sqrt(max(size(Y)));
param.beta_1 = 1/(5*std(Yn(setdiff(1:numel(Yn), ind_removed))));
param.max_iter = 100;
% param.err_tol = 0.001;
tic;
[L,S_horpca, obj_val] = horpca(Yn, param);
rmse_horpca(ind_outer) = norm(L-Y_gen)/sqrt(numel(Y));
mape_horpca(ind_outer) = sum(abs(L-Y_gen))/numel(Y);
time_horpca = toc
%% Envelope Analysis
if ind_outer == length(anom_list)
%     [fpr_horpca_en, recall_horpca_en] = analyze_envelope(S_horpca, X, ind_removed);
end
%% Top-K Analysis
[k_list(:,ind_outer), precision_horpca(:,ind_outer), recall_horpca(:,ind_outer),fpr_horpca(:,ind_outer)] = analyze_top_K(mahal_dist(S_horpca), X, ind_removed);
%% Visualize Decomposition
% plot_sensor_new(X, Yn, S, L, 5)
% plot_sensor_new(permute(X,[1,3,2,4]), permute(Yn,[1,3,2,4]), permute(S,[1,3,2,4]), permute(L,[1,3,2,4]), 5)
if ~nodisp
    plot_sensor_new(permute(X,[1,3,4,2]), permute(Yn,[1,3,4,2]), permute(S_horpca,[1,3,4,2]), permute(S_horpca,[1,3,4,2]), 3,10)
end

if ind_outer == length(anom_list)
    %     plot_sensor_new(permute(X,[1,3,2,4]), permute(Yn,[1,3,2,4]), permute(S,[1,3,2,4]), permute(L,[1,3,2,4]), 5)
end


if ~nodisp
    a_H = zeros(size(Y));
    for i=1:size(Yn,1)
        for j=1:size(Yn,2)
            for k=1:size(Yn,4)
                [~,a_H(i,j,:,k)] = LOF(squeeze(S_horpca(i,j,:,k))+randn(size(Yn,3))*0.001,5);
            end
        end
    end
    plot_sensor_new(permute(X,[1,3,4,2]), permute(Yn,[1,3,4,2]), permute(a_H,[1,3,4,2]), permute(S_horpca,[1,3,4,2]), 3,10)
    
    figure,
    plot(fpr_horpca_en,recall_horpca_en,'DisplayName','HORPCA')
    legend, grid
    title('ROC of the envelope')
    
    figure,
    plot(k_list, precision_horpca,'DisplayName','HORPCA');
    legend, grid
    title('Precision')
    
    figure,
    plot(k_list, recall_horpca,'DisplayName','HORPCA');
    legend, grid
    title('Recall')
end