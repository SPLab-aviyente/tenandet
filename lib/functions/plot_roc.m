function [auc, f, r] = plot_roc(fpr,  recall)
% auc = plot_roc(fpr, recall)
% Plots ROC curves of various methods of many experiments, including
% experiments on parameters like missing rate, anomaly amplitude, lambda
% and gamma.
% Inputs: 
%  fpr: False positive (FP) rates. FP/(FN+TN);
%  recall: Recalls. TP/(TP+FN)
n_exp = length(recall);
f = 0;
r = 0;
for i=1:n_exp
    f = f+permute(fpr{i},[2,3,4,5,1,6])./n_exp;
    r = r+permute(recall{i},[2,3,4,5,1,6])./n_exp;
end
len_gam = size(f,2);
len_lam = size(f,3);
len_n = size(f,4);
len_a = size(f,5);
auc = zeros(len_n, len_gam, len_lam, len_a, 10);

for l = 1:len_a
    for i = 1:len_n
        for j = 1:len_gam
            for k = 1:len_lam
%                 figure,
%                 plot(f(:,1,1,i,l,1), r(:,1,1,i,l,1),'Color','#0072BD','DisplayName','EE','LineWidth',7);
%                 hold on;
%                 plot(f(:,1,1,i,l,3), r(:,1,1,i,l,3),'Color','#EDB120','DisplayName','OCSVM','LineWidth',7)
%                 plot(f(:,1,1,i,l,2), r(:,1,1,i,l,2),'Color','#D95319','DisplayName','LOF','LineWidth',7);
%                 plot(f(:,1,k,i,l,4), r(:,1,k,i,l,4),'Color','#7E2F8E','DisplayName','HORPCA','LineWidth',7);
%                 plot(f(:,1,k,i,l,5), r(:,1,k,i,l,5),'--','Color','#77AC30','DisplayName','WHORPCA','LineWidth',7);
%                 plot(f(:,j,k,i,l,6), r(:,j,k,i,l,6),'--o','Color','#A2142F','DisplayName','CP','LineWidth',7,'MarkerSize',20);
%                 plot(f(:,j,k,i,l,6), r(:,j,k,i,l,7),'--+','Color','#4DBEEE','DisplayName','LOSS','LineWidth',7,'MarkerSize',20);
%                 plot(f(:,j,k,i,l,7), r(:,j,k,i,l,8),'--x','Color','#0072BD','DisplayName','GLOSS-EE','LineWidth',7,'MarkerSize',20)
%                 plot(f(:,j,k,i,l,8), r(:,j,k,i,l,9),'--x','Color','#EDB120','DisplayName','GLOSS-SVM','LineWidth',7,'MarkerSize',20)
%                 plot(f(:,j,k,i,l,9), r(:,j,k,i,l,10),'--x','Color','#D95319','DisplayName','GLOSS-LOF','LineWidth',7,'MarkerSize',20)
% %                 plot(f(:,j,k,i,l,10), r(:,j,k,i,l,10),'--x','Color','#A2142F','DisplayName','LOGSS','LineWidth',7,'MarkerSize',20)
%                 legend('location','southeast'), grid
% %                 title([num2str(round(100*24*n_missing(i)/prod(sizes))), '\lambda = ', num2str(lam_list_2(k)), '\gamma = ', num2str(gam_list_2(j)), ' 0.5 std'])
%                 ax = gcf;
%                 ax.CurrentAxes.FontSize = 26;
%                 ax.CurrentAxes.FontWeight = 'bold';
%                 for ii = 1:n_exp
                    auc(i,j,k,l,1) = trapz(f(:,1,1,i,l,1), r(:,1,1,i,l,1));
                    auc(i,j,k,l,2) = trapz(f(:,1,1,i,l,3), r(:,1,1,i,l,3));
                    auc(i,j,k,l,3) = trapz(f(:,1,1,i,l,2), r(:,1,1,i,l,2));
                    auc(i,j,k,l,4) = trapz(f(:,1,k,i,l,4), r(:,1,k,i,l,4));
                    auc(i,j,k,l,5) = trapz(f(:,1,k,i,l,5), r(:,1,k,i,l,5));
                    auc(i,j,k,l,6) = trapz(f(:,j,k,i,l,6), r(:,j,k,i,l,6));
                    auc(i,j,k,l,7) = trapz(f(:,j,k,i,l,7), r(:,j,k,i,l,7));
                    auc(i,j,k,l,8) = trapz(f(:,j,k,i,l,8), r(:,j,k,i,l,8));
                    auc(i,j,k,l,9) = trapz(f(:,j,k,i,l,9), r(:,j,k,i,l,9));
                    auc(i,j,k,l,10) = trapz(f(:,j,k,i,l,10), r(:,j,k,i,l,10));
%                 end
            end
        end
    end
end

end