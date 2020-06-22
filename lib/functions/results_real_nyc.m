function results_real_nyc(S_gloss, S_wlrt, S_lrt, S_whorpca, S_horpca, Y, a_L, ind_removed, param, w_mahal)

list_k = [20,100,500,1000,2000,5000,7000,14000,21000];
det_gloss = zeros(size(list_k));
det_glos = det_gloss;
det_lrt = det_gloss;
det_whorpca = det_gloss;
det_horpca = det_gloss;
det_ee = det_gloss;
det_lof = det_gloss;
if w_mahal
    for i=1:length(list_k)
        det_gloss(i) = sum(detect_real_events(mahal_dist(S_wlrt), ind_removed, list_k(i)));
        det_glos(i) = sum(detect_real_events(mahal_dist(S_gloss), ind_removed, list_k(i)));
        det_lrt(i) = sum(detect_real_events(mahal_dist(S_lrt), ind_removed, list_k(i)));
        det_whorpca(i) = sum(detect_real_events(mahal_dist(S_whorpca), ind_removed, list_k(i)));
        det_horpca(i) = sum(detect_real_events(mahal_dist(S_horpca), ind_removed, list_k(i)));
        det_ee(i) = sum(detect_real_events(mahal_dist(Y-repmat(mean(Y,3),1,1,size(Y,3),1)), ind_removed, list_k(i)));
        det_lof(i) = sum(detect_real_events(a_L, ind_removed, list_k(i)));
    end
else
    for i=1:length(list_k)
        det_gloss(i) = sum(detect_real_events((S_wlrt), ind_removed, list_k(i)));
        det_glos(i) = sum(detect_real_events((S_gloss), ind_removed, list_k(i)));
        det_lrt(i) = sum(detect_real_events((S_lrt), ind_removed, list_k(i)));
        det_whorpca(i) = sum(detect_real_events((S_whorpca), ind_removed, list_k(i)));
        det_horpca(i) = sum(detect_real_events((S_horpca), ind_removed, list_k(i)));
        det_ee(i) = sum(detect_real_events(mahal_dist(Y-repmat(mean(Y,3),1,1,size(Y,3),1)), ind_removed, list_k(i)));
        det_lof(i) = sum(detect_real_events(a_L, ind_removed, list_k(i)));
    end
end

figure,
plot(100*list_k/numel(Y), det_gloss,'DisplayName','WLOSS','LineWidth',3)
hold
plot(100*list_k/numel(Y), det_glos,'DisplayName','GLOSS','LineWidth',3)
plot(100*list_k/numel(Y), det_lrt,'DisplayName','LOSS','LineWidth',3)
plot(100*list_k/numel(Y), det_whorpca,'DisplayName','WHORPCA','LineWidth',3)
plot(100*list_k/numel(Y), det_horpca,'DisplayName','HORPCA','LineWidth',3)
plot(100*list_k/numel(Y), det_ee,'DisplayName','EE','LineWidth',3)
plot(100*list_k/numel(Y), det_lof,'DisplayName','LOF','LineWidth',3)
legend
title(['\lambda ',num2str(param.lambda), ', \gamma ',num2str(param.gamma)])
grid
ax = gcf;
ax.CurrentAxes.FontSize=19;
ax.CurrentAxes.FontWeight = 'bold';
% plot_nyc(Y, S, S_gh, S_lrt, S_horpca, a_L, 1, 75)
end