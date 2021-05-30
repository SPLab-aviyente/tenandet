
load('preliminary_results\logss__miss.mat')
n_exp = 10;
len_lam = length(lam_list_1);
len_gam = length(gam_list_1);
len_n = length(n_missing);
len_a = length(amp_anom);
auc = zeros(len_n,len_lam,len_gam,n_exp,4);
times = auc;
for i=1:n_exp
    for i_L = 1:len_lam
        for i_G = 1:len_gam
            for i_N = 1:len_n
                auc(i_N, i_L, i_G, i, 1) = trapz(fpr{i_N}(i,:,1,1,1,1), recall{i_N}(i,:,1,1,1,1));
                times(i_N, i_L, i_G, i, 1) = time_ee(i,1,i_N);
                auc(i_N, i_L, i_G, i, 2) = trapz(fpr{i_N}(i,:,1,i_L,1,4), recall{i_N}(i,:,1,i_L,1,4));
                times(i_N, i_L, i_G, i, 2) = time_horpca{i_N}(i, 1, i_L);
                auc(i_N, i_L, i_G, i, 3) = trapz(fpr{i_N}(i,:,i_G,i_L,1,6), recall{i_N}(i,:,i_G,i_L,1,6));
                times(i_N, i_L, i_G, i, 3) = times_loss{i_N}(i, 1, i_L, i_G);
                auc(i_N, i_L, i_G, i, 4) = trapz(fpr{i_N}(i,:,i_G,i_L,1,8), recall{i_N}(i,:,i_G,i_L,1,8));
                times(i_N, i_L, i_G, i, 4) = times_logss{i_N}(i, 1, i_L, i_G);
            end
        end
    end
end