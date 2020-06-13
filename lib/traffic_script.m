
[X, Y, Yn, days, weeks, sensors] = get_traffic_data;

gamma_range = 10.^(-2:0.2:0);
param.lambda = 1/sqrt(max(size(Y)));
param.theta = param.lambda/100;
param.alpha = param.lambda;
param.beta_1 = 1/(5*std(Yn(:)));
param.beta_2 = param.beta_1;
param.beta_3 = param.beta_1;
param.beta_4 = param.beta_1;
param.max_iter = 100;
param.err_tol = 0.01;

for i=1:length(gamma_range)
    param.gamma = 1/sqrt(max(size(Y)));
    t = cputime;
    [L,S,N,obj_val] = low_temp_sp_dec(Yn, param);
    t = cputime-t
    plot_sensor_new(X, Y, S, L, 5)
    t = cputime;
    [L,S, N,obj_val_gloss] = gloss(Yn, param);
    t = cputime-t
    plot_sensor_new(X, Y, S, L, 5)
    recall(i) = sum(X>0.5&(S~=0), 'all')/sum((X>0.5),'all');
    fpr(i) = sum(X<0.5&(S~=0), 'all')/sum((X<0.5),'all');
end
figure,
plot(fpr, recall,'LineWidth', 3)
grid
% thresh = 0.7;
% myDec = myTensorDecTree(Y,thresh,'LSA',[1,1,1,1]); % enter parameters
% myDec.prepareRoot; % prepare root node
% myDec.addLayers(1); % add scales
% myDec.cReconstruction(2);% 1 2 3 4
% S = Y-myDec.recData;
% 
% stds = 0.01:0.1:4;
% recall = zeros(1,length(stds));
% fpr = recall;
% for i = 1:length(stds)
%     num_stds = stds(i);
%     mean_S = mean(S,3);
%     std_S = std(S,[],3);
%     lbl_S = S<mean_S-num_stds*std_S | S>mean_S+num_stds*std_S;
%     recall(i) = sum((X>0.5)&lbl_S,'all')/sum((X>0.5),'all');
%     fpr(i) = sum((X<0.5)&lbl_S,'all')/sum((X<0.5),'all');
% end
% figure,
% plot(fpr, recall,'LineWidth', 3)
% grid
