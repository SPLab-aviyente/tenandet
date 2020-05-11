
sizes   = [50,30,15,10];
numAnom = 100;
lenAnom = 25;
amoAnom = 5;
modes   = [1,3];
[X,Y]   = gendata(sizes, numAnom, lenAnom, amoAnom);
% traffic_script
% Y = permute(Y, [1,3,2,4]);
% X = permute(X, [1,3,2,4]);

thresh = 0.7;
myDec = myTensorDecTree(Y,thresh,'LSA',[1,1,1,1]); % enter parameters
myDec.prepareRoot; % prepare root node
myDec.addLayers(1); % add scales
myDec.cReconstruction(2);% 1 2 3 4
S = Y-myDec.recData;
% 
num_stds = 2;
mean_S = mean(S,3);
std_S = std(S,[],3);
lbl_S = S<mean_S-num_stds*std_S | S>mean_S+num_stds*std_S;
% 
% stds = 0.01:0.1:4;
% recall = zeros(1,length(stds));
% fpr = recall;
% for i = 1:length(stds)
%     num_stds = stds(i);
%     mean_S = mean(S,3);
%     std_S = std(S,[],3);
%     lbl_S = S<mean_S-num_stds*std_S | S>mean_S+num_stds*std_S;
%     recall(i) = sum((X>0)&lbl_S,'all')/sum((X>0),'all');
%     fpr(i) = sum((X<0.5)&lbl_S,'all')/sum((X<0.5),'all');
% end
% figure,
% plot(fpr,recall,'LineWidth', 3)
% grid
plot_sensor_new(X, Y, S, lbl_S, 5)
gamma_range = 10.^(-2:0.2:0);
param.lambda = 1/sqrt(max(size(Y)));

param.beta_1 = 1/(10*std(Y(:)));
param.beta_2 = 1/3;
param.beta_3 = 1/3;
param.max_iter = 100;
param.err_tol = 0.0001;

for i=1:length(gamma_range)
    param.gamma = gamma_range(i);
    [L,S] = low_temp_sp_dec(Y, param);
%     plot_sensor(X, Y, S, L, 5)
    recall(i) = sum(X>0.5&(S~=0), 'all')/sum((X>0.5),'all');
    fpr(i) = sum(X<0.5&(S~=0), 'all')/sum((X<0.5),'all');
end
figure,
plot(fpr, recall,'LineWidth', 3)
grid
% [B, U, objVal]  = learnMappings(Y(:,:,:,1:5), X(:,:,:,1:5), modes);
% 
% checkSens = 5;
% plot_sensor(X(:,:,:,6:end), Y(:,:,:,6:end), B, U, checkSens, modes)