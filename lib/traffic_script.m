traffic_data = load('24-Oct-2018_data.mat');

X = create_labels(traffic_data);
X = reshape(X, 288, 365, 10);
X = reshape(X(:,1:364,:), 288, 7, 52, 10);
Y = reshape(traffic_data.station_counts, 288, 365, 10);
Y = reshape(Y(:,1:364,:), 288, 7, 52, 10);
Y = denan(Y);

% save('traffic_data', 'Y', 'X')
% subplot(1,2,1)
% hold on
% for i=1:7
% plot(X(:,i,1,1))
% end
% subplot(1,2,2)
% hold on
% for i=1:7
% plot(Y(:,i,1,1))
% end
gamma_range = 10.^(-2:0.2:0);
param.lambda = 1/sqrt(max(size(Y)));

param.beta_1 = 1/(10*std(Y(:)));
param.beta_2 = param.beta_1;
param.beta_3 = param.beta_1;
param.max_iter = 100;
param.err_tol = 0.0001;

for i=1:length(gamma_range)
    param.gamma = param.lambda;
    [L,S] = low_temp_sp_dec(Y, param);
%     plot_sensor(X, Y, S, L, 5)
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

function [X] = create_labels(data)
% [X] = createLabels(data)
%   Creates labels using the information provided in the data.
num_elements = size(data.station_counts,1);
num_sensors  = size(data.station_counts,2);
X = zeros(num_elements, num_sensors);
for i=1:length(data.station_ids)
    ind_start = round(minutes(data.inc_timestamp(i,:)-...
        datetime(data.station_times(1,:)))/5);
    ind_end = round(data.inc_duration(i)/5)+ind_start;
    [~,ind_station,~] = intersect(data.station_ids,data.inc_station_id(i));
    X(ind_start:ind_end, ind_station) = 1;
end
end