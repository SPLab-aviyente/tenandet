
load('M:\Documents\MATLAB\tensor_anomaly_detection\data\guanzhou_traffic\Data\tensor.mat')
Y = reshape(permute(tensor(:,1:56,:),[3,2,1]), 144, 7, 8, []);
X = zeros(size(Y));
Yn = Y;
Yn(:,:,5,5) = 0;

gamma_range = 10.^(-2:0.2:0);
param.lambda = 1/2/sqrt(max(size(Y)));
param.gamma = 1/2/sqrt(max(size(Y)));
param.alpha = param.lambda/1000;
param.beta_1 = 1/(5*std(Y(:)));
param.beta_2 = param.beta_1;
param.beta_3 = param.beta_1;
param.max_iter = 100;
param.err_tol = 0.0001;

[L,S,N] = low_temp_sp_dec(Yn, param);
plot_sensor_new(Y(:,:,:,:), Yn(:,:,:,:), L(:,:,:,:), S(:,:,:,:), 5)