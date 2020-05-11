covid_data = load('covid_19data.mat', 'tensor');

Y = permute(covid_data.tensor, [2,3,4,1]);
X = zeros(size(Y));
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

param.lambda = 1/sqrt(max(size(Y)));
param.gamma = 1/3;
param.beta_1 = 1/(10*std(Y(:)));
param.beta_2 = 1/3;
param.beta_3 = 1/3;
param.max_iter = 100;
param.err_tol = 0.0001;

[L,S] = low_temp_sp_dec(Y, param);
plot_sensor_new(X, Y, S, L, 1)