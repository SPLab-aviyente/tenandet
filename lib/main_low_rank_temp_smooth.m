
sizes   = [144,7,52,10];
numAnom = 100;
lenAnom = 25;
amoAnom = 3;
modes   = [1,3];
[X, Y, Yn]   = gendata(sizes, numAnom, lenAnom, amoAnom);

param.lambda = 1/3/sqrt(max(size(Y)));
param.gamma = 1/5/sqrt(max(size(Y)));
param.alpha = param.lambda;
param.beta_1 = 1/(5*std(Y(:)));
param.beta_2 = param.beta_1;
param.beta_3 = param.beta_1;
param.max_iter = 100;
param.err_tol = 0.001;

[L,S, N] = low_temp_sp_dec(Yn, param);
plot_sensor_new(Yn, Y, S, L, 5)