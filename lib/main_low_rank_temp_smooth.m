
sizes   = [50,30,15,10];
numAnom = 100;
lenAnom = 25;
amoAnom = 3;
modes   = [1,3];
[X,Y]   = gendata(sizes, numAnom, lenAnom, amoAnom);

param.lambda = 1/sqrt(max(size(Y)));
param.gamma = 1/3;
param.beta_1 = 1/(10*std(Y(:)));
param.beta_2 = 1/3;
param.beta_3 = 1/3;
param.max_iter = 100;
param.err_tol = 0.0001;

[L,S] = low_temp_sp_dec(Y, param);
plot_sensor_new(X, Y, S, S~=0, 5)