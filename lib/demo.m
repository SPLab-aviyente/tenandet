
sizes   = [50,30,15,10];
numAnom = 100;
lenAnom = 25;
amoAnom = 5;
modes   = [1,3];
[X,Y]   = gendata(sizes, numAnom, lenAnom, amoAnom);
% traffic_script
% Y = permute(Y, [1,3,2,4]);
% X = permute(X, [1,3,2,4]);

[B, U, objVal]  = learnMappings(Y(:,:,:,1:5), X(:,:,:,1:5), modes);

checkSens = 5;
plot_sensor(X(:,:,:,6:end), Y(:,:,:,6:end), B, U, checkSens, modes)