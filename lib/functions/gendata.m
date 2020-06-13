function [X, Y, Yn, removed_inds, mat_anomaly] = gendata(dims, number_anomaly, length_anomaly, amplitude_anomaly, varargin)
%gendata Generate four-mode simulation tensor for anomaly detection.
%   [X, Y, Yn] = gendata(dims, numAn, lengthAn, ampAn)
%
%   Inputs:
%   dims : Dimensionality of the tensor in vector form.
%   number_anomaly : Number of anomalies in the data.
%   length_anomaly : Length of anomalies.
%   amplitude_anomaly : Amplitude of anomalies.
%
%   Outputs:
%   X        : Labels of anomalies. 1 if there is an anomaly, -1 otherwise.
%   Y        : Data with anomalies.
%   Yn       : Noisy anomaly data.
%
%   example:
%   sizes   = [50,30,15,10];
%   number_anomaly = 1000;
%   length_anomaly = 25;
%   amplitude_anomaly = 3;
%   modes   = 1:3;
%   [X,Y]   = gendata(sizes, number_anomaly, length_anomaly, amplitude_anomaly);
rng(123);
if ~isempty(varargin)
    data = varargin{1};
    dims = size(data);
    Y = repmat(mean(data,3),[1,1,dims(3),1]);
    if length(varargin) ==2
        num_missing_days = varargin{2};
    else
        num_missing_days = 5000;
    end
else
    var=0.01;
    t=(1:dims(1))/dims(1);
    Y = zeros(dims);
    for k=1:dims(4)
        for j=1:dims(3)
            for i=1:dims(2)
                Y(:,i,j,k) = sin(pi*t*k/dims(4))+cos(pi*2*t*j/dims(3))+cos(2*pi*t+randn(1)*sqrt(var));
            end
        end
    end
end

mat_anomaly = zeros(number_anomaly,5);
for i=1:number_anomaly
    mat_anomaly(i,:)=[randi(dims(1)-length_anomaly)+[0, length_anomaly], randi(dims(2)),randi(dims(3)),randi(dims(4))];
end
X = zeros(size(Y));
Y = Y.*(0.3*randn(size(Y))+1);
[~,tempInd,~]=unique(mat_anomaly(:,3:5),'rows');
mat_anomaly = mat_anomaly(tempInd,:);
for i=1:size(mat_anomaly, 1)
    indCell=num2cell([(mat_anomaly(i,1):mat_anomaly(i,2))',repmat(mat_anomaly(i,3:5),mat_anomaly(i,2)-mat_anomaly(i,1)+1,1)],1);
    indVec=sub2ind(dims, indCell{:});
    Y(indVec)=amplitude_anomaly*Y(indVec)+amplitude_anomaly*log(mean(Y(indVec)))^2;
    X(indVec)=1;
end


Yn = Y;
% Y = Y-min(Yn(:));
% Yn = Yn-min(Yn(:));
days = randi(dims(3)*dims(2), num_missing_days, 1)-0.1;
weeks = floor(days/7)+1;
days = floor(mod(days, 7))+1;
removed_inds = [];
for i=1:num_missing_days
    rand_sens = randi(dims(4));
    if ~isempty(intersect([days(i),weeks(i),rand_sens], mat_anomaly(:,3:5),'rows'))
        temp = setdiff(1:dims(4),rand_sens);
        rand_sens = temp(randi(dims(4)-1));
    end
    indCell = num2cell([(1:dims(1))',repmat([days(i),weeks(i),rand_sens],dims(1),1)],1);
    indVec = sub2ind(dims,indCell{:});
    Yn(indVec) = 0;
    removed_inds = [removed_inds; indVec];
end
end