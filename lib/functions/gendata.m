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
    if length(varargin) == 2
        num_missing_days = varargin{2};
    else
        num_missing_days = 5000;
    end
    if length(varargin) == 3
        pt_vs_fb = varargin{3};
    else
        pt_vs_fb = false;
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
    pt_vs_fb = false;
end

mat_anomaly = zeros(number_anomaly,5);
[r,c,z] = ind2sub(dims(2:4),randperm(prod(dims(2:4)),number_anomaly));
mat_anomaly(:,3:5) = [r',c',z'];
for i=1:number_anomaly
    mat_anomaly(i,1:2) = randi(dims(1)-length_anomaly)+[0, length_anomaly];
end
X = zeros(size(Y));
Y = Y.*(0.3*randn(size(Y))+1);
Yn = Y;
ind_anoms = zeros(number_anomaly*length_anomaly,1);
for i=1:number_anomaly
    indCell = num2cell([(mat_anomaly(i,1):mat_anomaly(i,2))',repmat(mat_anomaly(i,3:5),mat_anomaly(i,2)-mat_anomaly(i,1)+1,1)],1);
    indVec = sub2ind(dims, indCell{:});
    if pt_vs_fb
        ind_anoms((i-1)*length_anomaly+1:i*length_anomaly) = indVec;
    end
    add_sub = sign(randn);
    Yn(indVec)=(amplitude_anomaly^add_sub)*Yn(indVec)+add_sub*amplitude_anomaly*log(mean(Yn(indVec)))^2;
    X(indVec)=1;
end

if pt_vs_fb
    removed_inds = setdiff(1:numel(Yn), ind_anoms);
    removed_inds = removed_inds(randperm(length(removed_inds),num_missing_days*size(Yn,1)));
    Yn(removed_inds)=0;
else
    [days,weeks,sens] = ind2sub(dims(2:4),randperm(prod(dims(2:4)),num_missing_days));
    removed_inds = [];
    for i=1:num_missing_days
        if ~isempty(intersect([days(i),weeks(i),sens(i)], mat_anomaly(:,3:5),'rows'))
            temp = setdiff(1:dims(4),sens(i));
            sens(i) = temp(randi(dims(4)-1));
        end
        indCell = num2cell([(1:dims(1))',repmat([days(i),weeks(i),sens(i)],dims(1),1)],1);
        indVec = sub2ind(dims,indCell{:});
        Yn(indVec) = 0;
        removed_inds = [removed_inds; indVec];
    end
end
end