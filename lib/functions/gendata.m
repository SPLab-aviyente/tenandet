function [X, Y] = gendata(dims, numAn, lengthAn, ampAn)
%gendata Generate four-mode simulation tensor for anomaly detection.
%   [X, Y] = gendata(dims, numAn, lengthAn, ampAn)
%
%   Inputs:
%   dims     : Dimensionality of the tensor in vector form.
%   numAn    : Number of anomalies in the data.
%   lengthAn : Length of anomalies.
%   ampAn    : Amplitude of anomalies.
%
%   Outputs:
%   X        : Labels of anomalies. 1 if there is an anomaly, -1 otherwise.
%   Y        : Data with anomalies.
%
%   example:
%   sizes   = [50,30,15,10];
%   numAnom = 1000;
%   lenAnom = 25;
%   amoAnom = 3;
%   modes   = 1:3;
%   [X,Y]   = gendata(sizes, numAnom, lenAnom, amoAnom);

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
indMatrix = zeros(numAn,5);
for i=1:numAn
    indMatrix(i,:)=[randi(dims(1)-lengthAn)+[0, lengthAn], randi(dims(2)),randi(dims(3)),randi(dims(4))];
end
X = -ones(size(Y));
[~,tempInd,~]=unique(indMatrix(:,3:5),'rows');
indMatrix = indMatrix(tempInd,:);
for i=1:size(indMatrix, 1)
    indCell=num2cell([(indMatrix(i,1):indMatrix(i,2))',repmat(indMatrix(i,3:5),indMatrix(i,2)-indMatrix(i,1)+1,1)],1);
    indVec=sub2ind(size(Y),indCell{:});
    Y(indVec)=ampAn;%Y(indVec)-2;
    X(indVec)=1;
end
end