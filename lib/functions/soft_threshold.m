function X = soft_threshold(X, sigma)
%softThresh Soft thresholding function
%   X = soft_threshold(X, sigma)
X(abs(X)<sigma) = 0;
temp = abs(X)-sigma;
temp(temp<0) = 0;
X = temp.*sign(X);
% X(abs(X)>sigma) = X(abs(X)>sigma) - sign(X(abs(X)>sigma)).*sigma;
end