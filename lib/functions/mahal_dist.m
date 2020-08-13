function [mahal_Y, est_Y] = mahal_dist(Y)

s = size(Y);
mahal_Y = zeros(s);
est_Y = mahal_Y;
% stand_Y = std(Y,[],3);
% est_Y = repmat(mean(Y,3),1,1,s(3),1);
% mahal_Y = (Y-est_Y)./repmat(stand_Y+0.1,1,1,s(3),1);

for i=1:s(1)
    for j=1:s(2)
        for k=1:s(4)
            [~,est_Y(i,j,:,k),mahal_Y(i,j,:,k)] = robustcov(squeeze(Y(i,j,:,k)));
            if sum(mahal_Y(i,j,:,k)~=0,'all')==0
                [~,est_Y(i,j,:,k),mahal_Y(i,j,:,k)] = robustcov(squeeze(Y(i,j,:,k))+sqrt(0.1)*randn(s(3),1));
            end
        end
    end
end
end