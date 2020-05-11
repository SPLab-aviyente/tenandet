function plot_sensor_new(X, Y, S, S_lbl, checkSens)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for i=1:size(Y,3)
    figure,
    subplot(4,1,1)
    plot(X(:,:,i,checkSens))
    subplot(4,1,2)
    plot(Y(:,:,i,checkSens))
    subplot(4,1,3)
    plot((S(:,:,i,checkSens)))
    subplot(4,1,4)
    plot((S_lbl(:,:,i,checkSens)))
end

end