function plot_covid(Y, S, L, counties)
% plot_covid(Y, S, L, counties)
%   Plot covid data with sparse and low rank parts.

for i=1:length(counties)
    figure,
    subplot(3,1,1)
    plot(Y(:,counties(i)))
    subplot(3,1,2)
    plot((S(:,counties(i))))
    subplot(3,1,3)
    plot((L(:,counties(i))))
end

end