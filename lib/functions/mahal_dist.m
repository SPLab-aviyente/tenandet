function mahal_Y = mahal_dist(Y)

s = size(Y);
mahal_Y = zeros(s);
stand_Y = std(Y,[],3);
for i=1:s(1)
    for j=1:s(2)
        for k=1:s(4)
            mahal_Y(i,j,:,k) = Y(i,j,:,k)/((stand_Y(i,j,k)+0.1));
        end
    end
end
end