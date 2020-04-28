function [X] = createLabel(data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
nElements = size(data.station_counts,1);
nSensors  = size(data.station_counts,2);
X = zeros(nElements, nSensors);
for i=1:length(data.inc_station_id)
    indE1         = round(minutes(data.inc_timestamp(i,:)-...
        datetime(data.station_times(1,:)))/5);
    indE2         = round(data.inc_duration(i)/5)+indE1;
    [~,indS,]     = intersect(data.station_ids,data.inc_station_id(i));
    X(indE1:indE2, indS) = 1;
end

end

