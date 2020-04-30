function [tensor, dates, county_names] = read_covid_csv(filename, varargin)
% [tensor, dates, county_names] = read_covid_csv(filename, varargin)
% reads covid19 data and produces tensor with size: 
%      number_counties x 7(days in a week) x number_weeks
if isempty(varargin{2})
    numweeks = 14;
else
    numweeks = varargin{2};
end
if isempty(varargin{1})
    state = 'MI';
else
    state = varargin{1};
end
if nargin ==0
    filename = 'covid_confirmed_usafacts.csv';
end
covid19 = readtable(filename);
mask = cellfun(@contains, table2cell(covid19(:,3)), repmat({state}, size(covid19,1),1));

county_names = table2cell(covid19(mask,2));
dates = table2cell(covid19(1,5:end));

mat = (table2cell(covid19(mask, 5:end)));
mat = cellfun(@str2num, mat);
% mat = mat(:,49:end); % Clear zeros for many states. 
% dates = dates(49:end);
tensor = permute(reshape(mat(:,1:numweeks*7),[],7,numweeks),[2,3,1]);
% save('covid_19data', 'dates', 'county_names', 'tensor');
end