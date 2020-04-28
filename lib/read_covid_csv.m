
covid19 = readtable('covid_confirmed_usafacts.csv');
mask = cellfun(@contains, table2cell(covid19(:,3)), repmat({'MI'}, size(covid19,1),1));

county_names=table2cell(covid19(mask,2));
dates = table2cell(covid19(1,5:end));

mat = (table2cell(covid19(mask, 5:end)));
mat = cellfun(@str2num, mat);
mat = mat(:,49:end);
dates = dates(49:end);
tensor = permute(reshape(mat,[],7,7),[2,3,1]);
save('covid_19data', 'dates', 'county_names', 'tensor');