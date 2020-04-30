
tensor = read_covid_csv('covid_confirmed_usafacts.csv', 'MI', 14);
tensor = tensor(:,8:14,:);
[deceased_tensor, dates, county_names] = read_covid_csv('covid_deaths_usafacts.csv', 'MI', 14);
deceased_tensor =deceased_tensor(:,8:end,:);
tensor = permute(cat(4, tensor, deceased_tensor),[4,1,2,3]);
save('covid_19data', 'dates', 'county_names', 'tensor');