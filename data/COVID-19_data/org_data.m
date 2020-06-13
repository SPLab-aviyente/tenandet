
tensor = read_covid_csv('covid_confirmed_usafacts.csv', 'MI', 17);
tensor = reshape(tensor(:,8:end,:), [], 83);
[deceased_tensor, dates, county_names] = read_covid_csv('covid_deaths_usafacts.csv', 'MI', 17);
deceased_tensor = reshape(deceased_tensor(:,8:end,:), [], 83);
tensor = permute(cat(3, tensor, deceased_tensor), [1,3,2]);
save('covid_19data', 'dates', 'county_names', 'tensor');