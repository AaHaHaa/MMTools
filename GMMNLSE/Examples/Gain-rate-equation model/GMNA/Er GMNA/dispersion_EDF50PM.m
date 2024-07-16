close all; clearvars;

data = readmatrix('EDF50PM.xlsx','NumHeaderLines',1);
wavelength = data(:,1)*1e-9; % m
D = data(:,2)*1e-6; % dispersion in ps/(nm*km) in data to s/m^2 in D

c = 299792458; % m/s
beta2 = -wavelength.^2./(2*pi*c).*D; % s^2/m

figure;
plot(wavelength*1e9,beta2*1e27,'linewidth',2);
xlabel('Wavelength (nm)');
ylabel('\beta_2 (fs^2/mm)');
set(gca,'fontsize',20);
xlim([min(wavelength),max(wavelength)]*1e9);

omega = 2*pi*c./wavelength; % 1/s, rad
beta1 = flipud(cumtrapz(flipud(omega),flipud(beta2))) + 4875e-12; % s/m
beta0 = flipud(cumtrapz(flipud(omega),flipud(beta1))) + 5.2e6; % 1/m

figure;
plot(wavelength*1e9,beta0,'linewidth',2);
xlabel('Wavelength (nm)');
ylabel('\beta_0 (1/m)');
set(gca,'fontsize',20);
xlim([min(wavelength),max(wavelength)]*1e9);

save('beta_EDF50PM.mat','wavelength','beta0');