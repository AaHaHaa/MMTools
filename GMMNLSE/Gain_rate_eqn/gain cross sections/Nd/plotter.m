clearvars; close all;

data = readmatrix('Nd (Martin thesis 2006).xlsx');
wavelength = data(:,1);

GSA011 = data(:,2); GSA011(wavelength>445.5e-9) = 0; % 2D5/2 + 2P1/2
GSA010 = data(:,2); GSA010(wavelength<445.5e-9 | wavelength>492.4e-9) = 0; % 2G11/2 + 2(D,P)3/2 + 2G9/2 + 2K15/2
GSA09 = data(:,2); GSA09(wavelength<492.4e-9 | wavelength>555.4e-9) = 0; % 4G9/2 + 4G7/2 + 2K13/2
GSA08 = data(:,2); GSA08(wavelength<555.4e-9 | wavelength>643e-9) = 0; % 2G7/2 + 4G5/2
GSA07 = data(:,2); GSA07(wavelength<643e-9 | wavelength>723.5e-9) = 0; % 4F9/2
GSA06 = data(:,2); GSA06(wavelength<723.5e-9 | wavelength>778.4e-9) = 0; % 4S3/2 + 4F7/2
GSA05 = data(:,2); GSA05(wavelength<778.4e-9 | wavelength>848.6e-9) = 0; % 2H9/2 + 4F5/2
GSA04 = data(:,2); GSA04(wavelength<848.6e-9 | wavelength>1000e-9) = 0; % 4F3/2

data2 = readmatrix('1.csv');
data2 = interp1(data2(:,1),data2(:,2),wavelength*1e9,'linear',0)/1225.69*0.793125e-24;
data2(wavelength<1652.15e-9 | wavelength>1696.6e-9) = 0;
GSA03 = smooth(data2,20); % 4I15/2

data3 = readmatrix('2.csv');
data3 = interp1(data3(:,1),data3(:,2),wavelength*1e9,'linear',0)*0.669461e-24;
data3(wavelength<1265e-9 | wavelength>1401e-9) = 0;
ESA49 = data3;
for i = 1:50
    ESA49 = smooth(ESA49,10); % 4G9/2 + 4G7/2
    ESA49(wavelength<1280e-9 | wavelength>1410e-9) = 0;
end

emi40 = data(:,3); emi40(wavelength>1000e-9) = 0; % 4I9/2
emi41 = data(:,3); emi41(wavelength<1000e-9 | wavelength>1200e-9) = 0; % 4I11/2
emi42 = data(:,3); emi42(wavelength<1200e-9) = 0; % 4I13/2

y = [GSA011,GSA010,GSA09,GSA08,GSA07,GSA06,GSA05,GSA04,GSA03,ESA49,emi40,emi41,emi42];

figure;
plot(wavelength,y);

dlmwrite('Nd.txt',[wavelength,y],'delimiter',',','precision','%15.9g');