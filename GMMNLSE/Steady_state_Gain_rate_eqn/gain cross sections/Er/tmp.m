close all;

data1 = readmatrix('Er (optiwave and Smirnov).xlsx');
lambda = data1(:,1)*1e9;
xa15 = data1(lambda<950,1)*1e9;
ya15 = data1(lambda<950,4);
xa03 = data1(lambda<865,1)*1e9;
ya03 = data1(lambda<865,2);
xa02 = data1(lambda>865 & lambda<1100,1)*1e9;
ya02 = data1(lambda>865 & lambda<1100,2);
xa01 = data1(lambda>1200,1)*1e9;
ya01 = data1(lambda>1200,2);
xe01 = data1(lambda>1200,1)*1e9;
ye01 = data1(lambda>1200,3);

data3 = readmatrix('absorption_ErZSG.csv');
lambda = data3(:,1);
data3(:,2) = data3(:,2)*4.802e-25/1.35496;
xa02 = data3(lambda>900 & lambda<1100,1);
ya02 = data3(lambda>900 & lambda<1100,2);
xa03 = data3(lambda>750 & lambda<900,1);
ya03 = data3(lambda>750 & lambda<900,2);
xa04 = data3(lambda>600 & lambda<750,1);
ya04 = data3(lambda>600 & lambda<750,2);
xa05 = data3(lambda>505 & lambda<600,1);
data3(lambda>520.2 & lambda<523.6,2) = smooth(data3(lambda>520.2 & lambda<523.6,2));
ya05 = data3(lambda>505 & lambda<600,2);
xa06 = data3(lambda>470 & lambda<505,1);
ya06 = data3(lambda>470 & lambda<505,2);
xa07 = data3(lambda>427.2 & lambda<470,1);
ya07 = data3(lambda>427.2 & lambda<470,2);
xa08 = data3(lambda<427.2,1);
ya08 = data3(lambda<427.2,2);

data2 = readmatrix('esa980.csv');
xa26 = data2(:,1);
ya26 = data2(:,2)*1e-24*max(ya02)/1.70532e-25;

data4 = readmatrix('2700a.csv');
xa12 = data4(:,1);
ya12 = data4(:,2)*1e-24;
data4 = readmatrix('2700e.csv');
xe21 = data4(:,1);
ye21 = data4(:,2)*1e-24;

data5 = readmatrix('4400a.csv');
xa23 = data5(:,1);
ya23 = data5(:,2)*1e-24;
data5 = readmatrix('4400e.csv');
xe32 = data5(:,1);
ye32 = data5(:,2)*1e-24;

x = linspace(min(350,min([data3(:,1);xa15;xa03;xa02;xa01;xe01;xa26;xa12;xe21;xa23;xe32])),max([xa15;xa03;xa02;xa01;xe01;xa26;xa12;xe21;xa23;xe32]),1000)';
ya15 = interp1(xa15,ya15,x,'linear',0);
ya08 = interp1(xa08,ya08,x,'linear',0);
ya07 = interp1(xa07,ya07,x,'linear',0);
ya06 = interp1(xa06,ya06,x,'linear',0);
ya05 = interp1(xa05,ya05,x,'linear',0);
ya04 = interp1(xa04,ya04,x,'linear',0);
ya03 = interp1(xa03,ya03,x,'linear',0);
ya02 = interp1(xa02,ya02,x,'linear',0);
ya01 = interp1(xa01,ya01,x,'linear',0);
ye01 = interp1(xe01,ye01,x,'linear',0);
ya26 = interp1(xa26,ya26,x,'linear',0);
ya12 = interp1(xa12,ya12,x,'linear',0);
ye21 = interp1(xe21,ye21,x,'linear',0);
ya23 = interp1(xa23,ya23,x,'linear',0);
ye32 = interp1(xe32,ye32,x,'linear',0);
y = [ya08,ya07,ya06,ya05,ya04,ya03,ya02,ya01,ya15,ya26,ya12,ya23,ye01,ye21,ye32];

y(y<0) = 0;

figure;
plot(x,y);

dlmwrite('Er.txt',[x*1e-9,y],'delimiter',',','precision','%15.9g');