close all;

data1 = readmatrix('1.csv');
data1(data1(:,1)<600,2) = 0;
% Adjustment below is to make Judd-Ofelt's Omega makes sense
lambda = data1(data1(:,1)>730 & data1(:,1)<940,1);
disp(trapz(lambda,data1(data1(:,1)>730 & data1(:,1)<940,2)));
lambda0 = 786.2;
tmpp = data1(data1(:,1)>730 & data1(:,1)<940,2); tmpp = tmpp/max(tmpp); x = 100;
tmpp = 1.2*tmpp.^((x+0.1)-x*exp(-(lambda-lambda0).^2/150^2));
data1(data1(:,1)>730 & data1(:,1)<940,2) = data1(data1(:,1)>730 & data1(:,1)<940,2).*tmpp;
disp(trapz(lambda,data1(data1(:,1)>730 & data1(:,1)<940,2)));

lambda = data1(data1(:,1)>1380,1);
%disp(trapz(lambda,data1(data1(:,1)>1380,2)));
lambda0 = 1635.68;
tmpp = data1(data1(:,1)>1380,2); tmpp = tmpp/max(tmpp);
tmpp = tmpp.^(-exp(-(lambda-lambda0).^2/80^2));
data1(data1(:,1)>1380,2) = data1(data1(:,1)>1380,2).*tmpp;
%disp(trapz(lambda,data1(data1(:,1)>1380,2)));

data5 = readmatrix('5.csv');
data5 = [data5; [463.5,0]];
data5(data5(:,1)<428.3,2) = 0;
data1 = [data5;data1];

% data1 is GSA
% Separate GSA into each transition
data1a = data1(data1(:,1)<=520,                   :);
data1b = data1(data1(:,1)>520  & data1(:,1)<=740, :);
data1c = data1(data1(:,1)>740  & data1(:,1)<=875, :);
data1d = data1(data1(:,1)>875  & data1(:,1)<=1400,:);
data1e = data1(data1(:,1)>1400,                   :);

data2 = readmatrix('2.csv');

data3 = readmatrix('3.csv');

data4 = readmatrix('4.csv');
data13 = readmatrix('thorlabs_e.csv');
cc = max(data4(:,2))/max(data13(:,2));
data13(:,2) = smooth(data13(:,2)*cc,10);
data12 = readmatrix('thorlabs_a.csv');
data12(:,2) = smooth(smooth(data12(:,2),60),60)*cc;
data12(data12(:,1)<1500,2) = data12(data12(:,1)<1500,2).*exp(linspace(-4,0,sum(data12(:,1)<1500))');
data12(data12(:,1)>1961,2) = data12(data12(:,1)>1961,2).*exp(linspace(0,-4,sum(data12(:,1)>1961))');

data6 = readmatrix('6.csv');

data7 = readmatrix('7.csv');

data8 = readmatrix('8.csv');

data9 = readmatrix('9.csv');

data10 = readmatrix('10.csv');

data11 = readmatrix('11.csv');

x = linspace(min([data1(1,1),data2(1,1),data3(1,1),data4(1,1),data5(1,1),data6(1,1),data7(1,1),data8(1,1),data9(1,1),data10(1,1),data11(1,1)]),...
             max([data1(end,1),data2(end,1),data3(end,1),data4(end,1),data5(end,1),data6(end,1),data7(end,1),data8(end,1),data9(end,1),data10(end,1),data11(end,1)]),3000)';
y1a = interp1(data1a(:,1),data1a(:,2),x,'linear',0)*1e-25;
y1b = interp1(data1b(:,1),data1b(:,2),x,'linear',0)*1e-25;
y1c = interp1(data1c(:,1),data1c(:,2),x,'linear',0)*1e-25;
y1d = interp1(data1d(:,1),data1d(:,2),x,'linear',0)*1e-25;
y1e = interp1(data1e(:,1),data1e(:,2),x,'linear',0)*1e-25;
y2 = interp1(data2(:,1),data2(:,2),x,'linear',0)*1e-25;
y3 = interp1(data3(:,1),data3(:,2),x,'linear',0)*1e-26;
y4 = interp1(data4(:,1),data4(:,2),x,'linear',0)*1e-25;
y6 = interp1(data6(:,1),data6(:,2),x,'linear',0)*1e-25;
y7 = interp1(data7(:,1),data7(:,2),x,'linear',0)*1e-25;
y8 = interp1(data8(:,1),data8(:,2),x,'linear',0)*1e-25;
y9 = interp1(data9(:,1),data9(:,2),x,'linear',0);
y10 = interp1(data10(:,1),data10(:,2),x,'linear',0)*1e-24;
y11 = interp1(data11(:,1),data11(:,2),x,'linear',0)*1e-24;
y12 = interp1(data12(:,1),data12(:,2),x,'linear',0)*1e-25;
y13 = interp1(data13(:,1),data13(:,2),x,'linear',0)*1e-25;
y = [y1a,y1b,y1c,y1d,y12,y2,y3,y8,y10,y6,y9,y7,y13,y11]; y(y<0)=0;

figure;plot(x,y);

dlmwrite('Tm.txt',[x*1e-9,y],'delimiter',',','precision','%15.9g');