close all;

data1 = readmatrix('a.csv');
xa01 = data1(:,1);
ya01 = data1(:,2);
data2 = readmatrix('e.csv');
xe10 = data2(:,1);
ye10 = data2(:,2);

data3 = readmatrix('GSA.csv');
w = data3(:,1);
data3(:,2) = data3(:,2)/0.679205*0.291508;
ya04 = data3(:,2); ya04(w<0.5815 | w>0.746) = 0;
ya03 = data3(:,2); ya03(w<0.746 | w>1) = 0;
ya02 = data3(:,2); ya02(w<1 | w>1.4) = 0;

data4 = readmatrix('ESA.csv');
ww = data4(:,1);
data4(:,2) = data4(:,2)/0.76894*0.0988206;

w = w*1e3;
x = linspace(min([w;xa01;xe10]),max([w;xa01;xe10]),5000)';
ya01 = smooth(interp1(xa01,ya01,x,'linear',0));
ya02 = smooth(smooth(interp1(w,ya02,x,'linear',0),20),20);
ya03 = smooth(smooth(interp1(w,ya03,x,'linear',0),20),20);
ya04 = smooth(smooth(interp1(w,ya04,x,'linear',0),20),20);
ya14 = smooth(smooth(interp1(ww,data4(:,2),x,'linear',0),20),20);
ye10 = smooth(smooth(interp1(xe10,ye10,x,'linear',0),20),20);

y = [ya04,ya03,ya02,ya01,ya14,ye10]*1e-24;

y(y<0) = 0;

figure;
plot(x,y);

dlmwrite('Ho.txt',[x*1e-9,y],'delimiter',',','precision','%15.9g');