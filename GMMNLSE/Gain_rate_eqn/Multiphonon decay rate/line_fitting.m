function [ linefitcoeff,corr ] = line_fitting( x,y )
%LINE_FITTING Summary of this function goes here
%   Detailed explanation goes here

linefitcoeff = polyfit(x,y,1);
yval = polyval(linefitcoeff,x);

plot(x,y,'b',x,yval,'r');
legend('data','fitted line','Location','southeast');

corr = corrcoef(y,yval);

end

