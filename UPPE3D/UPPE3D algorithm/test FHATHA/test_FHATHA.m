% This script verifies the Hankel transform scheme, FHATHA, by duplicating
% Fig.1 of the following paper:
%
%   Magni et al., "High-accuracy fast Hankel transform for optical beam
%   propagation," J. Opt. Soc. Am. A 9(11), 2031-2033 (1992)
%
% The function is f = sqrt(5/2/pi)*r.^2.

clearvars; close all;

addpath('../');

Nr = 500;

%% FHATHA
Nf = 10; % =2*pi*k_max*r_max
r_max = 1;
k_max = 2*pi*Nf/r_max;

[r,kr,...
 dr,dkr,...
 exp_prefactor,r2_prefactor,...
 ifftQ] = Hankel_info(Nr,r_max,k_max);

f = sqrt(5/2/pi)*r.^2;

f_H = FHATHA(f,...
             r_max,...
             r,kr,...
             dr,dkr,...
             exp_prefactor,r2_prefactor,...
             ifftQ);

figure;
plot(kr,real(f_H)*2*pi,'linewidth',2,'Color','b')

%% Analytic solution of the Hankel transform of f
kr = [linspace(kr(1),kr(2),100),kr(3:end)];
eta = kr*r_max;
F = sqrt(10*pi)*eta.^(-3).*(2*eta.*besselj(0,eta) + (eta.^2-4).*besselj(1,eta));

hold on;
plot(kr(1:5:end),F(1:5:end),'o','linewidth',1,'Color','r'); % reduce the sampling rate for visualization
hold off;
legend('FHATHA','Analytic');
set(gca,'fontsize',20);
xlabel('k_r');
print(gcf,'FHATHA_verification1.pdf','-dpdf');