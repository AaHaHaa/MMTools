% This script verifies the Hankel transform scheme, FHATHA, by duplicating
% Fig.1 of the following paper:
%
%   Magni et al., "High-accuracy fast Hankel transform for optical beam
%   propagation," J. Opt. Soc. Am. A 9(11), 2031-2033 (1992)
%
% The function is f = sqrt(5/2/pi)*r.^2.

clearvars; close all;

Nr = 500;

%% FHATHA
Nf = 10; % =2*pi*k_max*r_max
r_max = 1;
k_max = 2*pi*Nf/r_max;

[r,kr,...
 l0,exp_prefactor,...
 Q] = Hankel_info(Nr,r_max,k_max);

f = sqrt(5/2/pi)*r.^2;

f_H = FHATHA(f,...
             r_max,kr,...
             l0,exp_prefactor,...
             Q);

figure;
plot(kr,real(f_H)*2*pi,'linewidth',2,'Color','b')

%% Analytic solution of the Hankel transform of f
eta = kr*r_max;
F = sqrt(10*pi)*eta.^(-4).*(2*eta.^2.*besselj(0,eta) + (eta.^3-4*eta).*besselj(1,eta));

hold on;
plot(kr(1:5:end),F(1:5:end),'o','linewidth',1,'Color','r'); % reduce the sampling rate for visualization
hold off;
legend('FHATHA','Analytic');
set(gca,'fontsize',20);
xlabel('k_r');
print(gcf,'FHATHA_verification1.pdf','-dpdf');