% This script verifies the Hankel transform scheme, FHATHA, with a Gaussian
% function.
% Furthermore, it applies the inverse Hankel transform to the transformed
% signal to see whether the signal can be recovered.
%
% The signal in this example needs to drop to zero at the window edge;
% otherwise, there will be aliasing, as in DFT. Hankel transform results
% from 2D spatial Fourier transform.

clearvars; close all;

Nr = 500;

%% FHATHA
Nf = 10; % =2*pi*k_max*r_max
r_max = 1;
k_max = 2*pi*Nf/r_max;

[r,kr,...
 l0,exp_prefactor,...
 Q] = Hankel_info(Nr,r_max,k_max);

f = exp(-20*r.^2);

f_H = FHATHA(f,...
             r_max,kr,...
             l0,exp_prefactor,...
             Q);

figure;
plot(kr,real(f_H)*2*pi,'linewidth',2,'Color','b')

%% Analytic solution of the Hankel transform of f
eta = kr*r_max;
F = pi/20*exp(-kr.^2/80);

hold on;
plot(kr(1:5:end),F(1:5:end),'o','linewidth',1,'Color','r'); % reduce the sampling rate for visualization
hold off;
legend('FHATHA','Analytic');
set(gca,'fontsize',20);
xlabel('k_r');
print(gcf,'FHATHA_verification2.pdf','-dpdf');

%% Check: inverse Hankel transform
invf = FHATHA(f_H,...
              k_max,r,...
              l0,exp_prefactor,...
              Q);

figure;
plot(r,real(invf),'linewidth',2,'Color','b');
hold on;
plot(r(1:5:end),f(1:5:end),'o','linewidth',1,'Color','r');
hold off;
legend('inverse FHATHA','Analytic');
set(gca,'fontsize',20);
xlabel('r');
print(gcf,'invFHATHA_verification2.pdf','-dpdf');