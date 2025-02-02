% This script verifies the Hankel transform scheme, FHATHA, with a Gaussian
% function.
% Furthermore, it applies the inverse Hankel transform to the transformed
% signal to see whether the signal can be recovered.
%
% The signal in this example needs to drop to zero at the window edge;
% otherwise, there will be aliasing, as in DFT. Hankel transform results
% from 2D spatial Fourier transform.

clearvars; close all;

addpath('../');

Nr = 1000;

%% FHATHA
Nf = 10; % =2*pi*k_max*r_max
r_max = 1;
k_max = 2*pi*Nf/r_max;

[r,kr,...
 dr,dkr,...
 l0,exp_prefactor,r2_prefactor,...
 ifftQ] = Hankel_info(Nr,r_max,k_max);

f = exp(-20*r.^2);

f_H = FHATHA(f,...
             r_max,...
             r,kr,...
             dr,dkr,...
             l0,exp_prefactor,r2_prefactor,...
             ifftQ);

figure;
plot(kr,real(f_H)*2*pi,'linewidth',2,'Color','b')

%% Analytic solution of the Hankel transform of f
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
              k_max,...
              kr,r,...
              dkr,dr,...
              l0,exp_prefactor,r2_prefactor,...
              ifftQ);

figure;
plot(r,real(invf),'linewidth',2,'Color','b');
hold on;
plot(r(1:5:end),f(1:5:end),'o','linewidth',1,'Color','r');
hold off;
legend('inverse FHATHA','Analytic');
set(gca,'fontsize',20);
xlabel('r');
print(gcf,'invFHATHA_verification2.pdf','-dpdf');

%% Check many iterations
num = 100;
diff_E = zeros(num,1);
NRMSE = zeros(num,1);
for i = 1:num
    f_H = FHATHA(invf,...
                 r_max,...
                 r,kr,...
                 dr,dkr,...
                 l0,exp_prefactor,r2_prefactor,...
                 ifftQ);
    invf = FHATHA(f_H,...
                  k_max,...
                  kr,r,...
                  dkr,dr,...
                  l0,exp_prefactor,r2_prefactor,...
                  ifftQ);
    
    diff_E(i) = 2*pi*abs(trapz(r,abs(f).^2.*r) - trapz(r,abs(invf).^2.*r));
    NRMSE(i) = sqrt(sum(abs(f-invf).^2)/(Nr+1));
end

figure;
plot(r,real(invf),'linewidth',2,'Color','b');
hold on;
plot(r(1:5:end),f(1:5:end),'o','linewidth',1,'Color','r');
hold off;
legend('inverse FHATHA','Analytic');
set(gca,'fontsize',20);
xlabel('r');
print(gcf,'repeated_FHATHA.pdf','-dpdf');

figure;
yyaxis left;
plot(diff_E,'linewidth',2,'Color','b');
ylabel('Power deviation');
set(gca,'fontsize',20,'YColor','b');
yyaxis right;
plot(NRMSE,'linewidth',2,'Color',[0.8510,0.3255,0.0980]);
ylabel('RMSE');
xlabel('Iteration');
set(gca,'fontsize',20,'YColor',[0.8510,0.3255,0.0980]);
print(gcf,'repeated_FHATHA_error.pdf','-dpdf');