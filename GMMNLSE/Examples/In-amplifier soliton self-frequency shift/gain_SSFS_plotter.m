close all; clearvars;

load('SSFS.mat');

figure;
yyaxis left;
plot(all_duration,all_lambda1,'linewidth',2,'Color','b');
set(gca,'YColor','b');
xlabel('Chirped duration (fs)');
ylabel('Wavelength (nm)');
set(gca,'fontsize',20);
yyaxis right;
plot(all_duration,all_energy,'linewidth',2,'Color',[0.851,0.3255,0.098]);
ylabel('Energy (nJ)');