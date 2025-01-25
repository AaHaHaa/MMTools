close all; clearvars;

data_r = load('self_focusing_r.mat');
data_xy = load('self_focusing_xy.mat');

% Plot MFD
figure;
plot(data_r.prop_output.z*1e2,data_r.MFD,'linewidth',2,'Color','b');
hold on;
plot(data_xy.prop_output.z*1e2,data_xy.MFD,'linewidth',2,'Color','r');
hold off;
xlabel('Propagation distance (cm)');
ylabel('MFD (μm)');
set(gca,'fontsize',20);
legend('Hankel','Fourier');
print(gcf,'MFD.pdf','-dpdf');

% Plot the field in real xy-space
figure;
plot(data_r.r*1e6,abs(squeeze(data_r.prop_output.field(ceil(data_r.Nt/2),:,end))).^2,'linewidth',2,'Color','b');
hold on;
plot(data_xy.x*1e6,abs(squeeze(data_xy.prop_output.field(ceil(data_xy.Nt/2),:,ceil(data_xy.Nt/2),end))).^2,'linewidth',2,'Color','r');
plot(data_xy.x*1e6,abs(squeeze(data_xy.prop_output.field(ceil(data_xy.Nt/2),:,ceil(data_xy.Nt/2),end))).^2,'o','Color','r');
hold off;
xlabel('r (μm)');
xlim([0,40]);
ylabel('Intensity (W/m^2)');
set(gca,'fontsize',20);
legend('Hankel','Fourier');
print(gcf,'I.pdf','-dpdf');
