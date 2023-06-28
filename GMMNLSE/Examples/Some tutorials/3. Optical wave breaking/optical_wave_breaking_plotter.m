close all; clearvars;

load('optical_wave_breaking.mat');

%% Plot
% Time
figure;
h = plot(t/T0,abs(prop_output.fields(:,:,end)).^2);
xlim([-7 7]);
xlabel('t/T_0');
ylabel('Intensity');
title('Field');
set(h,'linewidth',2);
set(gca,'fontsize',14);

% Spectrum
figure;
h = plot((f-sim.f0)*T0,abs(fftshift(ifft(prop_output.fields(:,:,end)),1)).^2);
xlim([-20 20]);
xlabel('(\nu-\nu_0)*T_0');
ylabel('Intensity');
title('Spectrum');
set(h,'linewidth',2);
set(gca,'fontsize',14);

% Comparison of time
figure;
[x,y] = meshgrid(t/T0,z/LD);
pcolor(x,y,permute(abs(prop_output.fields(:,1,:)).^2,[3 1 2]));
shading interp;
xlim([-50 50]);
xlabel('t/T_0');
ylabel('z/LD');
title('Field during propagation');
set(gca,'fontsize',14);

% Comparison of spectra
figure;
[x,y] = meshgrid((f-sim.f0)*T0,z(2:end)/LD);
pcolor(x,y,permute(abs(fftshift(ifft(prop_output.fields(:,1,2:end)),1)).^2,[3 1 2]));
xlim([-20 20]);
shading interp;
xlabel('(\nu-\nu_0)*T_0');
ylabel('z/LD');
title('Spectrum during propagation');
set(gca,'fontsize',14);

% Movie
implay(Frame,20);