close all; clearvars;

load('Self_steepening.mat');

%% Plot
% Time
figure;
plot(t,abs(prop_output.fields(:,:,1)).^2,'linewidth',2);
hold on;
plot(t,abs(prop_output.fields(:,:,end)).^2,'linewidth',2);
hold off;
legend('Initial','Final');
xlim([-3,3]);
xlabel('t');
ylabel('Power');
title('Field');
set(gca,'fontsize',20);

% Spectrum
figure;
plot(f-sim.f0,abs(fftshift(ifft(prop_output.fields(:,:,end)),1)).^2,'linewidth',2);
xlim([-100,100]);
xlabel('\nu-\nu_0');
ylabel('PSD');
title('Spectrum');
set(gca,'fontsize',20);

% Comparison of time
figure;
[x,y] = meshgrid(t,prop_output.z);
pcolor(x,y,permute(abs(prop_output.fields(:,1,:)).^2,[3 1 2]));
shading interp; colormap(jet);
xlim([-3,3]);
xlabel('t');
ylabel('z');
title('Field during propagation');
set(gca,'fontsize',20);

% Comparison of spectra
figure;
[x,y] = meshgrid(f-sim.f0,prop_output.z(2:end));
pcolor(x,y,permute(abs(fftshift(ifft(prop_output.fields(:,1,2:end)),1)).^2,[3 1 2]));
xlim([-100,100]);
shading interp; colormap(jet);
xlabel('\nu-\nu_0');
ylabel('z/LD');
title('Spectrum during propagation');
set(gca,'fontsize',20);

% Movie
implay(Frame,20);

exportVideo = VideoWriter('Self_steepening');
exportVideo.FrameRate = 20;
open(exportVideo);
writeVideo(exportVideo, Frame);
close(exportVideo);