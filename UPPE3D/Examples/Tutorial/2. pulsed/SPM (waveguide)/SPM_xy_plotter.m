close all; clearvars;

load('SPM3D.mat');

%% Plot
% Time
figure;
h = plot(t,abs(output_field(:,:,end)).^2);
xlim([-0.2,0.2]);
xlabel('t');
ylabel('Power');
title('Field');
set(h,'linewidth',2);
set(gca,'fontsize',14);

% Spectrum
figure;
h = plot(f-sim.f0,abs(fftshift(ifft(output_field(:,:,end)),1)).^2);
xlim([-100,100]);
xlabel('\nu-\nu_0');
ylabel('PSD');
title('Spectrum');
set(h,'linewidth',2);
set(gca,'fontsize',14);

% Comparison of time
figure;
[x,y] = meshgrid(t,prop_output.z);
pcolor(x,y,permute(abs(output_field(:,1,:)).^2,[3 1 2]));
shading interp; colormap(jet);
xlim([-0.2,0.2]);
xlabel('t');
ylabel('z');
title('Field during propagation');
set(gca,'fontsize',14);

% Comparison of spectra
figure;
[x_f,y_z] = meshgrid(f-sim.f0,prop_output.z(2:end));
pcolor(x_f,y_z,permute(abs(fftshift(ifft(output_field(:,1,2:end)),1)).^2,[3 1 2]));
xlim([-100,100]);
shading interp; colormap(jet);
xlabel('\nu-\nu_0');
ylabel('z');
title('Spectrum during propagation');
set(gca,'fontsize',14);

% Show final real space
figure;
pcolor(x,x,abs(squeeze(prop_output.field(Nt/2,:,:,end))).^2); colormap(jet);colorbar;
shading interp;colormap(jet);colorbar;
xlabel('x (\mum)');
ylabel('y (\mum)');
xlim([-1,1]*10);
ylim([-1,1]*10);
daspect([1 1 1]); % make aspect ratio = 1
set(gca,'fontsize',20);
title('Final real space');
% Show final k space
figure;
pcolor(kx,kx,abs(fftshift(fft(fft(squeeze(prop_output.field(Nt/2,:,:,end)),[],1),[],2))).^2); colormap(jet);colorbar;
shading interp;colormap(jet);colorbar;
xlabel('k_x (2\pi/\mum)');
ylabel('k_y (2\pi/\mum)');
daspect([1 1 1]); % make aspect ratio = 1
set(gca,'fontsize',20);
title('Final k space')