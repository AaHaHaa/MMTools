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
[x,y] = meshgrid(f-sim.f0,prop_output.z(2:end));
pcolor(x,y,permute(abs(fftshift(ifft(output_field(:,1,2:end)),1)).^2,[3 1 2]));
xlim([-100,100]);
shading interp; colormap(jet);
xlabel('\nu-\nu_0');
ylabel('z');
title('Spectrum during propagation');
set(gca,'fontsize',14);

% The final spatial profile at the peak power
figure;
pcolor(xx,yy,abs(squeeze(prop_output.field(Nt/2,:,:,end))).^2);
shading interp;colormap(jet);colorbar;
xlabel('Length (\mum)');
ylabel('Length (\mum)');
set(gca,'fontsize',14);

figure;
pcolor(xx,yy,abs(fftshift(ifft(ifft(squeeze(prop_output.field(Nt/2,:,:,end)),[],1),[],2))).^2);
shading interp;colormap(jet);colorbar;
set(gca,'fontsize',14);