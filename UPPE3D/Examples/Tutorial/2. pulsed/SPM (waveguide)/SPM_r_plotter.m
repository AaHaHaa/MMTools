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
plot(r*1e6,abs(squeeze(prop_output.field(ceil(Nt/2),:,end))).^2,'linewidth',2,'Color','b');
xlabel('r (\mum)');
xlim([0,10]);
set(gca,'fontsize',20);
title('Final real space');
% Plot the 2D field with pcolor
% However, the Hankel transform doesn't sample at the origin r=0, so we
% need to find it first. This can be done with Hankel_f_at_0().
A0 = Hankel_f_at_0(prop_output.field(ceil(Nt/2),:,end),l0);
radialPcolor([0,r]*1e6,cat(2,abs(A0).^2,abs(squeeze(prop_output.field(ceil(Nt/2),:,end))).^2));
xlabel('x (\mum)');
ylabel('y (\mum)');
xlim([-10,10]);
ylim([-10,10]);
set(gca,'fontsize',20);
daspect([1 1 1]); % make aspect ratio = 1
title('Final real space');
% Show final k space
A_H = 2*pi*FHATHA(squeeze(prop_output.field(ceil(Nt/2),:,end)),...
                  r_max,kr,...
                  l0,exp_prefactor,...
                  Q);
figure;
plot(kr/1e6,abs(A_H).^2,'linewidth',2,'Color','r');
xlabel('k_r (2\pi/\mum)');
set(gca,'fontsize',20);
title('Final k space');
% Plot the 2D field with pcolor
A_H0 = Hankel_f_at_0(A_H,l0);
radialPcolor([0,kr]/1e6,cat(2,abs(A_H0).^2,abs(A_H).^2));
xlabel('k_x (2\pi/\mum)');
ylabel('k_y (2\pi/\mum)');
set(gca,'fontsize',20);
daspect([1 1 1]); % make aspect ratio = 1
title('Final k space');