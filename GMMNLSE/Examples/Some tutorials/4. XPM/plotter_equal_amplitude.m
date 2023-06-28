% This code plots the result from "vector_soliton_equal_amplitude".

close all; clearvars;

load('vector_soliton_equal_amplitude.mat');

% Time
figure;
plot(t/T0,abs(prop_output.fields(:,:,end)).^2,'linewidth',2);
xlim([-6 6]);
xlabel('t/T_0');
ylabel('Intensity');

% Spectrum
figure;
plot((f-sim.f0)*T0,abs(fftshift(ifft(prop_output.fields(:,:,end)),1)).^2,'linewidth',2);
xlim([-0.25 0.25]);
xlabel('(\nu-\nu_0)*T_0');
ylabel('Intensity');

% Comparison of time
figure;
subplot(2,1,1);
[x,y] = meshgrid(t/T0,prop_output.z/LD);
pcolor(x,y,permute(abs(prop_output.fields(:,1,:)).^2,[3 1 2]));
shading interp; colormap(jet);
xlim([-10 10]);
xlabel('t/T_0');
ylabel('z/LD');

subplot(2,1,2);
[x,y] = meshgrid(t/T0,prop_output.z/LD);
pcolor(x,y,permute(abs(prop_output.fields(:,2,:)).^2,[3 1 2]));
shading interp; colormap(jet);
xlim([-10 10]);
xlabel('t/T_0');
ylabel('z/LD');

% Comparison of spectra
figure;
subplot(2,1,1);
[x,y] = meshgrid((f-sim.f0)*T0,prop_output.z/LD);
pcolor(x,y,permute(abs(fftshift(ifft(prop_output.fields(:,1,:)),1)).^2,[3 1 2]));
shading interp; colormap(jet);
xlim([-0.5 0.5]);
xlabel('(\nu-\nu_0)*T_0');
ylabel('z/LD');

subplot(2,1,2);
[x,y] = meshgrid((f-sim.f0)*T0,prop_output.z/LD);
pcolor(x,y,permute(abs(fftshift(ifft(prop_output.fields(:,2,:)),1)).^2,[3 1 2]));
shading interp; colormap(jet);
xlim([-0.5 0.5]);
xlabel('(\nu-\nu_0)*T_0');
ylabel('z/LD');