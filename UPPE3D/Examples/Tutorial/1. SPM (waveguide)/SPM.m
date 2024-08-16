% This code demonstrates the SPM .
%
% This code uses adaptive-step RK4IP for the passive-fiber propagation.

close all; clearvars;

addpath('../../UPPE3D algorithm/','../../user_helpers/');

lambda0 = 1030e-9; % m

%% Setup fiber parameters
sim.lambda0 = lambda0; % the center wavelength
s%im.gpu_yes = false;

% Load default parameters like
%
% load fiber.betas and fiber.SR based on your multimode folder above
% sim.f0 = 3e5/sim.lambda0; THz
% sim.save_period = 0; Save only the field at the beginning and the end fiber
% sim.ellipticity = 0; Linear polarization
% sim.scalar = true; Use scalar propagation
% sim.adaptive_dz.threshold = 1e-5; the threshold of the adaptive-step method
% sim.gpu_yes = true; Use GPU
% sim.Raman_model = 1; Use the isotropic Raman model
% sim.pulse_centering = true; Always shift the pulse to the center of the time window
% sim.num_photon_noise_per_band = 0; Don't include photon shot noise
% sim.gpuDevice.Index = 1; Use the GPU device 1
% sim.progress_bar = true; Show the progress bar
% sim.progress_bar_name = ''; Empty name for the progress bar
% sim.cuda_dir_path = 'UPPE3D/cuda'; Where the cuda files are
% ......
%
% Please check this function for details.
[fiber,sim] = load_default_UPPE3D_propagate([],sim); % load default parameters

fiber.L0 = 0.01;
num_save = 10;
sim.save_period = fiber.L0/num_save;

%% Setup general parameters
Nt = 2^8; % the number of time points
time_window = 0.8; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

%% Initial condition
Nx = 2^8;
spatial_window = 100; % um
fiber.n = calc_index_profile(fiber,Nx,spatial_window,f);
x = (-Nx/2:Nx/2-1)*spatial_window/Nx; % um

loaded_data = load('mode1wavelength10300.mat','phi','x','epsilon');

% take only the central part
[loaded_data.xx,loaded_data.yy] = meshgrid(loaded_data.x,loaded_data.x);
[xx,yy] = meshgrid(x,x);
phi = interp2(loaded_data.xx,loaded_data.yy,loaded_data.phi,xx,yy,'linear',0); % downsampling
phi = phi/sqrt(sum(abs(phi(:)).^2)*mean(diff(x*1e-6))^2);

dx = mean(diff(x));
Aeff = 1/(sum(phi(:).^4)*(dx*1e-6)^2); % mode area; m^2
MFR = sqrt(Aeff/pi)*1e6; % um

% input pulse
tfwhm = 0.03; % ps
total_energy = 10; % nJ
initial_pulse = build_MMgaussian(tfwhm, time_window, total_energy, 1, Nt);
initial_pulse.field = recompose_into_space(sim.gpu_yes,phi,initial_pulse.fields,sim.cuda_dir_path); initial_pulse = rmfield(initial_pulse,'fields');
initial_pulse.dx = mean(diff(x*1e-6)); initial_pulse.dy = initial_pulse.dx;

figure;
pcolor(abs(squeeze(initial_pulse.field(Nt/2,:,:))).^2); colormap(jet);colorbar;
shading interp;colormap(jet);colorbar;
figure;
pcolor(abs(fftshift(ifft(ifft(squeeze(initial_pulse.field(Nt/2,:,:)),[],1),[],2))).^2); colormap(jet);colorbar;
shading interp;colormap(jet);colorbar;

%% Propagate
prop_output = UPPE3D_propagate(fiber,initial_pulse,sim);

z = prop_output.z;

spectrum3D = abs(fftshift(ifft(prop_output.field,[],1))).^2;
spectrum = squeeze(spectrum3D(:,Nx/2,Nx/2,:));

output_field = zeros(Nt,1,num_save+1);
for i = 1:num_save+1
    output_field(:,:,i) = decompose_into_modes(sim.gpu_yes,phi,prop_output.field(:,:,:,i), prop_output.dx, sim.cuda_dir_path);
end

energy = squeeze(sum(abs(output_field).^2,1))*dt/1e3;
energy3D = squeeze(sum(abs(prop_output.field).^2,[1,2,3]));

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

%% Save the data
save('SPM3D.mat');