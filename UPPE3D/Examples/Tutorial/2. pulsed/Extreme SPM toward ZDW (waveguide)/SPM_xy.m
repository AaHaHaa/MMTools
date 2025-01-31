% This code demonstrates the SPM toward the zero-dispersion wavelength 
% around 1300 nm in silica.
%
% This script employs the 3D-UPPE that uses full x-y dimension. For
% more-efficient modeling, pelase see its radially-symmetric version.
%
% Compared to the radially-symmetric script, the spectrum generated with 
% this script is more-distorted, if we use the GNLSE result as a reference
% point. This is due to the coarse k-space sampling in this example. To
% resolve this, spatial_window needs to be increase, as well as the spatial
% sampling point to maintain the spatial sampling density. However, this 
% overwhelms our 12-GB Nvidia Titan XP GPU. As a result, this example
% demonstrates the general difficulty in running a full-field simulation.
% CPU computations with RAM can increase the memory limit but the 
% simulation is order-of-magnitude slower than with GPU computations.

close all; clearvars;

addpath('../../../../UPPE3D algorithm/','../../../../user_helpers/');

% Add the path for the fiber profile
addpath('../../../../Fibers/Single-mode-fiber profile for examples/');

lambda0 = 1030e-9; % m

%% Setup fiber parameters
sim.lambda0 = lambda0; % the center wavelength
%im.gpu_yes = false;
sim.include_Raman = false;

% Load default parameters like
%
% load fiber.betas and fiber.SR based on your multimode folder above
% sim.f0 = 3e5/sim.lambda0; THz
% sim.save_period = 0; Save only the field at the beginning and the end fiber
% sim.ellipticity = 0; Linear polarization
% sim.scalar = true; Use scalar propagation
% sim.adaptive_dz.threshold = 1e-5; the threshold of the adaptive-step method
% sim.gpu_yes = true; Use GPU
% sim.include_Raman = 1; Consider the Raman
% sim.pulse_centering = true; Always shift the pulse to the center of the time window
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
Nx = 2^7;
spatial_window = 50; % um
fiber.n = calc_index_profile_xy(fiber,Nx,spatial_window,f);
x = (-Nx/2:Nx/2-1)*spatial_window/Nx; % um
kx = 2*pi*(-Nx/2:Nx/2-1)/spatial_window; % 2*pi/um

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
initial_condition = build_MMgaussian(tfwhm, time_window, total_energy, 1, Nt);
initial_condition.field = recompose_into_space(sim.gpu_yes,phi,initial_condition.fields,sim.cuda_dir_path); initial_condition = rmfield(initial_condition,'fields');
initial_condition.dx = mean(diff(x*1e-6)); initial_condition.dy = initial_condition.dx;

%% Plot the initial field
% Show initial real space
figure;
pcolor(x,x,abs(squeeze(initial_condition.field(Nt/2,:,:))).^2); colormap(jet);
shading interp;colormap(jet);colorbar;
xlabel('x (\mum)');
ylabel('y (\mum)');
xlim([-1,1]*10);
ylim([-1,1]*10);
daspect([1 1 1]); % make aspect ratio = 1
set(gca,'fontsize',20);
title('Initial real space');
% Show initial k space
figure;
pcolor(kx,kx,abs(fftshift(fft(fft(squeeze(initial_condition.field(Nt/2,:,:)),[],1),[],2))).^2); colormap(jet);
shading interp;colormap(jet);colorbar;
xlabel('k_x (2\pi/\mum)');
ylabel('k_y (2\pi/\mum)');
daspect([1 1 1]); % make aspect ratio = 1
set(gca,'fontsize',20);
title('Initial k space');

%% Propagate
prop_output = UPPE3D_propagate(fiber,initial_condition,sim);

output_field = zeros(Nt,1,num_save+1);
for i = 1:num_save+1
    output_field(:,:,i) = decompose_into_modes(sim.gpu_yes,phi,prop_output.field(:,:,:,i), prop_output.dx, sim.cuda_dir_path);
end

%% Results
MFD = squeeze(calcMFD_xy(squeeze(prop_output.field(floor(Nt/2)+1,:,:,:)),spatial_window*1e-6))*1e6; % um
energy = squeeze(sum(abs(prop_output.field).^2,[1,2,3]))*(dx^2*1e-12)*dt/1e3; % nJ

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
[x_t,y_z] = meshgrid(t,prop_output.z);
pcolor(x_t,y_z,permute(abs(output_field(:,1,:)).^2,[3 1 2]));
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

% MFD evolution
figure;
plot(prop_output.z,MFD,'linewidth',2,'Color','b');
xlabel('Propagation distance (m)');
ylabel('MFD (\mum)');
set(gca,'fontsize',20);

% Energy
figure;
plot(prop_output.z,energy,'linewidth',2,'Color','b');
xlabel('Propagation distance (m)');
ylabel('Power (nJ)');
set(gca,'fontsize',20);

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

%% Save the data
save('SPM3D.mat');