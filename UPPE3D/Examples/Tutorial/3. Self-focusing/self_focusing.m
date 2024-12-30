% This code simulates self-focusing of a Gaussian pulse.

close all; clearvars;

addpath('../../../UPPE3D algorithm/','../../../user_helpers/');

%% Setup fiber parameters
sim.lambda0 = 1030e-9; % the center wavelength
sim.gpuDevice.Index = 1;

% Load default parameters like
%
% load fiber.betas and fiber.SR based on your multimode folder above
% sim.f0 = 3e5/sim.lambda0; THz
% sim.save_period = 0; Save only the field at the beginning and the end fiber
% sim.scalar = true; Use scalar propagation
% sim.adaptive_dz.threshold = 1e-5; the threshold of the adaptive-step method
% sim.gpu_yes = true; Use GPU
% sim.include_Raman = 1; Consider the Raman
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.pulse_centering = true; Always shift the pulse to the center of the time window
% sim.gpuDevice.Index = 1; Use the GPU device 1
% sim.progress_bar = true; Show the progress bar
% sim.progress_bar_name = ''; Empty name for the progress bar
% sim.cuda_dir_path = 'UPPE3D/cuda'; Where the cuda files are
% ......
%
% Please check this function for details.
[fiber,sim] = load_default_UPPE3D_propagate([],sim); % load default parameters

fiber.L0 = 1e-2;
num_save = 100;
sim.save_period = fiber.L0/num_save;

%% Initial condition
MFD0 = 100e-6; % m
spatial_window = 1e-3; % m
tfwhm = 1; % ps
time_window = 10; % ps
energy = 9e3; % nJ
Nt = 2^8; % the number of time points
Nx = 2^8; % the number of spatial points
initial_condition = build_3Dgaussian(MFD0, spatial_window, tfwhm, time_window, energy, Nt, Nx);

%% Setup general parameters
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm
dx = spatial_window/Nx; % m
x = (-Nx/2:Nx/2-1)*dx;

%% Fiber waveguide structure
% Sellmeier coefficients
[a,b] = Sellmeier_coefficients(fiber.material);
% Calculate the index difference using the Sellmeier equation to generate n(lambda)
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));

fiber.n = n_from_Sellmeier(lambda/1e3).*ones(1,Nx,Nx);

%% Plot the initial field
% Show initial real space
figure;
pcolor(x*1e6,x*1e6,abs(squeeze(initial_condition.field(Nt/2,:,:))).^2); colormap(jet);colorbar;
shading interp;colormap(jet);colorbar;
xlabel('x (\mum)');
xlim([-1,1]*1e2);
ylim([-1,1]*1e2);
title('Initial real space');
% Show initial k space
figure;
pcolor(abs(fftshift(fft(fft(squeeze(initial_condition.field(Nt/2,:,:)),[],1),[],2))).^2); colormap(jet);colorbar;
shading interp;colormap(jet);colorbar;
title('Initial k space');

%% Propagate
prop_output = UPPE3D_propagate(fiber,initial_condition,sim);

%% Results
MFD = squeeze(calcMFD(squeeze(prop_output.field(Nt/2,:,:,:)),spatial_window))*1e6; % um
energy3D = squeeze(sum(abs(prop_output.field).^2,[1,2,3]));

%% Plot
% Show final real space
figure;
pcolor(x*1e6,x*1e6,abs(squeeze(prop_output.field(Nt/2,:,:,end))).^2); colormap(jet);colorbar;
shading interp;colormap(jet);colorbar;
xlabel('x (\mum)');
xlim([-1,1]*1e2);
ylim([-1,1]*1e2);
title('Final real space');
% Show final k space
figure;
pcolor(abs(fftshift(fft(fft(squeeze(prop_output.field(Nt/2,:,:,end)),[],1),[],2))).^2); colormap(jet);colorbar;
shading interp;colormap(jet);colorbar;
title('Final k space')

% Plot MFD
figure;
plot(prop_output.z*1e2,MFD,'linewidth',2);
xlabel('Propagation distance (cm)');
ylabel('MFD (\mum)');
set(gca,'fontsize',20);

% Movie
Frame = animator(prop_output.field,...
                 prop_output.z,MFD,...
                 Nt,dt,Nx,dx,lambda,...
                 fiber.L0);
implay(Frame(:),20);
exportVideo = VideoWriter('self-focusing');
exportVideo.FrameRate = 20;
open(exportVideo);
writeVideo(exportVideo, Frame(:));
close(exportVideo);