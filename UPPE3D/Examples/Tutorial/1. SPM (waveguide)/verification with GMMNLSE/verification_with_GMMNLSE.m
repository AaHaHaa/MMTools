% This code demonstrates the SPM for verification of the 3D-UPPE.
%
% This code won't run here. It needs to be in my GMMNLSE folder to work.
% I put it here simply for reader's information.
%
% This code uses adaptive-step RK4IP for the passive-fiber propagation.
%
%
% Highly recommend the user to read the following article to understand how
% optical wave breaking evolves. It contains detailed and good explanation.
% For example, the spectral peaks at nu-nu0 = +-15 result from the 
% four-wave mixing at the interference. The energy transfers from the pulse
% central frequency to the exterior part of the spectrum. See this clearly 
% under the log plot by setting "log_yes = true" in this code.
%
%
% Heidt, A. M.; Hartung, A. & Bartelt, H.,
% Generation of Ultrashort and Coherent Supercontinuum Light Pulses in All-Normal Dispersion Fibers,
% The Supercontinuum Laser Source: The Ultimate White Light, Springer New York, 247-280 (2016)

close all; clearvars;

addpath('../../../../../GMMNLSE/GMMNLSE algorithm/','../../../../../GMMNLSE/user_helpers/');

%% Setup fiber parameters
sim.lambda0 = 1030e-9; % the center wavelength
sim.include_Raman = false; % no Raman
sim.gpu_yes = false; % don't use GPU
fiber.MFD = 3.0941*2;
sim.pulse_centering = false;

% Load default parameters like 
%
% load fiber.betas and fiber.SR based on your multimode folder above
% sim.f0 = 3e5/sim.lambda0; THz
% sim.save_period = 0; Save only the fields at the beginning and the end fiber
% sim.ellipticity = 0; Linear polarization
% sim.scalar = true; Use scalar propagation
% sim.adaptive_dz.threshold = 1e-5; the threshold of the adaptive-step method
% sim.gpu_yes = true; Use GPU
% sim.include_Raman = 1; Consider the Raman
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.pulse_centering = true; Always shift the pulse to the center of the time window
% sim.gpuDevice.Index = 1; Use the GPU device 1
% sim.progress_bar = true; Show the progress bar
% sim.progress_bar_name = ''; Empty name for the progress bar
% sim.cuda_dir_path = 'GMMNLSE/cuda'; Where the cuda files are
% ......
%
% Please check this function for details.
[fiber,sim] = load_default_GMMNLSE_propagate(fiber,sim); % load default parameters

num_save = 10;
fiber.L0 = 0.01;
%fiber.betas(4:end) = 0; % no higher-order betas
sim.save_period = fiber.L0/num_save;

%% Setup general parameters
Nt = 2^8; % the number of time points
time_window = 0.8; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

%% calculate fiber betas from silica refractive index
% This is important to correctly simulate the broadband situations.
% Taylor-series coefficients is only good in narrowband situations.

% Sellmeier coefficients
material = 'fused silica';
[a,b] = Sellmeier_coefficients(material);
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
n_silica = n_from_Sellmeier(lambda/1e3);

fiber.betas = n_silica*2*pi./(lambda*1e-9);

%% Initial condition
tfwhm = 0.03; % ps
total_energy = 10; % nJ
initial_pulse = build_MMgaussian(tfwhm, time_window, total_energy, 1, Nt);

%% Propagate
prop_output = GMMNLSE_propagate(fiber,initial_pulse,sim);

%% Plot
% Time
figure;
h = plot(t,abs(prop_output.fields(:,:,end)).^2);
xlim([-0.2,0.2]);
xlabel('t');
ylabel('Power');
title('Field');
set(h,'linewidth',2);
set(gca,'fontsize',14);

% Spectrum
figure;
h = plot(f-sim.f0,abs(fftshift(ifft(prop_output.fields(:,:,end)),1)).^2);
xlim([-100,100]);
xlabel('\nu-\nu_0');
ylabel('PSD');
title('Spectrum');
set(h,'linewidth',2);
set(gca,'fontsize',14);

% Comparison of time
figure;
[x,y] = meshgrid(t,prop_output.z);
pcolor(x,y,permute(abs(prop_output.fields(:,1,:)).^2,[3 1 2]));
shading interp; colormap(jet);
xlim([-0.2,0.2]);
xlabel('t');
ylabel('z');
title('Field during propagation');
set(gca,'fontsize',14);

% Comparison of spectra
figure;
[x,y] = meshgrid(f-sim.f0,prop_output.z(2:end));
pcolor(x,y,permute(abs(fftshift(ifft(prop_output.fields(:,1,2:end)),1)).^2,[3 1 2]));
xlim([-100,100]);
shading interp; colormap(jet);
xlabel('\nu-\nu_0');
ylabel('z/LD');
title('Spectrum during propagation');
set(gca,'fontsize',14);

%% Save the data
save('verification_GMMNLSE.mat');