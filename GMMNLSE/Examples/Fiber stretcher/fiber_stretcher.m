% This code demonstrates the dechirping of a stretched pulse through a 
% fiber.
% We use prism and Treacy dechirpers.
% Here, we can see that only prism dechirper can dechirp the pulse well
% because prism dechirper has a different sign of TOD from the fiber such
% that it not only compensates GDD but also those from TOD. On the other
% hand, Treacy dechirper, although it can compensate GDD from the fiber,
% only adds more TOD to the pulse.
%
%  Prism dechirper: -GDD, -TOD
% Treacy dechirper: -GDD, +TOD
%
% This code uses adaptive-step RK4IP for the passive-fiber propagation.

close all; clearvars;

addpath('../../GMMNLSE algorithm/','../../user_helpers/');

%% Setup fiber parameters
sim.lambda0 = 1080e-9; % the center wavelength
sim.gpu_yes = false; % don't use GPU

% Load default parameters like 
%
% load fiber.betas and fiber.SR based on your multimode folder above
% sim.f0 = 3e5/sim.lambda0; THz
% sim.save_period = 0; Save only the fields at the beginning and the end fiber
% sim.ellipticity = 0; Linear polarization
% sim.scalar = true; Use scalar propagation
% sim.adaptive_deltaZ.threshold = 1e-5; the threshold of the adaptive-step method
% sim.gpu_yes = true; Use GPU
% sim.Raman_model = 1; Use the isotropic Raman model
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.pulse_centering = true; Always shift the pulse to the center of the time window
% sim.num_photon_noise_per_bin = 0; Don't include photon shot noise
% sim.gpuDevice.Index = 1; Use the GPU device 1
% sim.progress_bar = true; Show the progress bar
% sim.progress_bar_name = ''; Empty name for the progress bar
% sim.cuda_dir_path = 'GMMNLSE/cuda'; Where the cuda files are
% ......
%
% Please check this function for details.
[fiber,sim] = load_default_GMMNLSE_propagate([],sim); % load default parameters

num_save = 100;
fiber.L0 = 5;
sim.save_period = fiber.L0/num_save;

%% Setup general parameters
Nt = 2^12; % the number of time points
time_window = 50; % ps
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
tfwhm = 0.05; % ps
total_energy = 0.1; % nJ
initial_pulse = build_MMgaussian(tfwhm, time_window, total_energy, 1, Nt);

%% Propagate
prop_output = GMMNLSE_propagate(fiber,initial_pulse,sim);

%% Dechirped by a prism compressor
alpha = pi/3;
[dechirped_separation_prism,optimal_FWHM_prism,dechirped_field_prism,prism_height] = pulse_compressor( 'prism',[],sim.lambda0*1e9,t,prop_output.fields(:,:,end),alpha,'fused silica',false,true );

%% Dechirped by a Treacy compressor
incident_angle = 30; % deg
grating_spacing = 1e-3/1000; % m
[dechirped_separation_Treacy,optimal_FWHM_Treacy,dechirped_field_Treacy] = pulse_compressor( 'Treacy-t',incident_angle*pi/180,sim.lambda0*1e9,t,prop_output.fields(:,:,end),grating_spacing );

%% Plot
% Time
figure;
plot(t,abs([prop_output.fields(:,:,1),dechirped_field_prism,dechirped_field_Treacy]).^2,'linewidth',2);
xlim([-5,5]);
xlabel('Time (ps)');
ylabel('Power (W)');
set(gca,'fontsize',14);