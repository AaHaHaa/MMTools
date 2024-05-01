% This code demonstrates the pulse compression of a SPM-broadened pulse.
%
% This code uses adaptive-step RK4IP for the passive-fiber propagation.

close all; clearvars;

addpath('../../../GMMNLSE algorithm/','../../../user_helpers/');

%% Setup fiber parameters
sim.lambda0 = 1030e-9; % the center wavelength
sim.Raman_model = 0; % no Raman
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
fiber.L0 = 0.1;
fiber.betas = zeros(size(fiber.betas)); % no dispersion
sim.save_period = fiber.L0/num_save;

%% Setup general parameters
Nt = 2^10; % the number of time points
time_window = 10; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

%% Initial condition
tfwhm = 0.4; % ps
total_energy = 15; % nJ
initial_pulse = build_MMgaussian(tfwhm, time_window, total_energy, 1, Nt);

%% Propagate
prop_output = GMMNLSE_propagate(fiber,initial_pulse,sim);

%% Compress the pulse with a prism pair
alpha = pi/3;
%[dechirped_separation,optimal_FWHM,dechirped_field] = pulse_compressor( 'prism',[],sim.lambda0*1e9,t,prop_output.fields(:,:,end),alpha,'N-SF10',true,true );

[Strehl_ratio,dechirped_FWHM,transform_limited_FWHM] = analyze_field( t,f,prop_output.fields(:,:,end),'prism',alpha,'N-SF10' );

%% Save the data
save('pulse_compression.mat');