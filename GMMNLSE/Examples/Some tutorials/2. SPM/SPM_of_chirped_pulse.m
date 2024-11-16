% This code demonstrates the self-phase modulation of a chirped pulse.
%
% This code uses adaptive-step RK4IP for the passive-fiber propagation.

close all; clearvars;

addpath('../../../GMMNLSE algorithm/','../../../user_helpers/');

%% Setup fiber parameters
sim.lambda0 = 1030e-9; % the center wavelength
sim.include_Raman = false; % no Raman
sim.gpu_yes = false; % don't use GPU

% Load default parameters like 
%
% load fiber.betas and fiber.SR based on your multimode folder above
% sim.f0 = 3e5/sim.lambda0; THz
% sim.save_period = 0; Save only the fields at the beginning and the end fiber
% sim.ellipticity = 0; Linear polarization
% sim.scalar = true; Use scalar propagation
% sim.adaptive_dz.model = 1; Use adaptive-step method
% sim.adaptive_dz.threshold = 1e-5; the threshold of the adaptive-step method
% sim.gpu_yes = true; Use GPU
% sim.step_method = 'RK4IP'; Use RK4IP algorithm for pulse propagation
% sim.include_Raman = true; Consider Raman
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.rmc.model = 0; Don't include random mode coupling
% sim.pulse_centering = true; Always shift the pulse to the center of the time window
% sim.gpuDevice.Index = 1; Use the GPU device 1
% sim.progress_bar = true; Show the progress bar
% sim.progress_bar_name = ''; Empty name for the progress bar
% sim.cuda_dir_path = 'GMMNLSE/cuda'; Where the cuda files are
% ......
%
% Please check this function for details.
[fiber,sim] = load_default_GMMNLSE_propagate([],sim); % load default parameters

num_save = 100;
fiber.L0 = 1;
fiber.betas = zeros(size(fiber.betas)); % no dispersion
sim.save_period = fiber.L0/num_save;

% Or you can use the GUI version.
% Remember to modify according to the above values.
%[fiber,sim] = load_default_GMMNLSE_propagate_GUI(); % load default parameters with GUI

%% Setup general parameters
Nt = 2^10; % the number of time points
time_window = 10; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

%% Initial condition
tfwhm = 0.1; % ps
total_energy = 5; % nJ
initial_pulse = build_MMgaussian(tfwhm, time_window, total_energy, 1, Nt);

% Chirp the pulse to 1 ps
chirped_tfwhm = 1; % ps
func = calc_chirp;
omega = ifftshift(2*pi*f,1); % 2*pi*THz
[~,chirped_pulse] = func.General( chirped_tfwhm,omega,ifft(initial_pulse.fields),1 );
initial_pulse.fields = chirped_pulse;

%% Propagate
prop_output = GMMNLSE_propagate(fiber,initial_pulse,sim);

%% Plot
% Time
figure;
h = plot(t,abs(prop_output.fields(:,:,end)).^2);
xlim([-3,3]);
xlabel('t');
ylabel('Power');
title('Field');
set(h,'linewidth',2);
set(gca,'fontsize',14);

% Spectrum
figure;
h = plot(f-sim.f0,abs(fftshift(ifft(prop_output.fields(:,:,end)),1)).^2);
xlim([-20,20]);
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
xlim([-3,3]);
xlabel('t');
ylabel('z');
title('Field during propagation');
set(gca,'fontsize',14);

% Comparison of spectra
figure;
[x,y] = meshgrid(f-sim.f0,prop_output.z(2:end));
pcolor(x,y,permute(abs(fftshift(ifft(prop_output.fields(:,1,2:end)),1)).^2,[3 1 2]));
xlim([-20,20]);
shading interp; colormap(jet);
xlabel('\nu-\nu_0');
ylabel('z/LD');
title('Spectrum during propagation');
set(gca,'fontsize',14);

%% Visualize the nonlinear phase
Ef = fftshift(ifft(ifftshift(prop_output.fields(:,:,end),1)),1);
fitted_order = 3;
verbose = true;
[quardratic_phase,cubic_phase,fitted_param,quintic_phase] = characterize_spectral_phase( f,Ef,fitted_order,verbose );

%% Save the data
save('SPM_of_chirped_pulse.mat');