% This code runs the multimode Yb-doped fiber amplifier with the rate 
% equation gain, including ASE.

clearvars; close all;

%% Add the folders of multimode files and others
addpath('../../GMMNLSE algorithm//','../../user_helpers//'); % add where many GMMNLSE-related functions like  "GMMNLSE_propagate" is
fiber.MM_folder = '../../Fibers/YB1200-10_125DC-PM_wavelength1030nm/';
fiber.betas_filename = 'betas.mat';
fiber.S_tensors_filename = 'S_tensors_3modes.mat';

%% Gain info
gain_rate_eqn.gain_medium = 'Yb'; % specify the gain medium
gain_rate_eqn.base_medium = 'silica'; % specify the base medium
gain_rate_eqn.reuse_data = false; % For a ring or linear cavity, the pulse will enter a steady state eventually.
                                  % If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
gain_rate_eqn.linear_oscillator = false; % For a linear oscillator, there are pulses from both directions simultaneously, which will deplete the gain;
                                         % therefore, the backward-propagating pulses need to be taken into account.
gain_rate_eqn.core_diameter = 10; % um
gain_rate_eqn.cladding_diameter = 125; % um
gain_rate_eqn.core_NA = 0.3; % in fact, this is only used in single-mode
gain_rate_eqn.absorption_wavelength_to_get_N_total = 920; % nm
gain_rate_eqn.absorption_to_get_N_total = 1.7; % dB/m
gain_rate_eqn.pump_wavelength = 976; % nm
gain_rate_eqn.copump_power = 1; % W
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.downsampling_factor = 1; % an integer; downsample the eigenmode profiles to run faster
gain_rate_eqn.t_rep = 1/15e6; % assume 15 MHz here; s; the time required to finish a roundtrip (the inverse repetition rate of the pulse)
                              % This gain model solves the gain of the fiber under the steady-state condition; therefore, the repetition rate must be high compared to the lifetime of the doped ions.
gain_rate_eqn.ignore_ASE = false;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 20; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the iteration
gain_rate_eqn.verbose = true; % show the information(final pulse energy) during iterations of computing the gain

%% Field and simulation parameters
time_window = 50; % ps
N = 2^12; % the number of time points
dt = time_window/N;
t = (-N/2:N/2-1)'*dt; % ps

fiber.L0 = 1; % m; the length of the gain fiber
save_num = 50;
sim.save_period = fiber.L0/save_num;
sim.lambda0 = 1030e-9; % central wavelength; in "m"
%sim.progress_bar = false;
sim.gpuDevice.Index = 1;
% Rate-equation gain
sim.gain_model = 2;
sim.gpuDevice.Index = 1; % choose which GPU to use if you have multiple GPUs: 1,2,3...
[fiber,sim] = load_default_GMMNLSE_propagate(fiber,sim,'multimode');

fiber.betas(:,[2,3]) = (fiber.betas(:,[2,3])+fiber.betas(:,[3,2]))/2;

%% Initial pulse
total_energy = 0.1; % nJ
tfwhm = 1; % ps
input_field = build_MMgaussian(tfwhm, time_window, total_energy, length(sim.midx), N);
input_field.Power.ASE.forward = zeros(size(input_field.fields));
input_field.Power.ASE.backward = zeros(size(input_field.fields));

%% Gain parameters
% Precompute some parameters related to the gain to save the computational time
% Check "gain_info.m" for details.
f = ifftshift( (-N/2:N/2-1)'/N/dt + sim.f0 ); % in the order of "Omega" in the "GMMNLSE_propagate.m"
c = 299792.458; % nm/ps;
lambda = c./f; % nm
gain_rate_eqn = gain_info( fiber,sim,gain_rate_eqn,lambda );

%% Propagation
output_field = GMMNLSE_propagate(fiber,input_field,sim,gain_rate_eqn);
t_spent = datevec(output_field.seconds/3600/24);
fprintf('Running time: %2u:%3.1f\n',t_spent(5),t_spent(6));

%% Plot results
energy_rategain = permute(sum(trapz(abs(output_field.fields).^2),2)*dt/1e3,[3 2 1]);

% Energy
figure;
plot(output_field.z,energy_rategain);
xlabel('Propagation length (m)');
ylabel('Energy (nJ)');
title('Energy');

c = 299792458e-12; % m/ps
f = (-N/2:N/2-1)'/N/dt+c/sim.lambda0;
lambda = c./f*1e9;

c = 299792.458; % nm/ps
factor = c./lambda.^2; % change the spectrum from frequency domain into wavelength domain

% Field
figure;
subplot(2,1,1);
plot(t,abs(output_field.fields(:,:,end)).^2);
legend('mode 1','mode 2','mode 3');
xlabel('Time (ps)');
ylabel('Power (W)');
title('The final output field of YDFA (rate-equation gain)');

% Spectrum
subplot(2,1,2);
plot(lambda,abs(fftshift(ifft(output_field.fields(:,:,end)),1)).^2.*factor);
legend('mode 1','mode 2','mode 3');
xlabel('Wavelength (nm)');
ylabel('PSD (a.u.)');
title('The final output spectrum of YDFA (rate-equation gain)');
xlim([1010 1060]);

func = analyze_sim;
func.analyze_ASE(f,output_field.Power.ASE,output_field.z);

% =========================================================================
figure;
plot(output_field.z,permute(trapz(abs(output_field.fields).^2)*dt/1e3,[3 2 1]));
xlabel('Propagation length (m)');
ylabel('Energy (nJ)');
legend('LP_{01}','LP_{11a}','LP_{11b}');
title('rate-equation gain');

%% Save results
close all;
save('MM_YDFA_with_ASE.mat');