% This code runs the multimode Yb-doped fiber amplifier with the rate 
% equation gain and compare it with SM Gaussian and new gain models.
%
% The initial fields contain 3 modes with only one dominant field to check
% its consistency with other models under the single-mode configuration.
%
% Running the rate-equation model with only a fundamental mode is different
% from running it with multimode.
% With a fundamental mode, the spatial profile is assumed to be Gaussian
% such that many parameters can be pre-calculated; while for multimode,
% there's no assumption and the evolution of the spatial profiles needs to
% be computed.
%
% In the results, you'll see small discrepancies between running with a
% single mode or running with multimode (but only fundamental mode is
% nonzero). This is due to the fact that I use different models for two
% cases. In single-mode gain-rate-eqn model, the mode properties are
% calculated from "Whitley''s Gaussian-mode approximation," while those in
% multimode simulations are calculated from mode-field data generated from
% our mode solver, calculated from scripts in "BuildFiber/."

clearvars; close all;

%% Add the folders of multimode files and others
addpath('../../GMMNLSE algorithm/','../../user_helpers/'); % add where many GMMNLSE-related functions like  "GMMNLSE_propagate" is
fiber.MM_folder = '../../Fibers/YB1200-10_125DC-PM_wavelength1030nm/';
fiber.betas_filename = 'betas.mat';
fiber.S_tensors_filename = 'S_tensors_3modes.mat';

%% Gain info
gain_rate_eqn.gain_medium = 'Yb'; % specify the gain medium
gain_rate_eqn.core_diameter = 10; % um
gain_rate_eqn.cladding_diameter = 125; % um
gain_rate_eqn.core_NA = 0.08; % in fact, this is only used in single-mode
gain_rate_eqn.absorption_wavelength_to_get_N_total = 920; % nm
gain_rate_eqn.absorption_to_get_N_total = 1.7; % dB/m
gain_rate_eqn.pump_wavelength = 976; % nm
gain_rate_eqn.copump_power = 1; % W
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.reuse_data = false; % For a ring or linear cavity, the pulse will enter a steady state eventually.
                                  % If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
gain_rate_eqn.linear_oscillator = false; % For a linear oscillator, there are pulses from both directions simultaneously, which will deplete the gain;
                                         % therefore, the backward-propagating pulses need to be taken into account.
gain_rate_eqn.t_rep = 1/15e6; % assume 15 MHz here; s; the time required to finish a roundtrip (the inverse repetition rate of the pulse)
                             % This gain model solves the gain of the fiber under the steady-state condition; therefore, the repetition rate must be high compared to the lifetime of the doped ions.
gain_rate_eqn.downsampling_factor = 1; % an integer; downsample the eigenmode profiles to run faster
gain_rate_eqn.export_N2 = true; % whether to export N2, the ion density in the upper state or not
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 10; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the iteration
gain_rate_eqn.verbose = true; % show the information(final pulse energy) during iterations of computing the gain

%% Field and simulation parameters
time_window = 30; % ps
N = 2^12; % the number of time points
dt = time_window/N;
t = (-N/2:N/2-1)'*dt; % ps

fiber.L0 = 2; % m; the length of the gain fiber
fiber.MFD = 10.528; % um; the mode-field diameter
fiber.t_rep = gain_rate_eqn.t_rep; % for calculating the saturation intensity
save_num = 200;
sim.save_period = fiber.L0/save_num;
%sim.progress_bar = false;

sim.lambda0 = 1030e-9; % center wavelength; in "m"

% Some of the parameters are actually found by running the rate-equation
% simulations first and are set into SM and new-gain models, such as MFD,
% dB_gain, and saturation_intensity for a better fit.
%
% Note here that the saturation intensity calculated from h*f/(sigma*tau),
% where f is the center frequency,
%       sigma is the sum of the emission and absorption cross sections,
%       tau is the lifetime of the higher energy level for population inversion,
% is 22, which is pretty close to 30 I set here for a better fit between
% SM/new-gain models and the rate-equation-gain models.

% Gaussian gain
sim_GaussianGain = sim;
sim_GaussianGain.gain_model = 1;
fiber_GaussianGain = fiber;
fiber_GaussianGain.dB_gain = 55.5; % the small-signal dB gain (found after running rate-equation-gain model)
fiber_GaussianGain.saturation_intensity = 30; % J/m^2 (found after running rate-equation-gain model)
[fiber_GaussianGain,sim_GaussianGain] = load_default_GMMNLSE_propagate(fiber_GaussianGain,sim_GaussianGain,'multimode');

% Rate-equation gain with single-mode scheme
sim_SMrategain = sim;
sim_SMrategain.gain_model = 2;
sim_SMrategain.midx = 1; % choose only the first mode
[fiber_SMrategain,sim_SMrategain] = load_default_GMMNLSE_propagate(fiber,sim_SMrategain,'multimode');

% Rate-equation gain with multimode scheme
sim_MMrategain = sim;
sim_MMrategain.gain_model = 2;
[fiber_MMrategain,sim_MMrategain] = load_default_GMMNLSE_propagate(fiber,sim_MMrategain,'multimode');

fiber = {fiber_GaussianGain fiber_SMrategain fiber_MMrategain};
sim = {sim_GaussianGain sim_SMrategain sim_MMrategain};

%% Initial pulse
total_energy = 0.01; % nJ
tfwhm = 1; % ps; pulse duration
SM_field = build_MMgaussian(tfwhm, time_window, total_energy, 1, N);

MM_field = SM_field;
MM_field.fields = [SM_field.fields(:,1) zeros(N,2)];

input_field = [MM_field SM_field MM_field];

%% Gain parameters
% We need some parameters of gain before computations.
f = ifftshift( (-N/2:N/2-1)'/N/dt + sim{1}.f0 ); % in the order of "omegas" in the "GMMNLSE_propagate.m"
c = 299792.458; % nm/ps;
lambda = c./f; % nm

% SM
SM_gain_rate_eqn = gain_info( fiber{2},sim{2},gain_rate_eqn,lambda );
% MM
MM_gain_rate_eqn = gain_info( fiber{3},sim{3},gain_rate_eqn,lambda );

gain_rate_eqn = {[],SM_gain_rate_eqn,MM_gain_rate_eqn};

%% Propagation
t_end = zeros(1,3);
model_name = {'Gaussian gain','SM_rate-eqn gain','MM_rate-eqn gain'};
output_field = cell(1,3);
for i = 1:3
    output_field{i} = GMMNLSE_propagate(fiber{i},input_field(i),sim{i},gain_rate_eqn{i});
    t_spent = datevec(output_field{i}.seconds/3600/24);
    fprintf('Running time for %s: %2u:%3.1f\n',model_name{i},t_spent(5),t_spent(6));
end

%% Plot results
energy_GaussianGain  = permute(sum(trapz(abs(output_field{1}.fields).^2),2)*dt/1e3,[3 2 1]);
energy_SMrategain = permute(sum(trapz(abs(output_field{2}.fields).^2),2)*dt/1e3,[3 2 1]);
energy_MMrategain = permute(sum(trapz(abs(output_field{3}.fields).^2),2)*dt/1e3,[3 2 1]);

% Energy
figure;
plot(output_field{1}.z,[energy_GaussianGain energy_SMrategain energy_MMrategain]);
legend('Gaussian gain','SMrategain','MMrategain');
xlabel('Propagation length (m)');
ylabel('Energy (nJ)');
title('Energy');

c = 299792458e-12; % m/ps
f = (-N/2:N/2-1)'/N/dt+c/sim{1}.lambda0;
lambda = c./f*1e9;

% -------------------------------------------------------------------------
% Gaussian gain
% -------------------------------------------------------------------------
% Field
figure;
subplot(2,1,1);
plot(t,abs(output_field{1}.fields(:,:,end)).^2);
legend('mode 1','mode 2','mode 3');
xlabel('Time (ps)');
ylabel('Power (W)');
title('The final output field of YDFA (Gaussian gain)');

% Spectrum
subplot(2,1,2);
plot(lambda,abs(ifftshift(ifft(output_field{1}.fields(:,:,end)),1)).^2);
legend('mode 1','mode 2','mode 3');
xlabel('Wavelength (nm)');
ylabel('PSD (a.u.)');
title('The final output spectrum of YDFA (Gaussian gain)');
xlim([1000 1080]);

% -------------------------------------------------------------------------
% SM rategain
% -------------------------------------------------------------------------
% Field
figure;
subplot(2,1,1);
plot(t,abs(output_field{2}.fields(:,:,end)).^2);
xlabel('Time (ps)');
ylabel('Power (W)');
title('The final output field of YDFA (SM rate-equation gain)');

% Spectrum
subplot(2,1,2);
plot(lambda,abs(ifftshift(ifft(output_field{2}.fields(:,:,end)),1)).^2);
xlabel('Wavelength (nm)');
ylabel('PSD (a.u.)');
title('The final output spectrum of YDFA (SM rate-equation gain)');
xlim([1000 1080]);

% -------------------------------------------------------------------------
% MM rategain
% -------------------------------------------------------------------------
% Field
figure;
subplot(2,1,1);
plot(t,abs(output_field{3}.fields(:,:,end)).^2);
legend('mode 1','mode 2','mode 3');
xlabel('Time (ps)');
ylabel('Power (W)');
title('The final output field of YDFA (MM rate-equation gain)');

% Spectrum
subplot(2,1,2);
plot(lambda,abs(ifftshift(ifft(output_field{3}.fields(:,:,end)),1)).^2);
legend('mode 1','mode 2','mode 3');
xlabel('Wavelength (nm)');
ylabel('PSD (a.u.)');
title('The final output spectrum of YDFA (MM rate-equation gain)');
xlim([1000 1080]);

%% Save results
save('MM_YDFA_check_single_mode.mat');