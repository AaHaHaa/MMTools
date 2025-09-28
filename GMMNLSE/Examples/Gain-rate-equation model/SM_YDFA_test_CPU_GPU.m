% This code runs the single-mode Yb-doped fiber amplifier with the gain 
% rate equation and compare it with the SM Gaussian gain model.
%
% The purpose of this code is to see the performance difference between CPU
% and GPU.

clearvars; close all;

addpath('../../GMMNLSE algorithm/','../../user_helpers/');

%% Gain info
gain_rate_eqn.gain_medium = 'Yb'; % specify the gain medium
gain_rate_eqn.base_medium = 'silica'; % specify the base medium
gain_rate_eqn.reuse_data = false; % For a ring or linear cavity, the pulse will enter a steady state eventually.
                                  % If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
gain_rate_eqn.linear_oscillator = false; % For a linear oscillator, there are pulses from both directions simultaneously, which will deplete the gain;
                                         % therefore, the backward-propagating pulses need to be taken into account.
gain_rate_eqn.core_diameter = 6; % um
gain_rate_eqn.cladding_diameter = 125; % um
gain_rate_eqn.core_NA = 0.12;
gain_rate_eqn.absorption_wavelength_to_get_N_total = 920; % nm
gain_rate_eqn.absorption_to_get_N_total = 0.55; % dB/m
gain_rate_eqn.pump_wavelength = 976; % nm
gain_rate_eqn.copump_power = 0.4; % W
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.t_rep = 1/15e6; % assume 15 MHz here; s; the time required to finish a roundtrip (the inverse repetition rate of the pulse)
                              % This gain model solves the gain of the fiber under the steady-state condition; therefore, the repetition rate must be high compared to the lifetime of the doped ions.
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 10; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the iteration
gain_rate_eqn.verbose = true; % show the information(final pulse energy) during iterations of computing the gain

%% Field and simulation parameters
time_window = 50; % ps
N = 2^13; % the number of time points
dt = time_window/N;
t = (-N/2:N/2-1)'*dt; % ps

fiber.L0 = 2; % m; the length of the gain fiber
save_num = 100; % the number of saved data
sim.save_period = fiber.L0/save_num;

sim.lambda0 = 1030e-9; % central wavelength; in "m"

% SM Gaussian gain
sim_GaussianGain = sim;
sim_GaussianGain.gain_model = 1;
fiber_GaussianGain = fiber;
fiber_GaussianGain.dB_gain = 24.5; % the small-signal dB/m gain
fiber_GaussianGain.saturation_intensity = 1000; % J/m^2
[fiber_GaussianGain,sim_GaussianGain] = load_default_GMMNLSE_propagate(fiber_GaussianGain,sim_GaussianGain);

% Rate-equation gain
sim_rategain = sim;
sim_rategain.gain_model = 2;
[fiber_rategain,sim_rategain] = load_default_GMMNLSE_propagate(fiber,sim_rategain);

fiber = {fiber_GaussianGain fiber_rategain};
sim = {sim_GaussianGain sim_rategain};

%% Initial pulse
total_energy = 0.01; % nJ
tfwhm = 1; % ps
input_field = build_MMgaussian(tfwhm, time_window, total_energy, 1, N);

%% Gain parameters
% We need some parameters of gain before computations.
f = ifftshift( (-N/2:N/2-1)'/N/dt + sim{1}.f0 ); % in the order of "Omega" in the "GMMNLSE_propagate.m"
c = 299792.458; % nm/ps;
lambda = c./f; % nm
gain_rate_eqn = gain_info( fiber{2},sim{2},gain_rate_eqn,lambda );

%% Propagation
t_end = zeros(1,2);
model_name = {'Gaussian gain','rate-eqn gain'};
output_field = cell(1,2);

% GPU
disp('GPU:');
profile on
for i = 1:2
    output_field{i} = GMMNLSE_propagate(fiber{i},input_field,sim{i},gain_rate_eqn);
    t_spent = datevec(output_field{i}.seconds/3600/24);
    fprintf('Running time for %s: %2u:%3.1f\n',model_name{i},t_spent(5),t_spent(6));
end
profile off
profsave(profile('info'),'profile_results_GPU');

% CPU
disp('CPU:');
profile on
for i = 1:2
    sim{i}.gpu_yes = false;
    
    output_field{i} = GMMNLSE_propagate(fiber{i},input_field,sim{i},gain_rate_eqn);
    t_spent = datevec(output_field{i}.seconds/3600/24);
    fprintf('Running time for %s: %2u:%3.1f\n',model_name{i},t_spent(5),t_spent(6));
end
profile off
profsave(profile('info'),'profile_results_CPU');

%% Plot results
energy_GaussianGain   = permute(sum(trapz(abs(output_field{1}.fields).^2),2)*dt/1e3,[3 2 1]);
energy_rategain = permute(sum(trapz(abs(output_field{2}.fields).^2),2)*dt/1e3,[3 2 1]);

% Energy
figure;
plot(output_field{1}.z,[energy_GaussianGain energy_rategain]);
legend('Gaussian gain','rategain');
xlabel('Propagation length (m)');
ylabel('Energy (nJ)');
title('Energy');

c = 299792458e-12; % m/ps
f = (-N/2:N/2-1)'/N/dt+c/sim{1}.lambda0;
lambda = c./f*1e9;

c = 299792.458; % nm/ps
factor = c./lambda.^2; % change the spectrum from frequency domain into wavelength domain

% -------------------------------------------------------------------------
% Gaussian gain
% -------------------------------------------------------------------------
% Field
figure;
subplot(2,1,1);
plot(t,abs(output_field{1}.fields(:,:,end)).^2);
xlabel('Time (ps)');
ylabel('Power (W)');
title('The final output field of YDFA (Gaussian gain)');

% Spectrum
subplot(2,1,2);
plot(lambda,abs(fftshift(ifft(output_field{1}.fields(:,:,end)),1)).^2.*factor);
xlabel('Wavelength (nm)');
ylabel('PSD (a.u.)');
title('The final output spectrum of YDFA (Gaussian gain)');
xlim([1010 1080]);

% -------------------------------------------------------------------------
% rategain
% -------------------------------------------------------------------------
% Field
figure;
subplot(2,1,1);
plot(t,abs(output_field{2}.fields(:,:,end)).^2);
xlabel('Time (ps)');
ylabel('Power (W)');
title('The final output field of YDFA (rate-equation gain)');

% Spectrum
subplot(2,1,2);
plot(lambda,abs(fftshift(ifft(output_field{2}.fields(:,:,end)),1)).^2.*factor);
xlabel('Wavelength (nm)');
ylabel('PSD (a.u.)');
title('The final output spectrum of YDFA (rate-equation gain)');
xlim([1010 1080]);

%% Save results
save('SM_YDFA.mat');