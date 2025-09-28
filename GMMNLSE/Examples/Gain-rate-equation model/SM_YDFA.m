% This code runs the single-mode Yb-doped fiber amplifier with the gain 
% rate equation and compare it with the Gaussian gain model.

clearvars; close all;

addpath('../../GMMNLSE algorithm/','../../user_helpers/');

%% Gain info
% Please find details of all the parameters in "gain_info.m".
% Note that the use of single spatial mode is different from multi-spatial modes.
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
gain_rate_eqn.copump_power = 1; % W
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.t_rep = 1/15e6; % assume 15 MHz here; s; the time required to finish a roundtrip (the inverse repetition rate of the pulse)
                              % This gain model solves the gain of the fiber under the steady-state condition; therefore, the repetition rate must be high compared to the lifetime of the doped ions.
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 10; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the iteration
gain_rate_eqn.verbose = true; % show the information(final pulse energy) during iterations

%% Field and simulation parameters
time_window = 50; % ps
Nt = 2^13; % the number of time points
dt = time_window/Nt;
t = (-Nt/2:Nt/2-1)'*dt; % ps

fiber.L0 = 2; % m; the length of fiber length
fiber.MFD = 6; % um; mode-field diameter of the fiber
fiber.t_rep = gain_rate_eqn.t_rep; % for calculating saturation intensity for Gaussian gain model
save_num = 100; % the number of saved data
sim.save_period = fiber.L0/save_num;
sim.gpu_yes = false;

sim.lambda0 = 1030e-9; % central wavelength; in "m"

% Gaussian gain
% 100 dB/m is set intentionally to match the initial growth rate of GaussianGain.
% Typically it should be 30 dB/m for a good amplifier. Instead, the saturation energy reflects the varying pump power.
sim_GaussianGain = sim;
sim_GaussianGain.gain_model = 1;
fiber_GaussianGain = fiber;
fiber_GaussianGain.dB_gain = 50; % the small-signal dB gain
fiber_GaussianGain.saturation_intensity = 30; % J/m^2
fiber_GaussianGain.saturation_energy = fiber_GaussianGain.saturation_intensity*(pi*(fiber.MFD*1e-6/2)^2)*1e9; % ~0.85 nJ here
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
input_field = build_MMgaussian(tfwhm, time_window, total_energy, 1, Nt);

%% Gain parameters
% Precompute some parameters related to the gain to save the computational time
% Check "gain_info.m" for details.
f = ifftshift( (-Nt/2:Nt/2-1)'/Nt/dt + sim{1}.f0 ); % in the order of "Omega" in the "GMMNLSE_propagate.m"
c = 299792.458; % nm/ps;
lambda = c./f; % nm
gain_rate_eqn = gain_info( fiber{2},sim{2},gain_rate_eqn,lambda );

%% Propagation
if sim_GaussianGain.gpu_yes
    reset(gpuDevice);
end
t_end = zeros(1,2);
model_name = {'SM gain','rate-eqn gain'};
output_field = cell(1,2);
for i = 1:2
    t_start = tic;
    output_field{i} = GMMNLSE_propagate(fiber{i},input_field,sim{i},gain_rate_eqn);
    t_end(i) = toc(t_start);
    t_spent = datevec(t_end(i)/3600/24);
    fprintf('Running time for %s: %2u:%3.1f\n',model_name{i},t_spent(5),t_spent(6));
end

%% Plot results
energy_GaussianGain   = permute(sum(trapz(abs(output_field{1}.fields).^2),2)*dt/1e3,[3 2 1]);
energy_rategain = permute(sum(trapz(abs(output_field{2}.fields).^2),2)*dt/1e3,[3 2 1]);

% Energy
distance = (0:save_num)*sim{1}.save_period;
figure;
plot(distance,[energy_GaussianGain energy_rategain]);
legend('GaussianGain','rategain');
xlabel('Propagation length (m)');
ylabel('Energy (nJ)');
title('Energy');

c = 299792458e-12; % m/ps
f = (-Nt/2:Nt/2-1)'/Nt/dt+c/sim{1}.lambda0;
lambda = c./f*1e9;

c = 299792.458; % nm/ps
factor = c./lambda.^2; % change the spectrum from frequency domain into wavelength domain

% -------------------------------------------------------------------------
% GaussianGain
% -------------------------------------------------------------------------
% Field
figure;
subplot(2,1,1);
plot(t,abs(output_field{1}.fields(:,:,end)).^2);
xlabel('Time (ps)');
ylabel('Power (W)');
title('The final output field of YDFA (SM gain)');

% Spectrum
subplot(2,1,2);
plot(lambda,abs(fftshift(ifft(output_field{1}.fields(:,:,end)),1)).^2.*factor);
xlabel('Wavelength (nm)');
ylabel('PSD (a.u.)');
title('The final output spectrum of YDFA (SM gain)');
xlim([1000 1100]);

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
xlim([1000 1100]);

%% Save results
%save('SM_YDFA.mat');