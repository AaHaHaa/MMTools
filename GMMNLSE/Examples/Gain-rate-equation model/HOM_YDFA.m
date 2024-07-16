% This code runs the higher-order-mode Yb-doped fiber amplifier with the 
% rate equation gain.
%
% This simulation considers only LP11a to check if the code works correctly
% with only one higher-order mode. Since there is only one spatial mode,
% the code will run with single-mode scheme, which employs RK4IP, rather
% than MPA in multimode scenarios.

clearvars; close all;

%% Add the folders of multimode files and others
addpath('../../GMMNLSE algorithm/','../../user_helpers/'); % add where many GMMNLSE-related functions like  "GMMNLSE_propagate" is
fiber.MM_folder = '../../Fibers/YB1200-10_125DC-PM_wavelength1030nm/';
fiber.betas_filename = 'betas.mat';
fiber.S_tensors_filename = 'S_tensors_3modes.mat';

%% Gain info
gain_rate_eqn.gain_medium = 'Yb'; % specify the gain medium
gain_rate_eqn.reuse_data = false; % For a ring or linear cavity, the pulse will enter a steady state eventually.
                                  % If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
gain_rate_eqn.linear_oscillator = false; % For a linear oscillator, there are pulses from both directions simultaneously, which will deplete the gain;
                                         % therefore, the backward-propagating pulses need to be taken into account.
gain_rate_eqn.core_diameter = 10; % um
gain_rate_eqn.cladding_diameter = 125; % um
gain_rate_eqn.core_NA = 0.08; % in fact, this is only used in single-mode
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
gain_rate_eqn.max_iterations = 10; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the iteration
gain_rate_eqn.verbose = true; % show the information(final pulse energy) during iterations of computing the gain

%% Field and simulation parameters
time_window = 50; % ps
Nt = 2^12; % the number of time points
dt = time_window/Nt;
t = (-Nt/2:Nt/2-1)'*dt; % ps

fiber.L0 = 1; % m; the length of the gain fiber
save_num = 50;
sim.save_period = fiber.L0/save_num;
sim.lambda0 = 1030e-9; % central wavelength; in "m"
sim.midx = 2; % consider only LP11a (single spatial mode)
sim.gpu_yes = false; % Because of the single-mode simulation, using CPU is faster
%sim.progress_bar = false;

sim.gain_model = 2;
[fiber,sim] = load_default_GMMNLSE_propagate(fiber,sim,'multimode');

%% Initial pulse
total_energy = 0.1; % nJ
tfwhm = 1; % ps
input_field = build_MMgaussian(tfwhm, time_window, total_energy, length(sim.midx), Nt);
input_field.Power.ASE.forward = zeros(Nt,length(sim.midx));
input_field.Power.ASE.backward = zeros(Nt,length(sim.midx));

%% Gain parameters
% Precompute some parameters related to the gain to save the computational time
% Check "gain_info.m" for details.
f = ifftshift( (-Nt/2:Nt/2-1)'/Nt/dt + sim.f0 ); % in the order of "omegas" in the "GMMNLSE_propagate.m"
c = 299792.458; % nm/ps;
lambda = c./f; % nm
gain_rate_eqn = gain_info( fiber,sim,gain_rate_eqn,lambda );
    
%% Propagation
prop_output = GMMNLSE_propagate(fiber,input_field,sim,gain_rate_eqn);

%% Plot results
energy = permute(sum(trapz(abs(prop_output.fields).^2),2)*dt/1e3,[3 2 1]);

% Energy
figure;
plot(prop_output.z,energy);
xlabel('Propagation length (m)');
ylabel('Energy (nJ)');

c = 299792458e-12; % m/ps
f = (-Nt/2:Nt/2-1)'/Nt/dt+c/sim.lambda0;
lambda = c./f*1e9;

c = 299792.458; % nm/ps
factor = c./lambda.^2; % change the spectrum from frequency domain into wavelength domain

% Field
figure;
subplot(2,1,1);
plot(t,abs(prop_output.fields(:,:,end)).^2);
xlabel('Time (ps)');
ylabel('Power (W)');

% Spectrum
subplot(2,1,2);
plot(lambda,abs(fftshift(ifft(prop_output.fields(:,:,end)),1)).^2.*factor);
xlabel('Wavelength (nm)');
ylabel('PSD (a.u.)');
xlim([1010 1050]);