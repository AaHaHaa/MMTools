% This file serves the purpose of letting people know how to add random
% mode coupling to a GMMNLSE simulation.
%
% 6 modes of pulses are sent into a multimode fiber.
% The input pulse are arranged into several modes with amplitudes that
% allows only a fundamental mode at the output through random mode
% coupling through the multimode fiber. The random mode coupling overall
% acts as a transmission matrix of a linear system.
%
%   Input pulse --> [Multimode fiber: a transmission matrix] --> Output pulse

close all; clearvars;

%% Load folders
addpath('../../GMMNLSE algorithm/','../../user_helpers/');

% Add multimode fiber informaiton
fiber.MM_folder = '../../Fibers/step_wavelength1550nm_diameter15um/';
fiber.betas_filename = 'betas.mat';
fiber.S_tensors_filename = 'S_tensors_6modes.mat';

%% GPU
sim.gpu_yes = true; % true(1) or false(0); Use GPU or not

%% Propagation information
% fiber.betas and fiber.SR will be loaded with
% "load_default_GMMNLSE_propagate.m" function below.

% Propagation length
fiber.L0 = 1; % fiber length (m)
sim.save_period = fiber.L0/100; % save the fields every...

% Simulation details
%sim.dz = 1e-6; % the simulation z step
sim.lambda0 = 1550e-9; % the center wavelength
sim.midx = [1 2 3]; % included modes;
                    % Because I have six modes of all my betas and S_tensors,
                    % I use this command to pick only the 1st to the 3rd modes.
num_modes = length(sim.midx)*2; % 2 is because of polarizations

sim.scalar = false; % true or false; consider polarized modes

% Load default parameters like 
%
% loading fiber.betas and fiber.SR based on your multimode folder above
% sim.include_Raman = true; Consider Raman
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.gpu_yes = true; Use GPU (default to true)
% ......
%
% Please check this function for details.
[fiber,sim] = load_default_GMMNLSE_propagate(fiber,sim,'multimode');

%% Preset some parameters
N = 2^13;
tfwhm = 10; % ps
time_window = 100; % ps
dt = time_window/N;
t = (-N/2:N/2-1)'*dt; % ps

%% Random mode coupling
addpath('../../Random_Mode_Coupling/');

sim.rmc.model = true;

% Strong coupling for all modes
sim.rmc.varn = 5e-2; % variance of the refractive index
sim.rmc.stdQ_polarizedmode = 10; % the coupling strength between strongly coupled polarization modes
sim.dz = 100e-6; % m
save_points = int32(fiber.L0/sim.dz);
sim.rmc.matrices = create_rmc_matrices(fiber,sim,num_modes,save_points);

% Calculate the effective transmission matrix of the entire multimode fiber
transmission_matrix = calc_effective_rmc_matrix(fiber,sim,N,dt,...
                                             sim.rmc.matrices);
% Use this transmission matrix to determine the mode amplitudes for the
% input pulse to obtain a fundamental mode at the fiber output
mode_amplitude = transmission_matrix\[1;zeros(num_modes-1,1)];

%% Initial condition/pulse
total_energy = 1e-9; % nJ; weak pulse for linear propagation
input_field = build_MMgaussian(tfwhm, time_window, total_energy, num_modes, N, {'ifft',0}, mode_amplitude);

output_field = GMMNLSE_propagate(fiber, input_field, sim);

save('RMC.mat','output_field');