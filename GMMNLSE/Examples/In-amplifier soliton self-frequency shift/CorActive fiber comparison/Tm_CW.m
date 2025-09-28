% This code computes amplification performance of a CW seed, and compares
% with the CorActive's datasheet.
%
% To do this simulation, I converge the CW power into a pulse in a time
% window. Nonlinear and dispersion are ignored. Pulse duration and time
% window size aren't important as long as the pulse is correctly covered by
% the time window.
%
% Vary the copump_power from 1 to 10 W for plotting with the plotter.

close all; clearvars;

addpath('../../../GMMNLSE algorithm/','../../../user_helpers/');

%% Gain info
% Please find details of all the parameters in "gain_info.m" if not specified here.
% Note that the use of single spatial mode is different from multi-spatial modes.
% "reuse_data" and "linear_oscillator_model" are activated and some related parameters are set.
gain_rate_eqn.gain_medium = 'Tm'; % specify the gain medium
gain_rate_eqn.base_medium = 'silica'; % specify the base medium
gain_rate_eqn.core_diameter = 10; % um
gain_rate_eqn.cladding_diameter = 128; % um
gain_rate_eqn.core_NA = 0.22;
gain_rate_eqn.absorption_wavelength_to_get_N_total = 790; % nm
gain_rate_eqn.absorption_to_get_N_total = 5; % dB/m
gain_rate_eqn.pump_wavelength = 790; % nm
gain_rate_eqn.copump_power = 10; % W;
gain_rate_eqn.counterpump_power = 0; % W; it's set below
gain_rate_eqn.reuse_data = false; % For a ring or linear cavity, the pulse will enter a steady state eventually.
                                  % If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
gain_rate_eqn.linear_oscillator = false; % For a linear oscillator, there are pulses from both directions simultaneously, which will deplete the gain;
                                         % therefore, the backward-propagating pulses need to be taken into account.
gain_rate_eqn.t_rep = 1/10e6; % Assume 10 MHz here; s; the time required to finish a roundtrip (the inverse repetition rate of the pulse)
                              % This gain model solves the gain of the fiber under the steady-state condition; therefore, the repetition rate must be high compared to the lifetime of the doped ions.
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 50; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-3; % the tolerance for the above iterations
gain_rate_eqn.verbose = true; % show the information(final pulse energy) during iterations of computing the gain

%% Setup fiber parameters
% load_default_GMMNLSE_propagate loads MFD=6.2um.

% General parameters
sim.lambda0 = 2000e-9; % m
sim.f0 = 2.99792458e-4/sim.lambda0; % THz
sim.dz = 2000e-6; % m
sim.save_period = 0.05;
sim.gpu_yes = false; % For only the fundamental mode, running with CPU is faster if the number of points is lower than 2^(~18).

% -------------------------------------------------------------------------
% -------------------------------- Arm (6um) ------------------------------
% -------------------------------------------------------------------------
% Gain fiber
sim_Gain = sim;
sim_Gain.gain_model = 2;
fiber_Gain.L0 = 3; % the length of the gain fiber
fiber_Gain.MFD = 1e10; % ignore nonlinearity

% Load default parameters like 
%
% loading fiber.betas and fiber.SR based on your multimode folder above
% sim.include_Raman = true; Consider Raman
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.gpu_yes = true; Use GPU (default to true)
% ......
%
% Please check this function for details.
[fiber_Gain,sim_Gain] = load_default_GMMNLSE_propagate(fiber_Gain,sim_Gain,'single_mode');

%% Setup general cavity parameters
max_rt = 100; % maximum number of roundtrips
Nt = 2^9; % the number of points
time_window = 10; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % 

%% Setup initial conditions
input_power = 1; % W
tfwhm = 0.1; % ps; not important
total_energy = input_power*gain_rate_eqn.t_rep*1e9; % nJ
prop_output = build_MMgaussian(tfwhm, time_window, total_energy,1,Nt);

%% Load gain parameters
gain_rate_eqn = gain_info( fiber_Gain,sim_Gain,gain_rate_eqn,ifftshift(lambda,1) );

%% Run the simulation
prop_output = GMMNLSE_propagate(fiber_Gain, prop_output, sim_Gain, gain_rate_eqn);

output_power = trapz(t,abs(prop_output.fields(:,:,end)).^2/1e3*1e-9/gain_rate_eqn.t_rep);

save(sprintf('Tm_CW_%uW.mat',gain_rate_eqn.copump_power),'output_power');