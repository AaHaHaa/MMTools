% This code finds the ASE evolutions in a single-mode Yb-doped fiber 
% amplifier with the gain rate equation.

clearvars; close all;

addpath('../../GMMNLSE algorithm/','../../user_helpers/');

%% Gain info
% Please find details of all the parameters in "gain_info.m".
% Note that the use of single spatial mode is different from multi-spatial modes.
gain_rate_eqn.cross_section_filename = 'Liekki Yb_AV_20160530.txt';
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
gain_rate_eqn.tau = 840e-6; % lifetime of Yb in F_(5/2) state (Paschotta et al., "Lifetme quenching in Yb-doped fibers"); in "s"
gain_rate_eqn.export_N2 = false; % whether to export N2, the ion density in the upper state or not
gain_rate_eqn.ignore_ASE = false;
gain_rate_eqn.max_iterations = 50; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-3; % the tolerance for the iteration
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

% Rate-equation gain
sim.gain_model = 2;
[fiber,sim] = load_default_GMMNLSE_propagate(fiber,sim);

%% Initial pulse
total_energy = 0; % nJ
tfwhm = 1; % ps
input_field = build_MMgaussian(tfwhm, time_window, total_energy, 1, Nt);
input_field.Power.ASE.forward = zeros(Nt,1);
input_field.Power.ASE.backward = zeros(Nt,1);

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
func = analyze_sim;
func.analyze_ASE(fftshift(f,1),prop_output.Power.ASE,prop_output.z);