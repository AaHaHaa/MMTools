% This code runs the simulation with two solitons in LP11a and LP11b with
% three cases:
%   1. no random mode coupling
%   2. random mode coupling occurs only for polarization modes (two mode groups)
%   3. random mode coupling occurs for all the modes (there's only one mode group)
%
%
% Because two solitons have different frequencies, they'll collide with
% each other during propagations. This code shows what will happen after
% the collision with and without random mode coupling for all the modes or
% only the polarization modes.
%
% Reference:
% Mecozzi et al., Nonlinear propagation in multi-mode fibers in the strong
% coupling regime, Opt. Express (2012)
%
% *Because I used a mode solver to get all the fiber parameters instead of
%  using those from the paper, some parameters will be slightly different
%  from those of the paper. For example, this is why I chose 80 soliton 
%  periods for the propagation distance instead of 150 in the paper.

close all; clearvars;

addpath('../../GMMNLSE algorithm/','../../user_helpers/');

fiber.MM_folder = '../../Fibers/step_wavelength1550nm_diameter15um/';
fiber.betas_filename = 'betas.mat';
fiber.S_tensors_filename = 'S_tensors_6modes.mat';

sim.progress_bar_name = 'soliton collision';
sim.include_Raman = false; % no Raman term
sim.scalar = 0; % consider polarized fields
sim.lambda0 = 1550e-9; % the central wavelength
sim.pulse_centering = false;
sim.midx = [2 3]; % choose only the LP11a and LP11b modes
num_spatial_modes = length(sim.midx);
num_modes = num_spatial_modes*2; % consider polarization modes

tau = 30; % ps
tfwhm = tau*asech(1/sqrt(2))*2; % fwhm of the two solitons
T = 300; % separation between two solitons initially
delta_omega = 0.005; % THz*2*pi; the frequency difference between two solitons
                     % If it's too large, the time for coupling between two pulses will be too short. I've tried 0.7/tau.
delta_f = delta_omega/2/pi;

c = 299792458;
N = 2^10; % the number of time grid points
time_window = 1000;
dt = time_window/N;
t = (-N/2:N/2-1)'*dt; % ps

[fiber,sim] = load_default_GMMNLSE_propagate(fiber,sim,'multimode');
fiber.betas(:,2) = fiber.betas(:,1);

% Used to build solitons
Aeff = get_Aeff(fiber.SR);

% Propagate throught 100km
num_save = 100;
fiber.L0 = 3000e3; % m
sim.save_period = fiber.L0/num_save;

%% Input solitons
soliton = build_MMsoliton(tfwhm*ones(1,num_spatial_modes), fiber.betas(3,:), 1./Aeff*0.7106853, sim.lambda0, time_window, num_spatial_modes, N, {'ifft',[-delta_f/2 delta_f/2]}, ones(1,num_spatial_modes),[-T/2 T/2]);

% Characteristic lengths
[ LD,LNL ] = characteristic_lengths( abs(soliton.fields).^2,t,c./sim.lambda0*1e-12,fiber.betas(3,:),Aeff );
z0 = pi./2*LD; % soliton period(length)

soliton.fields = [soliton.fields(:,1) zeros(size(t)) ... % extend to include polarization modes
                  soliton.fields(:,2) zeros(size(t))];

%% no random mode coupling
sim_noRMC = sim;

output_field = GMMNLSE_propagate(fiber, soliton, sim_noRMC);
save('soliton_collision_noRMC.mat','output_field','t');

%% Random mode coupling between polarization modes only
addpath('../../Random_Mode_Coupling/');

sim_RMC = sim;

sim_RMC.rmc.model = true;

% No coupling between spatial modes
sim_RMC.rmc.varn = 0; % variance of the refractive index
sim_RMC.rmc.stdQ_polarizedmode = 1e-2; % the coupling strength between strongly coupled polarization modes
sim_RMC.dz = 1e2; % m
save_points = int32(fiber.L0/sim_RMC.dz);
sim_RMC.rmc.matrices = create_rmc_matrices(fiber,sim_RMC,num_modes,save_points);

output_field = GMMNLSE_propagate(fiber, soliton, sim_RMC);
save('soliton_collision_strong_RMC_1.mat','output_field','t');

%% Strong random mode coupling between all modes
addpath('../../Random_Mode_Coupling/');

sim_RMC = sim;

sim_RMC.rmc.model = true;

% Strong coupling for all modes
sim_RMC.rmc.varn = 1e-5; % variance of the refractive index
sim_RMC.rmc.stdQ_polarizedmode = 1e-2; % the coupling strength between strongly coupled polarization modes
sim_RMC.dz = 1e2; % m
save_points = int32(fiber.L0/sim_RMC.dz);
sim_RMC.rmc.matrices = create_rmc_matrices(fiber,sim_RMC,num_modes,save_points);

output_field = GMMNLSE_propagate(fiber, soliton, sim_RMC);
save('soliton_collision_strong_RMC_2.mat','output_field','t');