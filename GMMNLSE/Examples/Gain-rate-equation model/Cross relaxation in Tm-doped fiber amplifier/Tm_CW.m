% This code computes the output ASE powers, both backward and forward, from
% a Tm-doped fiber amplifier, pumped by a 793-nm diode laser.
% 10% reflection from fiber end facets is considered.
% Due to high doping concentration, cross relaxation due to upconversion
% takes effect, resulting in higher slope efficiency for the one-to-two
% process than the one-to-one process that exhibits a 41% Stokes limit.
%
% This code should demonstrate 55% slope efficiency.
%
% Run several of this script with various pump powers to see the slope
% efficiency with "plotter_efficiency.m".
%
% For reference of this simulation, check p.21 of
% Jollivet et al., "HIGH PERFORMANCE LARGE-MODE AREA DOUBLE-CLAD FIBERS FOR
% KW POWER SCALING OF FIBER LASERS FROM 1 TO 2 MICRONS," 7th International
% Workshop on Specialty Optical Fibers (2022)
% https://indico.cern.ch/event/1131431/contributions/5044488/attachments/2527541/4473695/1130%20Thu%20HallA%20Jollivet.pdf

close all; clearvars;

addpath('../../../GMMNLSE algorithm/','../../../user_helpers/');

%% Gain info
% Please find details of all the parameters in "gain_info.m" if not specified here.
% Note that the use of single spatial mode is different from multi-spatial modes.
% "reuse_data" and "linear_oscillator_model" are activated and some related parameters are set.
gain_rate_eqn.gain_medium = 'Tm'; % specify the gain medium
gain_rate_eqn.base_medium = 'silica'; % specify the base medium
gain_rate_eqn.core_diameter = 25; % um
gain_rate_eqn.cladding_diameter = 400; % um
gain_rate_eqn.core_NA = 0.09;
gain_rate_eqn.absorption_wavelength_to_get_N_total = 793; % nm
gain_rate_eqn.absorption_to_get_N_total = 5.2; % dB/m
gain_rate_eqn.pump_wavelength = 793; % nm
gain_rate_eqn.copump_power = 30; % W; it's set below
gain_rate_eqn.counterpump_power = 0; % W; it's set below
gain_rate_eqn.reuse_data = false; % For a ring or linear cavity, the pulse will enter a steady state eventually.
                                  % If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
gain_rate_eqn.linear_oscillator = false; % For a linear oscillator, there are pulses from both directions simultaneously, which will deplete the gain;
                                         % therefore, the backward-propagating pulses need to be taken into account.
gain_rate_eqn.t_rep = 1/10e6; % Assume 10 MHz here; s; the time required to finish a roundtrip (the inverse repetition rate of the pulse)
                              % This gain model solves the gain of the fiber under the steady-state condition; therefore, the repetition rate must be high compared to the lifetime of the doped ions.
gain_rate_eqn.ignore_ASE = false;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 50; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-3; % the tolerance for the above iterations
gain_rate_eqn.verbose = true; % show the information(final pulse energy) during iterations of computing the gain

%% Setup fiber parameters
% load_default_GMMNLSE_propagate loads MFD=6.2um.

% General parameters
sim.lambda0 = 1985e-9; % m
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
fiber_Gain.L0 = 2; % the length of the gain fiber

% Load default parameters like 
%
% loading fiber.betas and fiber.SR based on your multimode folder above
% sim.Raman_model = 1; Use isotropic Raman model
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.gpu_yes = true; Use GPU (default to true)
% ......
%
% Please check this function for details.
[fiber_Gain,sim_Gain] = load_default_GMMNLSE_propagate(fiber_Gain,sim_Gain,'single_mode');

fiber_cavity = [fiber_Gain fiber_Gain];

%% Setup general cavity parameters
max_rt = 100; % maximum number of roundtrips
Nt = 2^9; % the number of points
time_window = 10; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % 
OC = 0.9; % output coupler

%% calculate fiber betas from silica refractive index
% Sellmeier coefficients
material = 'fused silica';
[a,b] = Sellmeier_coefficients(material);
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
n_silica = n_from_Sellmeier(lambda/1e3);

fiber.betas = n_silica*2*pi./(lambda*1e-9);

%% Setup initial conditions
tfwhm = 0.1; % ps
total_energy = 0; % nJ
prop_output = build_MMgaussian(tfwhm, time_window, total_energy,1,Nt);
prop_output.Power.ASE.forward = zeros(Nt,1);
prop_output.Power.ASE.backward = zeros(Nt,1);

%% Load gain parameters
gain_rate_eqn = gain_info( fiber_Gain,sim_Gain,gain_rate_eqn,ifftshift(lambda,1) );

%% Run the simulation
func = analyze_sim;

% Initialize some parameters
rt_num = 0;
output_power = zeros(1,max_rt);
while rt_num < max_rt
    rt_num = rt_num +1;
    
    t_iteration_start = tic;
    cprintf('*[1 0.5 0.31]','Iteration %d', rt_num);
    
    % ASE
    prop_output.Power.ASE.backward = prop_output.Power.ASE.forward(:,:,end)*(1-OC);
    prop_output.Power.ASE.forward = prop_output.Power.ASE.backward(:,:,1)*(1-OC);
    
    prop_output = GMMNLSE_propagate(fiber_Gain, prop_output, sim_Gain, gain_rate_eqn);
    
    % ---------------------------------------------------------------------
    % Finish propagation
    % ---------------------------------------------------------------------
    output_power(rt_num) = (trapz(f,prop_output.Power.ASE.forward(:,:,end))+trapz(f,prop_output.Power.ASE.backward(:,:,1)))*OC; % W
    
    if rt_num ~= 1
        close(fig_gain);
    end
    fig_gain = func.analyze_gain(prop_output.z,[],prop_output.Power.pump,prop_output.population);

    % ---------------------------------------------------------------------
    % Break if converged
    if rt_num > 1
        if abs(output_power(rt_num)/output_power(rt_num-1) - 1) < 1e-3
            cprintf('blue','ASE powers have converged!\n');
            break;
        end
    end
end
output_power = output_power(rt_num);
absorbed_pump_power = prop_output.Power.pump.forward(1) - prop_output.Power.pump.forward(end);

%% Finish the simulation and save the data
close(fig_gain);

save(sprintf('Tm_CW_%uW.mat',gain_rate_eqn.copump_power),'output_power','absorbed_pump_power');