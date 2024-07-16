close all; clearvars;

addpath('../../../GMMNLSE algorithm','../../../user_helpers');

%% Setup fiber parameterss
% Please find details of all the parameters in "load_default_GMMNLSE_propagate.m".
% Only necessary parameters are set here; otherwise, default is used.
sim.lambda0 = 1060e-9; % m
sim.f0 = 2.99792458e-4/sim.lambda0;
%sim.gpu_yes = false;
sim.save_period = 0.1; % m

% -------------------------------------------------------------------------

% Gain fiber
sim.gain_model = 2;
sim.progress_bar_name = 'Gain';
fiber.L0 = 5.7; % m; fiber length
fiber.MFD = sqrt(410/pi)*2; % um; mode-field diameter
fiber.n2 = 3.07e-20; % m^2/W; nonlinear refractive index

% Load default parameters like 
%
% loading fiber.betas and fiber.SR based on your multimode folder above
% sim.Raman_model = 1; Use isotropic Raman model
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.gpu_yes = true; Use GPU (default to true)
% ......
%
% Please check this function for details.
[fiber,sim] = load_default_GMMNLSE_propagate(fiber,sim,'single_mode'); % load default parameters for "fiber" and "sim"

%% Gain info
% Please find details of all the parameters in "gain_info.m" if not specified here.
% Note that the use of single spatial mode is different from multi-spatial modes.
% Activating "reuse_data" or "linear_oscillator_model" requires setting other parameters.
% Check the example or "gain_info.m".
gain_rate_eqn.gain_medium = 'Yb'; % specify the gain medium
gain_rate_eqn.base_medium = 'silica'; % specify the base medium
gain_rate_eqn.core_diameter = 20; % um
gain_rate_eqn.cladding_diameter = 400; % um
gain_rate_eqn.core_NA = 0.06;
gain_rate_eqn.absorption_wavelength_to_get_N_total = 915; % nm
gain_rate_eqn.absorption_to_get_N_total = 0.5; % dB/m
gain_rate_eqn.pump_wavelength = 976; % nm
gain_rate_eqn.copump_power = 0; % W
gain_rate_eqn.counterpump_power = 14*0.9; % W
gain_rate_eqn.reuse_data = false; % For a ring or linear cavity, the pulse will enter a steady state eventually.
                                  % If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
gain_rate_eqn.linear_oscillator = false; % For a linear oscillator, there are pulses from both directions simultaneously, which will deplete the gain;
                                         % therefore, the backward-propagating pulses need to be taken into account.
gain_rate_eqn.t_rep = 1/212e6; % Assume 212 MHz here; s; the time required to finish a roundtrip (the inverse repetition rate of the pulse)
                               % This gain model solves the gain of the fiber under the steady-state condition; therefore, the repetition rate must be high compared to the lifetime of the doped ions.
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 50; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the above iterations
gain_rate_eqn.verbose = true; % show the information(final pulse energy) during iterations of computing the gain

%% Setup general parameters
Nt = 2^12; % the number of time points
time_window = 10; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

% In my Yb cross-section repository, there are two data:
% 1. Yb_Gen_VIII_Cross_Section (Nufern)
% 2. Liekki Yb_AV_20160530
% I realized that to fit the figures from Lindberg's paper, it's required
% to use Nufern's data, rather than Liekki's. They have slightly different
% shape which results in different spectral shapes.
%
% This line of code needs to be before calling gain_info() so that gain_info() can read the corresponding cross-sectional data.
gain_rate_eqn.cross_section_filename = 'Yb_Gen_VIII_Cross_Section (Nufern).txt';

% Precompute some parameters related to the gain to save the computational time
% Check "gain_info.m" for details.
gain_rate_eqn = gain_info( fiber,sim,gain_rate_eqn,ifftshift(lambda,1) );

%% calculate fiber betas from silica refractive index
% This is important to correctly simulate the broadband situations.
% Taylor-series coefficients is only good in narrowband situations.

% Sellmeier coefficients
material = 'fused silica';
[a,b] = Sellmeier_coefficients(material);
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
n_silica = n_from_Sellmeier(lambda/1e3);

fiber.betas = n_silica*2*pi./(lambda*1e-9);

%% Setup initial conditions
tfwhm = 0.173; % ps
total_energy = 375/212*0.8; % nJ

pulse_lambda0 = 1041e-9;
f_now = c/sim.lambda0*1e-12;
f_pulse = c/pulse_lambda0*1e-12;
freq_shift = f_pulse - f_now;
prop_output = build_MMgaussian(tfwhm, time_window, total_energy, 1, Nt, {'ifft',freq_shift});

%% Run the cavity simulation
prop_output = GMMNLSE_propagate(fiber, prop_output, sim, gain_rate_eqn);

%% Finish the simulation and save the data
% Energy of the output field
energy = sum(trapz(abs(prop_output.fields(:,:,end)).^2))*prop_output.dt/10^3; % energy in nJ

fprintf('counterpump_%4.1fW_output_E_%4.1fW\n',gain_rate_eqn.counterpump_power,energy/gain_rate_eqn.t_rep/1e9);

factor_correct_unit = time_window^2/1e3; % to make the spectrum of the correct unit "nJ/THz"
                                         % "/1e3" is to make pJ into nJ
spectrum = abs(fftshift(ifft(prop_output.fields(:,:,end)),1)).^2*factor_correct_unit; % in frequency domain
factor = c/1e3./lambda.^2; % change the spectrum from frequency domain into wavelength domain

figure;
plot(lambda,spectrum.*factor,'linewidth',2);
xlabel('Wavelength (nm)'); ylabel('Spectrum (nJ/nm)');
set(gca,'fontsize',16);
xlim([1020,1100]);
print(sprintf('counterpump_%4.1fW.jpg',gain_rate_eqn.counterpump_power),'-djpeg');

%% Save data
save(sprintf('counterpump_%4.1fW.mat',gain_rate_eqn.counterpump_power));