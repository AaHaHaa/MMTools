% This code solves the 4-um-core Nd-doped gain-managed fiber amplifier 
% (GMNA).

close all; clearvars;

addpath('../../../../GMMNLSE algorithm/','../../../../user_helpers/');

%% Setup fiber parameterss
% Please find details of all the parameters in "load_default_GMMNLSE_propagate.m".
% Only necessary parameters are set here; otherwise, defaults are used.
sim.lambda0 = 920e-9;
sim.f0 = 2.99792458e-4/sim.lambda0;
sim.gpu_yes = false;
sim.save_period = 0.6;

% -------------------------------------------------------------------------

% Gain fiber
sim.gain_model = 2;
sim.progress_bar_name = 'Gain (4um)';
fiber.L0 = 3; % m; fiber length
fiber.betas = [9.8810e6; 4.8872e3; 0.0315; 2.0457e-5; 1.2737e-9];
fiber.MFD = 4.7; % um; mode-field diameter

% Load default parameters like 
%
% loading fiber.betas and fiber.SR based on your multimode folder above
% sim.include_Raman = true; Consider Raman
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
gain_rate_eqn.gain_medium = 'Nd'; % specify the gain medium
gain_rate_eqn.base_medium = 'silica'; % specify the base medium
gain_rate_eqn.core_diameter = 4.7; % um
gain_rate_eqn.cladding_diameter = 77; % um
gain_rate_eqn.core_NA = 0.16;
gain_rate_eqn.absorption_wavelength_to_get_N_total = 800; % nm
gain_rate_eqn.absorption_to_get_N_total = 0.6; % dB/m
gain_rate_eqn.pump_wavelength = 808; % nm
gain_rate_eqn.copump_power = 2; % W
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.reuse_data = false; % For a ring or linear cavity, the pulse will enter a steady state eventually.
                                  % If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
gain_rate_eqn.linear_oscillator = false; % For a linear oscillator, there are pulses from both directions simultaneously, which will deplete the gain;
                                         % therefore , the backward-propagating pulses need to be taken into account.
gain_rate_eqn.t_rep = 1/100e6; % Assume 30 MHz here; s; the time required to finish a roundtrip (the inverse repetition rate of the pulse)
                              % This gain model solves the gain of the fiber under the steady-state condition; therefore, the repetition rate must be high compared to the lifetime of the doped ions.
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 50; % For counterpumping or considering ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the iteration
gain_rate_eqn.verbose = true; % show the information(final pulse energy) during iterations of computing the gain

%% Setup general parameters
Nt = 2^13; % the number of time points
time_window = 100; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

% Precompute some parameters related to the gain to save the computational time
% Check "gain_info.m" for details.
gain_rate_eqn = gain_info( fiber,sim,gain_rate_eqn,ifftshift(lambda,1) );

%% calculate fiber betas from silica refractive index
% This is important to correctly simulate the broadband situations.
% Taylor-series coefficients is only good in narrowband situations.

% Sellmeier coefficients
material = 'silica';
[a,b] = Sellmeier_coefficients(material);
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
n_silica = n_from_Sellmeier(lambda/1e3);

fiber.betas = n_silica*2*pi./(lambda*1e-9);

%% Setup initial conditions
tfwhm = 0.5; % ps
total_energy = 1; % nJ

pulse_lambda0 = 920e-9;
f_now = c/sim.lambda0*1e-12;
f_pulse = c/pulse_lambda0*1e-12;
freq_shift = f_pulse - f_now;
prop_output = build_MMgaussian(tfwhm, time_window, total_energy, 1, Nt, {'ifft',freq_shift});

%% Run the simulation
prop_output = GMMNLSE_propagate(fiber, prop_output, sim, gain_rate_eqn);

%% Finish the simulation and save the data
% Energy of the output field
energy = squeeze(sum(trapz(abs(prop_output.fields).^2,1),2)*prop_output.dt/10^3); % energy in nJ

[Strehl_ratio,dechirped_FWHM,transform_limited_FWHM,peak_power] = analyze_field( t,f,prop_output.fields(:,:,end),'Treacy-t',pi/6,1e-6 );

func = analyze_sim;
fig_gain = func.analyze_gain(prop_output.z,[],prop_output.Power.pump,prop_output.population);

% Assume only two levels participate
Ntmp = permute(cat(4,1-sum(prop_output.population,4),prop_output.population),[1,2,3,5,6,7,8,4])*gain_rate_eqn.N_total;
gain = pi*(gain_rate_eqn.core_diameter/2)^2*fftshift(permute(gain_rate_eqn.overlap_factor.signal.*sum(gain_rate_eqn.plusminus.*gain_rate_eqn.cross_sections.*Ntmp(:,:,:,:,:,:,:,gain_rate_eqn.N_idx),8),[5,3,1,2,4]),1)*1e6;
g2 = gain; %g2(g2<0) = 0;
figure;
h = plot(lambda,10*log10(exp(1))*g2); set(h,'linewidth',2);
hold on;
plot(lambda,zeros(Nt,1),'linewidth',4,'LineStyle','--','Color','k');
hold off;
set(gca,'fontsize',20);
xlabel('Wavelength (nm)');
ylabel('Gain (dB/m)');
xlim([880,1000]);
%title('Gain spectrum');
print(gcf,'gain_spectrum.pdf','-dpdf');
set(h,'linewidth',6);
set(gca,'fontsize',30);
ylim([-5,10]);
xlabel(''); ylabel('');
print(gcf,'gain_spectrum (inset).pdf','-dpdf');