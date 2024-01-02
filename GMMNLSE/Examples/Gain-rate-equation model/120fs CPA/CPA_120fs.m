% This code demonstrates the nonlinear CPA.
%
% For details, please refer to
% Kuznetsova et al., "Chirped-pulse amplification near the gain-narrowing
% limit of Yb-doped fiber using a reflection grism compressor," Appl. Phys.
% B 88, 515-518 (2007).
%
% Other recommended papers of nonlinear CPA, please read
% 1. Zhou et al., "Compensation of nonlinear phase shifts with third-order
%    dispersion in short-pulse fiber amplifiers," Opt. Express 13, 4869-
%    4877 (2005).
% 2. Kuznetsova and Wise, "Scaling of femtosecond Yb-doped fiber amplifiers
%    to tens of microjoule pulse energy via nonlinear chirped pulse 
%    amplification," Opt. Lett. 32, 2671-2673 (2007).

close all; clearvars;

addpath('../../../GMMNLSE algorithm/','../../../user_helpers/');

%% Setup fiber parameterss
% Please find details of all the parameters in "load_default_GMMNLSE_propagate.m".
% Only necessary parameters are set here; otherwise, defaults are used.
sim.lambda0 = 1040e-9;
sim.f0 = 2.99792458e-4/sim.lambda0;

% -------------------------------------------------------------------------

% Gain fiber
sim_Gain = sim;
sim_Gain.gain_model = 2;
sim_Gain.progress_bar_name = 'Gain (6um)';
fiber_Gain.L0 = 1.5;
fiber_Gain.MFD = 35;

% passive fiber
fiber_passive.L0 = 10;

% Load default parameters like 
%
% loading fiber.betas and fiber.SR based on your multimode folder above
% sim.Raman_model = 1; Use isotropic Raman model
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.gpu_yes = true; Use GPU (default to true)
% ......
%
% Please check this function for details.
[fiber_Gain,sim_Gain] = load_default_GMMNLSE_propagate(fiber_Gain,sim_Gain,'single_mode'); % load default parameters for "fiber" and "sim"

sim_Gain.save_period = fiber_Gain.L0/10;

[fiber_passive,sim_passive] = load_default_GMMNLSE_propagate(fiber_passive,sim,'single_mode'); % load default parameters for "fiber" and "sim"
sim_passive.save_period = fiber_passive.L0/10;

%% Gain info
% Please find details of all the parameters in "gain_info.m" if not specified here.
% Note that the use of single spatial mode is different from multi-spatial modes.
% Activating "reuse_data" or "linear_oscillator_model" requires setting other parameters.
% Check the example or "gain_info.m".
gain_rate_eqn.cross_section_filename = 'Liekki Yb_AV_20160530.txt';
gain_rate_eqn.core_diameter = 35; % um
gain_rate_eqn.cladding_diameter = 255; % um
gain_rate_eqn.core_NA = 0.07;
gain_rate_eqn.absorption_wavelength_to_get_N_total = 915; % nm
gain_rate_eqn.absorption_to_get_N_total = 2.8; % dB/m
gain_rate_eqn.pump_wavelength = 976; % nm
gain_rate_eqn.copump_power = 0; % W
gain_rate_eqn.counterpump_power = 4.5; % W
gain_rate_eqn.reuse_data = false; % For a ring or linear cavity, the pulse will enter a steady state eventually.
                                  % If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
gain_rate_eqn.linear_oscillator = false; % For a linear oscillator, there are pulses from both directions simultaneously, which will deplete the gain;
                                         % therefore, the backward-propagating pulses need to be taken into account.
gain_rate_eqn.t_rep = 1/5e6; % Assume 5 MHz here; s; the time required to finish a roundtrip (the inverse repetition rate of the pulse)
                             % This gain model solves the gain of the fiber under the steady-state condition; therefore, the repetition rate must be high compared to the lifetime of the doped ions.
gain_rate_eqn.tau = 840e-6; % lifetime of Yb in F_(5/2) state (Paschotta et al., "Lifetme quenching in Yb-doped fibers"); in "s"
gain_rate_eqn.export_N2 = true; % whether to export N2, the ion density in the upper state or not
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 50; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the above iterations
gain_rate_eqn.verbose = true; % show the information(final pulse energy) during iterations of computing the gain

%% Setup general parameters
Nt = 2^16; % the number of time points
time_window = 200; % ps
dt = time_window/Nt;
f = sim_Gain.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

% Precompute some parameters related to the gain to save the computational time
% Check "gain_info.m" for details.
gain_rate_eqn = gain_info( fiber_Gain,sim_Gain,gain_rate_eqn,ifftshift(lambda,1) );

%% Setup initial conditions
ffwhm = 14/3; % THz
total_energy = 0.2*0.5; % nJ
initial_pulse = build_MMspectralGaussian(ffwhm, time_window, total_energy, 1, Nt, {'ifft',0}, 1,0,3);

%% Run the simulations
prop_output0 = GMMNLSE_propagate(fiber_passive, initial_pulse, sim_passive); 
prop_output = GMMNLSE_propagate(fiber_Gain, prop_output0, sim_Gain, gain_rate_eqn); 

%% Finish the simulation and save the data
N2 = prop_output.N2;
pump = prop_output.Power.pump.forward;
output_field = prop_output.fields(:,:,end);

% -----------------------------------------------------------------
% Energy of the output field
energy = squeeze(sum(trapz(abs(prop_output.fields).^2,1),2)*prop_output.dt/10^3); % energy in nJ

%[~,~,~,~,~,cb] = calc_spectrogram(t,f,prop_output.fields(:,:,end),[-1,1]*300,[1000,1250],100,100);
%colormap(whitejet_lower(512)); set(cb,'Color','[0 0 0]');

[Strehl_ratio,dechirped_FWHM,transform_limited_FWHM,peak_power] = analyze_field( t,f,prop_output.fields(:,:,end),'Treacy-t',pi/6,1e-6,true,false );

gain = pi*(gain_rate_eqn.core_diameter/2)^2*fftshift(permute(gain_rate_eqn.overlap_factor.signal*((gain_rate_eqn.cross_sections.emission+gain_rate_eqn.cross_sections.absorption).*prop_output.N2*gain_rate_eqn.N_total-gain_rate_eqn.cross_sections.absorption*gain_rate_eqn.N_total),[5,3,1,2,4]),1)*1e6;
g2 = gain; %g2(g2<0) = 0;
figure;
h = plot(lambda,10*log10(exp(1))*g2); set(h,'linewidth',2); clear h;
set(gca,'fontsize',20);
xlabel('Wavelength (nm)');
ylabel('Gain (1/m)');
xlim([1000,1200]);
title('Gain spectrum');

figure;
plot(prop_output.z,energy,'linewidth',2);
xlabel('Propagation distance (m)');
ylabel('Pulse energy (nJ)');
set(gca,'fontsize',20);

nonlinear_phase = accumulated_nonlinear_phase( fiber_Gain.L0,1/fiber_Gain.SR,sim_Gain.f0,prop_output.fields,prop_output.z,dt );
fprintf('N_NL = %6.4f*pi\n',nonlinear_phase/pi);

figure;
factor = c/1e3./lambda.^2; % change the spectrum from frequency domain into wavelength domain
spectrum = squeeze(abs(fftshift(ifft(prop_output.fields),1)).^2).*factor;
spectrum = spectrum./max(spectrum);
plot(lambda,spectrum,'linewidth',2);
xlabel('Wavelength (nm)');
ylabel('PDS (norm.)');
set(gca,'fontsize',20);
xlim([980,1100]);

filtered_lambda = 1040;
filtered_pulse = edgepass_filter('lowpass', prop_output, sim.f0, filtered_lambda, 0.3,1,true);
analyze_field( t,f,filtered_pulse.fields(:,:,end),'Treacy-t',pi/6,1e-6,true,false );

% Save the final output field
%save('CPA_120fs.mat','-v7.3');