% This code runs a CPA with narrowband transformation to speed up
% simulations. Please compare with "CPA_cs1.m" that applies no speedup from
% the narrowband transformation.
%
% If this is the first time you see scaled Fourier transform and its
% narrowband transformation, please check the reference below and run
% another "SPM_of_chirped_pulse.m" example that is simpler and has
% comparisons.
%
% The pulse is initially stretched by a long passive fiber. Narrowband
% transformation cannot be applied to this stretching process, because it
% begins with a transform-limited pulse that narrowband transformation is
% not applicable to.
%
% This code uses adaptive-step RK4IP with newly-proposed narrowband 
% transformation based on scaled Fourier transform.
%
% Narrowband transformation is applied to speed up the simulation, which is
% controlled by sim.cs.cs, the "scaled factor." It is a method to simulate
% pulse propagation for pulses with fewer sampling points below the
% numerical Nyquist limit. It only applies to pulses with smooth profiles
% in both spectral and temporal domains, such as a highly-chirped pulse in
% a CPA system whose temporal profile is a smooth long pulse with a smooth
% dominantly-quadratic phase, along with its smooth spectral profile.
%
% For scaled Fourier transform and its narrowband transformation, please
% find details in
% Chen and Wise, "Field-based treatment of transient gain in short-pulse
% optical amplifiers," Optica 12(6), 879-889 (2025).

close all; clearvars;

addpath('../../GMMNLSE algorithm/','../../user_helpers/');

%% Setup fiber parameterss
% Please find details of all the parameters in "load_default_GMMNLSE_propagate.m".
% Only necessary parameters are set here; otherwise, defaults are used.
sim.lambda0 = 1040e-9;
sim.f0 = 2.99792458e-4/sim.lambda0;
sim.gpu_yes = false;

% -------------------------------------------------------------------------

% Gain fiber
sim_Gain = sim;
sim_Gain.gain_model = 2;
sim_Gain.progress_bar_name = 'Gain (6um)';
fiber_Gain.L0 = 1.5; % m; fiber length
fiber_Gain.MFD = 35; % um; mode-field diameter

% Passive fiber
fiber_passive.L0 = 0.1; % m; fiber length

% Load default parameters like 
%
% loading fiber.betas and fiber.SR based on your multimode folder above
% sim.include_Raman = true; Consider Raman
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.gpu_yes = true; Use GPU (default to true)
% ......
%
% Please check this function for details.
[fiber_Gain,sim_Gain] = load_default_GMMNLSE_propagate(fiber_Gain,sim_Gain,'single_mode'); % load default parameters for "fiber" and "sim"

sim_Gain.save_period = fiber_Gain.L0/10;

[fiber_passive,sim_passive] = load_default_GMMNLSE_propagate(fiber_passive,sim,'single_mode'); % load default parameters for "fiber" and "sim"
sim_passive.save_period = fiber_passive.L0/10;
sim_passive.gpu_yes = true;

%% Gain info
% Please find details of all the parameters in "gain_info.m" if not specified here.
% Note that the use of single spatial mode is different from multi-spatial modes.
% Activating "reuse_data" or "linear_oscillator_model" requires setting other parameters.
% Check the example or "gain_info.m".
gain_rate_eqn.gain_medium = 'Yb'; % specify the gain medium
gain_rate_eqn.base_medium = 'silica'; % specify the base medium
gain_rate_eqn.core_diameter = 35; % um
gain_rate_eqn.cladding_diameter = 255; % um
gain_rate_eqn.core_NA = 0.07;
gain_rate_eqn.absorption_wavelength_to_get_N_total = 915; % nm
gain_rate_eqn.absorption_to_get_N_total = 2.8; % dB/m
gain_rate_eqn.pump_wavelength = 976; % nm
gain_rate_eqn.copump_power = 5; % W
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.reuse_data = false; % For a ring or linear cavity, the pulse will enter a steady state eventually.
                                  % If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
gain_rate_eqn.linear_oscillator = false; % For a linear oscillator, there are pulses from both directions simultaneously, which will deplete the gain;
                                         % therefore, the backward-propagating pulses need to be taken into account.
gain_rate_eqn.t_rep = 1/5e6; % Assume 5 MHz here; s; the time required to finish a roundtrip (the inverse repetition rate of the pulse)
                             % This gain model solves the gain of the fiber under the steady-state condition; therefore, the repetition rate must be high compared to the lifetime of the doped ions.
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 50; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the above iterations
gain_rate_eqn.verbose = true; % show the information(final pulse energy) during iterations of computing the gain

%% Setup general parameters
Nt = 2^22; % the number of time points
time_window = 10000; % ps
dt = time_window/Nt;
f = sim_Gain.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

% Precompute some parameters related to the gain to save the computational time
% Check "gain_info.m" for details.
cs = 32; % this scaled factor for the scaled Fourier transform needs to be transfered into gain_info()
sim_Gain.cs.cs = cs;
gain_rate_eqn = gain_info( fiber_Gain,sim_Gain,gain_rate_eqn,ifftshift(lambda,1) );

%% Setup initial conditions
ffwhm = 14/3; % THz
total_energy = 0.1; % nJ
initial_pulse = build_MMspectralGaussian(ffwhm, time_window, total_energy, 1, Nt, {'ifft',0}, 1,0,3);

%% Run the simulation: pulse temporal stretching
prop_stretch = GMMNLSE_propagate(fiber_passive, initial_pulse, sim_passive);

%% Apply the narrowband transformation
scaledFT_func = narrowband_scaledFT();

transformed_At = scaledFT_func.convert(prop_stretch.fields(:,:,end),cs);
transformed_t = linspace(-floor(Nt/cs/2), floor((Nt/cs-1)/2), Nt/cs)*(dt*cs);

% See the transformed field, which should have a narrower spectrum but the
% same temporal profile.
figure;
plot(transformed_t,abs(transformed_At).^2,'linewidth',2,'Color','b');
xlabel('Time (ps)');
ylabel('Power (W)');
set(gca,'fontsize',20);

hold on;
plot(t,abs(prop_stretch.fields(:,:,end)).^2,'linewidth',2,'Color','r');
hold off;
legend('Scaled-FT','Original');
title('Before and after transformation','FontSize',12);

transformed_f = sim_Gain.f0+(-Nt/cs/2:Nt/cs/2-1)'/(Nt*dt); % THz

abs2_transformed_Af = abs(fftshift(ifft(transformed_At,[],1),1)).^2;
norm_factor = max(abs2_transformed_Af);

figure;
plot(transformed_f,abs2_transformed_Af/norm_factor,'linewidth',2,'Color','b');
xlabel('Frequency (THz)');
ylabel('PSD');
set(gca,'fontsize',20);

hold on;

plot(f,abs(fftshift(ifft(prop_stretch.fields(:,:,end),[],1),1)).^2/norm_factor,'linewidth',2,'Color','r');
hold off;
legend('Scaled-FT','Original');
title('Before and after transformation','FontSize',12);

% ===========================================
% See Sec. 3 in the supplement of our Optica paper, reference above, for
% details regarding how to apply the scaled factor to physical quantities
% to obtain with the correct (physical-useful) units.
scaled_transformed_f = sim.f0+(-Nt/cs/2:Nt/cs/2-1)'/(Nt*dt)*cs; % THz

figure;
plot(scaled_transformed_f,abs2_transformed_Af/norm_factor,'linewidth',2,'Color','b');
xlabel('Frequency (THz)');
ylabel('PSD');
set(gca,'fontsize',20);

hold on;
plot(f,abs(fftshift(ifft(initial_pulse.fields(:,:,end),[],1),1)).^2/norm_factor*cs,'linewidth',2,'Color','r');
hold off;
legend('Scaled-FT','Original');
title({'Before and after transformation','(Plot w.r.t. original f in x-axis + normalization)'},'FontSize',12);
% ===========================================

% -------------------------------------------------------------------------
% Try to convert the transformed field back to the original scale
At = scaledFT_func.recover(transformed_At,cs);

figure;
plot(transformed_t,abs(transformed_At).^2,'linewidth',2,'Color','b');
xlabel('Time (ps)');
ylabel('Power (W)');
set(gca,'fontsize',20);

hold on;
plot(t,abs(At).^2,'linewidth',2,'Color','r');
hold off;
legend('Scaled-FT','Recovered');
title({'After transformation and then recovered','(Recovered one should be the same as that without transformation)'},'FontSize',12);

figure;
plot(transformed_f,abs(fftshift(ifft(transformed_At,[],1),1)).^2/norm_factor,'linewidth',2,'Color','b');
xlabel('Frequency (THz)');
ylabel('PSD');
set(gca,'fontsize',20);

hold on;
plot(f,abs(fftshift(ifft(At,[],1),1)).^2/norm_factor,'linewidth',2,'Color','r');
hold off;
legend('Scaled-FT','Recovered');
title({'After transformation and then recovered','(Recovered one should be the same as that without transformation)'},'FontSize',12);

%% Run the simulation: pulse amplification
prop_output = GMMNLSE_propagate(fiber_Gain, prop_stretch, sim_Gain, gain_rate_eqn);

%% Finish the simulation and save the data
% Energy of the output field
energy = squeeze(sum(trapz(abs(prop_output.fields).^2,1),2)*prop_output.dt/10^3); % energy in nJ

%[~,~,~,~,~,cb] = calc_spectrogram(t,f,prop_output.fields(:,:,end),[-1,1]*300,[1000,1250],100,100);
%colormap(whitejet_lower(512)); set(cb,'Color','[0 0 0]');

[Strehl_ratio,dechirped_FWHM,transform_limited_FWHM,peak_power] = analyze_field( t,f,prop_output.fields(:,:,end),'Treacy-t',pi/6,1e-6,true,false );

% Compute the gain
% Note that the original (non-transformed) cross sections are used,
% gain_rate_eqn.cross_sections0
Ntmp = permute(cat(4,1-prop_output.population,prop_output.population),[1,2,3,5,6,7,8,4])*gain_rate_eqn.N_total;
gain = pi*(gain_rate_eqn.core_diameter/2)^2*fftshift(permute(gain_rate_eqn.overlap_factor.signal.*sum(gain_rate_eqn.plusminus.*gain_rate_eqn.cross_sections0.*Ntmp(:,:,:,:,:,:,:,gain_rate_eqn.N_idx),8),[5,3,1,2,4]),1)*1e6;
g2 = gain; %g2(g2<0) = 0;
figure;
h = plot(lambda,10*log10(exp(1))*g2); set(h,'linewidth',2); clear h;
set(gca,'fontsize',20);
xlabel('Wavelength (nm)');
ylabel('Gain (dB/m)');
xlim([1000,1200]);
title('Gain spectrum');

figure;
plot(prop_output.z,energy,'linewidth',2);
xlabel('Propagation distance (m)');
ylabel('Pulse energy (nJ)');
set(gca,'fontsize',20);

nonlinear_phase = accumulated_nonlinear_phase( fiber_Gain.L0,1/fiber_Gain.SR,sim_Gain.f0,prop_output.fields,prop_output.z );
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