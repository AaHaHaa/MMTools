clearvars; close all;

tfwhm = 10000; % ps

%% Add the folders of multimode files and others
addpath('../../../../MMTools/GMMNLSE/GMMNLSE algorithm/','../../../../MMTools/GMMNLSE/user_helpers/'); % add where many GMMNLSE-related functions like  "GMMNLSE_propagate" is
fiber.MM_folder = '../../Fibers/GRIN_168_400_wavelength1030nm/';
fiber.betas_filename = 'betas.mat';
fiber.S_tensors_filename = 'S_tensors_21modes.mat';

%% Gain info
gain_rate_eqn.gain_medium = 'Yb'; % specify the gain medium
gain_rate_eqn.base_medium = 'silica'; % specify the base medium
gain_rate_eqn.core_diameter = 84; % um
gain_rate_eqn.cladding_diameter = 400; % um
gain_rate_eqn.core_NA = 0.2; % in fact, this is only used in single-mode
gain_rate_eqn.absorption_wavelength_to_get_N_total = 975; % nm
gain_rate_eqn.absorption_to_get_N_total = 40; % dB/m
gain_rate_eqn.pump_wavelength = 976; % nm
gain_rate_eqn.copump_power = 50; % W
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.reuse_data = false; % For a ring or linear cavity, the pulse will enter a steady state eventually.
                                  % If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
gain_rate_eqn.linear_oscillator = false; % For a linear oscillator, there are pulses from both directions simultaneously, which will deplete the gain;
                                         % therefore, the backward-propagating pulses need to be taken into account.
gain_rate_eqn.downsampling_factor = 5; % an integer; downsample the eigenmode profiles to run faster
gain_rate_eqn.t_rep = 1/10e3; % assume 20 kHz here; s; the time required to finish a roundtrip (the inverse repetition rate of the pulse)
                              % This gain model solves the gain of the fiber under the steady-state condition; therefore, the repetition rate must be high compared to the lifetime of the doped ions.
gain_rate_eqn.tau = 840e-6; % lifetime of Yb in F_(5/2) state (Paschotta et al., "Lifetme quenching in Yb-doped fibers"); in "s"
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 10; % For counterpumping or considering ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the iteration
gain_rate_eqn.verbose = true; % show the information(final pulse energy) during iterations of computing the gain

%% Field and simulation parameters
time_window = tfwhm*5; % ps
Nt = 2^9; % the number of time points
dt = time_window/Nt;
t = (-Nt/2:Nt/2-1)'*dt; % ps

fiber.L0 = 10; % m; the length of the gain fiber
save_num = 100;
sim.save_period = fiber.L0/save_num;
sim.lambda0 = 1030e-9; % central wavelength; in "m"
%sim.progress_bar = false;
sim.gain_model = 2;
[fiber,sim] = load_default_GMMNLSE_propagate(fiber,sim,'multimode');

%% Initial pulse
peak_power = 1; % kW
total_energy = peak_power*tfwhm; % nJ
input_field = build_MMgaussian(tfwhm, time_window, total_energy, length(sim.midx), Nt, {'ifft',0}, [ones(1,6),ones(1,length(sim.midx)-6)*0.01]);
rand_phase = rand(1,length(sim.midx));
input_field.fields = input_field.fields.*exp(1i*2*pi*rand_phase);

%% Gain parameters
% Precompute some parameters related to the gain to save the computational time
% Check "gain_info.m" for details.
f = ifftshift( (-Nt/2:Nt/2-1)'/Nt/dt + sim.f0 ); % in the order of "omegas" in the "GMMNLSE_propagate.m"
c = 299792.458; % nm/ps;
lambda = c./f; % nm
gain_rate_eqn = gain_info( fiber,sim,gain_rate_eqn,lambda );

%% Propagation
output_field = GMMNLSE_propagate(fiber,input_field,sim,gain_rate_eqn);
%t_spent = datevec(output_field.seconds/3600/24);
%fprintf('Running time for %s: %2u:%3.1f\n','rate-eqn gain',t_spent(5),t_spent(6));

%% Plot results
energy_rategain = permute(sum(trapz(abs(output_field.fields).^2),2)*dt/1e3,[3 2 1]);

c = 299792458e-12; % m/ps
f = (-Nt/2:Nt/2-1)'/Nt/dt+c/sim.lambda0;
lambda = c./f*1e9;

c = 299792.458; % nm/ps
factor = c./lambda.^2; % change the spectrum from frequency domain into wavelength domain

% Field
figure;
subplot(2,1,1);
plot(t,abs(output_field.fields(:,:,end)).^2);
xlabel('Time (ps)');
ylabel('Power (W)');
title('The final output field of YDFA (rate-equation gain)');

% Spectrum
subplot(2,1,2);
plot(lambda,abs(fftshift(ifft(output_field.fields(:,:,end)),1)).^2.*factor);
xlabel('Wavelength (nm)');
ylabel('PSD (a.u.)');
title('The final output spectrum of YDFA (rate-equation gain)');
%xlim([900 1200]);

% =========================================================================
each_energy = permute(trapz(abs(output_field.fields).^2)*dt/1e3,[3,2,1]);
figure;
subplot(2,1,1);
plot(output_field.z,each_energy./sum(each_energy,2),'linewidth',2);
ylabel('Energy ratio');
subplot(2,1,2);
plot(output_field.z,squeeze(max(sum(abs(output_field.fields).^2,2),[],1)/1e3),'linewidth',2);
ylabel('Energy (nJ)');
xlabel('Propagation length (m)');

close all;
save(sprintf('MM_YDFA%u_%ups.mat',round(rand(1)*1000),tfwhm));