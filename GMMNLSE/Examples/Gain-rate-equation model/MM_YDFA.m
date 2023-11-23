% This code runs the multimode Yb-doped fiber amplifier with the rate 
% equation gain and compare it with Gaussian-gain models.

clearvars; close all;

%% Add the folders of multimode files and others
addpath('../../GMMNLSE algorithm/','../../user_helpers/'); % add where many GMMNLSE-related functions like  "GMMNLSE_propagate" is
fiber.MM_folder = '../../Fibers/YB1200-10_125DC-PM_wavelength1030nm/';
fiber.betas_filename = 'betas.mat';
fiber.S_tensors_filename = 'S_tensors_3modes.mat';

%% Gain info
gain_rate_eqn.cross_section_filename = 'Liekki Yb_AV_20160530.txt';
gain_rate_eqn.reuse_data = false; % For a ring or linear cavity, the pulse will enter a steady state eventually.
                                  % If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
gain_rate_eqn.linear_oscillator = false; % For a linear oscillator, there are pulses from both directions simultaneously, which will deplete the gain;
                                         % therefore, the backward-propagating pulses need to be taken into account.
gain_rate_eqn.core_diameter = 10; % um
gain_rate_eqn.cladding_diameter = 125; % um
gain_rate_eqn.core_NA = 0.08; % in fact, this is only used in single-mode
gain_rate_eqn.absorption_wavelength_to_get_N_total = 920; % nm
gain_rate_eqn.absorption_to_get_N_total = 1.7; % dB/m
gain_rate_eqn.pump_wavelength = 976; % nm
gain_rate_eqn.copump_power = 1; % W
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.downsampling_factor = 1; % an integer; downsample the eigenmode profiles to run faster
gain_rate_eqn.t_rep = 1/15e6; % assume 15MHz here; s; the time required to finish a roundtrip (the inverse repetition rate of the pulse)
                             % This gain model solves the gain of the fiber under the steady-state condition; therefore, the repetition rate must be high compared to the lifetime of the doped ions.
gain_rate_eqn.tau = 840e-6; % lifetime of Yb in F_(5/2) state (Paschotta et al., "Lifetme quenching in Yb-doped fibers"); in "s"
gain_rate_eqn.export_N2 = true; % whether to export N2, the ion density in the upper state or not
gain_rate_eqn.ignore_ASE = true;%false;
gain_rate_eqn.max_iterations = 10; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the iteration
gain_rate_eqn.verbose = true; % show the information(final pulse energy) during iterations of computing the gain

%% Field and simulation parameters
time_window = 50; % ps
Nt = 2^12; % the number of time points
dt = time_window/Nt;
t = (-Nt/2:Nt/2-1)'*dt; % ps

fiber.L0 = 1; % m; the length of the gain fiber
save_num = 50;
sim.save_period = fiber.L0/save_num;
sim.lambda0 = 1030e-9; % central wavelength; in "m"
%sim.progress_bar = false;

% ----------------------------- Gaussian gain -----------------------------
% dB_gain and saturation_intensity are adjusted to have a good bit with the rate-equation-gain model
sim_GaussianGain = sim;
sim_GaussianGain.gain_model = 1;
fiber_GaussianGain = fiber;
fiber_GaussianGain.gain_doped_diameter = 10; % um
fiber_GaussianGain.dB_gain = 68; % the small-signal dB gain
fiber_GaussianGain.saturation_intensity = 26; % J/m^2
[fiber_GaussianGain,sim_GaussianGain] = load_default_GMMNLSE_propagate(fiber_GaussianGain,sim_GaussianGain,'multimode');
% -------------------------- Rate-equation gain ---------------------------
sim_rategain = sim;
sim_rategain.gain_model = 2;
[fiber_rategain,sim_rategain] = load_default_GMMNLSE_propagate(fiber,sim_rategain,'multimode');

fiber = {fiber_GaussianGain,fiber_rategain};
sim = {sim_GaussianGain,sim_rategain};

%% Initial pulse
total_energy = 0.1; % nJ
tfwhm = 1; % ps
input_field = build_MMgaussian(tfwhm, time_window, total_energy, length(sim{1}.midx), Nt);
input_field.Power.ASE.forward = zeros(Nt,length(sim{1}.midx));
input_field.Power.ASE.backward = zeros(Nt,length(sim{1}.midx));

%% Gain parameters
% Precompute some parameters related to the gain to save the computational time
% Check "gain_info.m" for details.
f = ifftshift( (-Nt/2:Nt/2-1)'/Nt/dt + sim{1}.f0 ); % in the order of "omegas" in the "GMMNLSE_propagate.m"
c = 299792.458; % nm/ps;
lambda = c./f; % nm
gain_rate_eqn = gain_info( fiber{2},sim{2},gain_rate_eqn,lambda );
    
%% Propagation
t_end = zeros(1,length(fiber));
model_name = {'Gaussian gain','rate-eqn gain'};
output_field = cell(1,length(fiber));
for i = 1:length(fiber)
    output_field{i} = GMMNLSE_propagate(fiber{i},input_field,sim{i},gain_rate_eqn);
    t_spent = datevec(output_field{i}.seconds/3600/24);
    fprintf('Running time for %s: %2u:%3.1f\n',model_name{i},t_spent(5),t_spent(6));
end

%% Plot results
energy_GaussianGain  = permute(sum(trapz(abs(output_field{1}.fields).^2),2)*dt/1e3,[3 2 1]);
energy_rategain = permute(sum(trapz(abs(output_field{2}.fields).^2),2)*dt/1e3,[3 2 1]);

% Energy
figure;
plot(output_field{1}.z,[energy_GaussianGain,energy_rategain]);
legend('Gaussian gain','rate-eqn gain');
xlabel('Propagation length (m)');
ylabel('Energy (nJ)');

c = 299792458e-12; % m/ps
f = (-Nt/2:Nt/2-1)'/Nt/dt+c/sim{1}.lambda0;
lambda = c./f*1e9;

c = 299792.458; % nm/ps
factor = c./lambda.^2; % change the spectrum from frequency domain into wavelength domain

% -------------------------------------------------------------------------
% Gaussian gain
% -------------------------------------------------------------------------
% Field
figure;
subplot(2,1,1);
plot(t,abs(output_field{1}.fields(:,:,end)).^2);
legend('mode 1','mode 2','mode 3');
xlabel('Time (ps)');
ylabel('Power (W)');
title('The final output field of YDFA (Gaussian gain)');

% Spectrum
subplot(2,1,2);
plot(lambda,abs(fftshift(ifft(output_field{1}.fields(:,:,end)),1)).^2.*factor);
legend('mode 1','mode 2','mode 3');
xlabel('Wavelength (nm)');
ylabel('PSD (a.u.)');
title('The final output spectrum of YDFA (Gaussian gain)');
xlim([1010 1050]);

% -------------------------------------------------------------------------
% rategain
% -------------------------------------------------------------------------
% Field
figure;
subplot(2,1,1);
plot(t,abs(output_field{2}.fields(:,:,end)).^2);
legend('mode 1','mode 2','mode 3');
xlabel('Time (ps)');
ylabel('Power (W)');
title('The final output field of YDFA (rate-equation gain)');

% Spectrum
subplot(2,1,2);
plot(lambda,abs(fftshift(ifft(output_field{2}.fields(:,:,end)),1)).^2.*factor);
legend('mode 1','mode 2','mode 3');
xlabel('Wavelength (nm)');
ylabel('PSD (a.u.)');
title('The final output spectrum of YDFA (rate-equation gain)');
xlim([1010 1050]);

% =========================================================================
figure;
for i = 1:length(fiber)
    subplot(1,length(fiber),i);
    plot(output_field{i}.z,permute(trapz(abs(output_field{i}.fields).^2)*dt/1e3,[3 2 1]));
    xlabel('Propagation length (m)');
    ylabel('Energy (nJ)');
    ylim([0,14]);
    legend('LP_{01}','LP_{11a}','LP_{11b}');
    switch i
        case 1
            title_string = 'Gaussian gain';
        case 2
            title_string = 'rate-equation gain';
    end
    title(title_string);
end

%% Save results
save('MM_YDFA_10um.mat');