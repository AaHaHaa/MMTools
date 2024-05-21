% This code runs the single-mode Yb-doped fiber amplifier with the gain 
% rate equation and studies the destabilization of ASE to the pulse.
% This occurs when the seed is too weak with a strong pump power.
%
% If the ASE is too strong, due to its noisy nature, the iterations in the
% rate-eqn-gain model will not converge. The simulation will stop after 30
% iterations.
% Not only the ASE powers but also the pulse energy don't converge during
% iterations.
% We have included ASE into the pulse, so a strong noisy ASE makes the
% pulse noisy as well.

clearvars; close all;

addpath('../../GMMNLSE algorithm/','../../user_helpers/');

%% Gain info
% Please find details of all the parameters in "gain_info.m".
% Note that the use of single spatial mode is different from multi-spatial modes.
gain_rate_eqn.gain_medium = 'Yb'; % specify the gain medium
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
gain_rate_eqn.copump_power = 2; % W
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.t_rep = 1/10e6; % assume 5 MHz here; s; the time required to finish a roundtrip (the inverse repetition rate of the pulse)
                             % This gain model solves the gain of the fiber under the steady-state condition; therefore, the repetition rate must be high compared to the lifetime of the doped ions.
gain_rate_eqn.export_N2 = false; % whether to export N2, the ion density in the upper state or not
gain_rate_eqn.ignore_ASE = false;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 5; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the iteration
gain_rate_eqn.verbose = true; % show the information(final pulse energy) during iterations

%% Field and simulation parameters
time_window = 50; % ps
Nt = 2^13; % the number of time points
dt = time_window/Nt;
t = (-Nt/2:Nt/2-1)'*dt; % ps

fiber.L0 = 2; % m; the length of fiber length
save_num = 100; % the number of saved data
sim.save_period = fiber.L0/save_num;
sim.gpu_yes = false;

sim.lambda0 = 1030e-9; % central wavelength; in "m"

% Rate-equation gain
sim.gain_model = 2;
[fiber,sim] = load_default_GMMNLSE_propagate(fiber,sim);

%% Initial pulse
total_energy = 0.1; % nJ
tfwhm = 1; % ps
input_field = build_MMgaussian(tfwhm, time_window, total_energy, 1, Nt);
input_field.Power.ASE.forward = zeros(size(input_field.fields));
input_field.Power.ASE.backward = zeros(size(input_field.fields));

%% Gain parameters
% Precompute some parameters related to the gain to save the computational time
% Check "gain_info.m" for details.
f = ifftshift( (-Nt/2:Nt/2-1)'/Nt/dt + sim.f0 ); % in the order of "omegas" in the "GMMNLSE_propagate.m"
c = 299792.458; % nm/ps;
lambda = c./f; % nm
gain_rate_eqn_with_ASE = gain_info( fiber,sim,gain_rate_eqn,lambda );

% Prepare for the one without ASE
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn_without_ASE = gain_info( fiber,sim,gain_rate_eqn,lambda );

%% Propagations
output_with_ASE = GMMNLSE_propagate(fiber,input_field,sim,gain_rate_eqn_with_ASE);

% Run another case without ASE
output_without_ASE = GMMNLSE_propagate(fiber,input_field,sim,gain_rate_eqn_without_ASE);

%% Plot results
energy_with_ASE = permute(sum(trapz(abs(output_with_ASE.fields).^2),2)*dt/1e3,[3 2 1]);
energy_without_ASE = permute(sum(trapz(abs(output_without_ASE.fields).^2),2)*dt/1e3,[3 2 1]);

% Energy
distance = (0:save_num)*sim.save_period;
figure;
h = plot(distance,[energy_with_ASE,energy_without_ASE],'linewidth',2');
set(h(1),'Color','b'); set(h(2),'Color','r');
l = legend('with ASE','without ASE'); set(l,'location','northwest');
xlabel('Propagation distance (m)');
ylabel('Energy (nJ)');
set(gca,'fontsize',20);

c = 299792458e-12; % m/ps
f = (-Nt/2:Nt/2-1)'/Nt/dt+c/sim.lambda0;
lambda = c./f*1e9;

c = 299792.458; % nm/ps
factor = c./lambda.^2; % change the spectrum from frequency domain into wavelength domain
factor_correct_unit = time_window^2/1e3; % to make the spectrum of the correct unit "nJ/THz"
                                         % "/1e3" is to make pJ into nJ

% -------------------------------------------------------------------------
% with ASE
% -------------------------------------------------------------------------
% Field
figure;
plot(t,abs(output_with_ASE.fields(:,:,end)).^2/1e3,'linewidth',2,'Color','b');
xlabel('Time (ps)');
ylabel('Power (kW)');
set(gca,'fontsize',20);

% Spectrum
figure;
plot(lambda,abs(fftshift(ifft(output_with_ASE.fields(:,:,end)),1)).^2.*factor*factor_correct_unit,'linewidth',2,'Color','r');
xlabel('Wavelength (nm)');
ylabel('PSD (nJ/nm)');
xlim([850 1350]);
set(gca,'fontsize',20);

% -------------------------------------------------------------------------
% without ASE
% -------------------------------------------------------------------------
% Field
figure;
plot(t,abs(output_without_ASE.fields(:,:,end)).^2/1e3,'linewidth',2,'Color','b');
xlabel('Time (ps)');
ylabel('Power (kW)');
set(gca,'fontsize',20);

% Spectrum
figure;
plot(lambda,abs(fftshift(ifft(output_without_ASE.fields(:,:,end)),1)).^2.*factor*factor_correct_unit,'linewidth',2,'Color','r');
xlabel('Wavelength (nm)');
ylabel('PSD (nJ/nm)');
xlim([850 1350]);
set(gca,'fontsize',20);

%% ASE
func = analyze_sim;
func.analyze_ASE(f,output_with_ASE.Power.ASE,output_with_ASE.z);