close all; clearvars;

addpath('../../GMMNLSE algorithm/','../../user_helpers/');
fiber.MM_folder = '1060XP_wavelength1030nm/';
fiber.betas_filename = 'betas.mat';
fiber.S_tensors_filename = 'S_tensors_1modes.mat';

%% Setup fiber parameters
sim.lambda0 = 1030e-9; % the central wavelength
sim.scalar = false; % polarized fields
sim.gain_model = 2;
fiber.L0 = 4; % m
sim.save_period = fiber.L0/100;

% Load default parameters like 
%
% loading fiber.betas and fiber.SR based on your multimode folder above
% sim.Raman_model = 1; Use isotropic Raman model
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.gpu_yes = true; Use GPU (default to true)
% ......
%
% Please check this function for details.
[fiber,sim] = load_default_GMMNLSE_propagate(fiber,sim,'multimode');

%% Gain info
% Please find details of all the parameters in "gain_info.m".
% Note that the usage of single spatial mode is different from multi-spatial modes.
gain_rate_eqn.cross_section_filename = 'Liekki Yb_AV_20160530.txt';
gain_rate_eqn.core_diameter = 6; % um
gain_rate_eqn.cladding_diameter = 125; % um
gain_rate_eqn.core_NA = 0.12;
gain_rate_eqn.absorption_wavelength_to_get_N_total = 920; % nm
gain_rate_eqn.absorption_to_get_N_total = 0.55; % dB/m
gain_rate_eqn.pump_wavelength = 976; % nm
gain_rate_eqn.copump_power = 4; % W
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.downsampling_factor = 1; % an integer; downsample the eigenmode profiles to run faster
gain_rate_eqn.reuse_data = false; % For a ring or linear cavity, the pulse will enter a steady state eventually.
                                  % If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
gain_rate_eqn.linear_oscillator = false; % For a linear oscillator, there are pulses from both directions simultaneously, which will deplete the gain;
                                         % therefore, the backward-propagating pulses need to be taken into account.
gain_rate_eqn.t_rep = 24e-9; % Assume 24 MHz here; s; the time required to finish a roundtrip (the inverse repetition rate of the pulse)
                             % This gain model solves the gain of the fiber under the steady-state condition; therefore, the repetition rate must be high compared to the lifetime of the doped ions.
gain_rate_eqn.tau = 840e-6; % lifetime of Yb in F_(5/2) state (Paschotta et al., "Lifetme quenching in Yb-doped fibers"); in "s"
gain_rate_eqn.export_N2 = false; % whether to export N2, the ion density in the upper state or not
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.max_iterations = 50; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the above iterations
gain_rate_eqn.verbose = false; % show the information(final pulse energy) during iterations

%% Setup general parameters
Nt = 2^13; % the number of time points
time_window = 50; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

%% Set up betas here because I computed them as a function of frequency in betas, instead of Taylor-series coefficients
fiber.betas = myPchip(fiber.beta_f,fiber.betas,f,6,sim.cuda_dir_path);

%% Initial condition
tfwhm = 0.5; % ps
total_energy = 0.1; % nJ
initial_pulse = build_MMgaussian(tfwhm, time_window, total_energy, 1, Nt);
random_polarization = rand(Nt,1)+1i*rand(Nt,1); random_polarization = random_polarization./abs(random_polarization);
initial_pulse.fields = [initial_pulse.fields initial_pulse.fields/100.*random_polarization];

%% Random mode coupling
addpath('../../Random_Mode_Coupling/');

sim.rmc.model = true;

% weaker coupling for spatial modes
sim.rmc.stdQ_polarizedmode = 10; % the coupling strength between polarization modes
sim.deltaZ = 100e-6; % m
save_points = int32(fiber.L0/sim.deltaZ);
sim.rmc.matrices = create_rmc_matrices(fiber,sim,2,save_points);
effective_matrix = calc_effective_rmc_matrix(fiber,sim,Nt,dt,...
                                             sim.rmc.matrices,...
                                             initial_pulse.fields);
fraction = effective_matrix\[1;0];

% Update gain_rate_eqn based on some settings
% sim.rmc.model affects it too, so it needs to be after sim.rmc.model is
% set.
gain_rate_eqn = gain_info( fiber,sim,gain_rate_eqn,ifftshift(lambda,1) );

%% Propagate
% circular polarization
sim.ellipticity = 1;
prop_output_circular = GMMNLSE_propagate(fiber,initial_pulse,sim,gain_rate_eqn);

%% Results
energy_circular = permute(trapz(abs(prop_output_circular.fields).^2,1)*dt/1e3,[3 2 1]);

[phi_circular,theta_circular] = calc_ellipticity( prop_output_circular.fields(:,:,end),1);

% total_field
total_field_circular = sum(abs(prop_output_circular.fields(:,:,end)).^2,2);

%% Plot Results
c = 299792458; % m/s

% circular polarization
fig = figure('Name','circular');
fp = get(fig,'position');
screen_size = get(0,'ScreenSize');
original_top = screen_size(4)-fp(2)-fp(4);
set(fig,'position',[fp(1) screen_size(4)-original_top-fp(4)*7/4 fp(3)*5/4 fp(4)*7/4]);
subplot(3,2,1);
h = plot(prop_output_circular.z,energy_circular(:,1));
set(h,'linewidth',2);
xlabel('Propagation distance (m)'); ylabel('Energy (nJ)');
title('\sigma_+ polarization');
subplot(4,2,2);
h = plot(prop_output_circular.z,energy_circular(:,2));
set(h,'linewidth',2);
xlabel('Propagation distance (m)'); ylabel('Energy (nJ)');
title('\sigma_- polarization');
subplot(3,2,3:4);
yyaxis left
plot(t,abs(prop_output_circular.fields(:,:,end)).^2);
hold on; plot(t,total_field_circular,'r'); hold off;
ylabel('Intensity (W)');
legend('\sigma_+','\sigma_-','total field');
yyaxis right
h = plot(t,[phi_circular 270*ones(Nt,1) 90*ones(Nt,1)]);
set(h(1),'linewidth',2);
xlim([-10 0]);
xlabel('time (ps)'); ylabel('the phase difference (deg)');
title('Final fields');
subplot(3,2,5:6);
spectrum_circular = abs(fftshift(ifft(prop_output_circular.fields(:,:,end)),1)).^2;
h = plot(lambda(lambda>0),spectrum_circular(lambda>0,:)*(Nt*dt)^2*c/1e6./lambda(lambda>0).^2);
set(h,'linewidth',2);
xlim([950,1200]);
ylabel('Intensity (nJ/nm)');
xlabel('Wavelength (nm)');
title('Spectrum');

%% Save the result
save('circular_stability.mat','t','Nt','dt','lambda','sim','prop_output_circular');