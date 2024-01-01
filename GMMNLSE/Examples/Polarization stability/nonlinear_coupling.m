% This code tests if a field of only one polarization will couple into
% another orthogonal polarization. Linear, circular, and elliptical
% (ellipticity=0.1) polarizations are considered here.
%
% The result shows that the polarization is maintained if the input is
% linear or circular. For a general elliptical polarization, there is
% always coupling between different orthogonal polarization modes.

close all; clearvars;

addpath('../../GMMNLSE algorithm/','../../user_helpers/');

%% Setup fiber parameters
sim.lambda0 = 1025e-9; % the central wavelength
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
[fiber,sim] = load_default_GMMNLSE_propagate(fiber,sim,'single_mode');

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
gain_rate_eqn.reuse_data = false; % For a ring or linear cavity, the pulse will enter a steady state eventually.
                                  % If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
gain_rate_eqn.linear_oscillator = false; % For a linear oscillator, there are pulses from both directions simultaneously, which will deplete the gain;
                                         % therefore, the backward-propagating pulses need to be taken into account.
gain_rate_eqn.t_rep = 24e-9; % Assume 24 MHz here; s; the time required to finish a roundtrip (the inverse repetition rate of the pulse)
                             % This gain model solves the gain of the fiber under the steady-state condition; therefore, the repetition rate must be high compared to the lifetime of the doped ions.
gain_rate_eqn.tau = 840e-6; % lifetime of Yb in F_(5/2) state (Paschotta et al., "Lifetme quenching in Yb-doped fibers"); in "s"
gain_rate_eqn.export_N2 = false; % whether to export N2, the ion density in the upper state or not
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
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

gain_rate_eqn = gain_info( fiber,sim,gain_rate_eqn,ifftshift(lambda,1) );

%% Initial condition
tfwhm = 0.5; % ps
total_energy = 0.1; % nJ
initial_pulse = build_MMgaussian(tfwhm, time_window, total_energy, 1, Nt);
random_polarization = rand(Nt,1)+1i*rand(Nt,1); random_polarization = random_polarization./abs(random_polarization);
initial_pulse.fields = [initial_pulse.fields initial_pulse.fields/100.*random_polarization];

%% Propagate
% linear polarization
sim.ellipticity = 0;
prop_output_linear = GMMNLSE_propagate(fiber,initial_pulse,sim,gain_rate_eqn);

% circular polarization
sim.ellipticity = 1;
prop_output_circular = GMMNLSE_propagate(fiber,initial_pulse,sim,gain_rate_eqn);

% elliptical polarization
sim_e = repmat(sim,9,1);
prop_output_elliptical = repmat(prop_output_linear,9,1);
parfor i = 1:9
    sim_e(i).ellipticity = i/10;
    prop_output_elliptical(i) = GMMNLSE_propagate(fiber,initial_pulse,sim_e(i),gain_rate_eqn);
end

%% Results
energy_linear = permute(trapz(abs(prop_output_linear.fields).^2,1)*dt/1e3,[3 2 1]);
energy_circular = permute(trapz(abs(prop_output_circular.fields).^2,1)*dt/1e3,[3 2 1]);
energy_elliptical = cell(9,1);
for i = 1:9
    energy_elliptical{i} = permute(trapz(abs(prop_output_elliptical(i).fields).^2,1)*dt/1e3,[3 2 1]);
end

[phi_linear,theta_linear] = calc_ellipticity( prop_output_linear.fields(:,:,end),0);
[phi_circular,theta_circular] = calc_ellipticity( prop_output_circular.fields(:,:,end),1);
phi_elliptical = cell(9,1); theta_elliptical = cell(9,1);
for i = 1:9
    [phi_elliptical{i},theta_elliptical{i}] = calc_ellipticity( prop_output_elliptical(i).fields(:,:,end),sim_e(i).ellipticity);
end

% Transform into the circular basis
P = @(r)[1  1i*r;...
         r   -1i]/sqrt(1+r^2); % the transformation matrix from basis (x,y) into one with ellipticity "r"
P_circ = @(r) P(1)/P(r); % finally transform into the circular basis

prop_output_linear_circ = prop_output_linear;
prop_output_elliptical_circ = prop_output_elliptical;
prop_output_linear_circ.fields = gather(permute(pagefun(@mtimes, P(1),gpuArray(permute(prop_output_linear.fields,[2 1 3]))),[2 1 3]));
for i = 1:9
    prop_output_elliptical_circ(i).fields = gather(permute(pagefun(@mtimes, P_circ(sim_e(i).ellipticity),gpuArray(permute(prop_output_elliptical(i).fields,[2 1 3]))),[2 1 3]));
end

% Under circular basis
energy_linear_circ = permute(trapz(abs(prop_output_linear_circ.fields).^2,1)*dt/1e3,[3 2 1]);
ratio_linear_circ = energy_linear_circ(:,1)./energy_linear_circ(:,2);
energy_elliptical_circ = cell(9,1);
ratio_elliptical_circ = cell(9,1);
for i = 1:9
    energy_elliptical_circ{i} = permute(trapz(abs(prop_output_elliptical_circ(i).fields).^2,1)*dt/1e3,[3 2 1]);
    ratio_elliptical_circ{i} = energy_elliptical_circ{i}(:,1)./energy_elliptical_circ{i}(:,2);
end
ratio_circular = energy_circular(:,1)./energy_circular(:,2);

% total_field
total_field_linear = sum(abs(prop_output_linear.fields(:,:,end)).^2,2);
total_field_circular = sum(abs(prop_output_circular.fields(:,:,end)).^2,2);
total_field_elliptical = cell(9,1);
for i = 1:9
    total_field_elliptical{i} = sum(abs(prop_output_elliptical(i).fields(:,:,end)).^2,2);
end

%% Plot Results
c = 299792458; % m/s

% linear polarization
fig = figure('Name','linear');
fp = get(fig,'position');
screen_size = get(0,'ScreenSize');
original_top = screen_size(4)-fp(2)-fp(4);
set(fig,'position',[fp(1) screen_size(4)-original_top-fp(4)*7/4 fp(3)*5/4 fp(4)*7/4]);
subplot(4,3,1);
h = plot(prop_output_linear.z,energy_linear(:,1));
set(h,'linewidth',2);
xlabel('Propagation distance (m)'); ylabel('Energy (nJ)');
title('x polarization');
subplot(4,3,4);
h = plot(prop_output_linear.z,energy_linear(:,2));
set(h,'linewidth',2);
xlabel('Propagation distance (m)'); ylabel('Energy (nJ)');
title('y polarization');
subplot(4,3,2);
h = plot(prop_output_linear.z,energy_linear_circ(:,1));
set(h,'linewidth',2);
xlabel('Propagation distance (m)'); ylabel('Energy (nJ)');
title('\sigma_+ polarization');
subplot(4,3,5);
h = plot(prop_output_linear.z,energy_linear_circ(:,2));
set(h,'linewidth',2);
xlabel('Propagation distance (m)'); ylabel('Energy (nJ)');
title('\sigma_- polarization');
subplot(4,3,[3 6]);
h = plot(prop_output_linear.z,ratio_linear_circ);
set(h,'linewidth',2);
xlabel('Propagation distance (m)'); ylabel('Ratio');
title('Energy ratio (\sigma_+/\sigma_-)');
subplot(4,3,7:9);
yyaxis left
plot(t,abs(prop_output_linear_circ.fields(:,:,end)).^2);
hold on; plot(t,total_field_linear,'r'); hold off;
ylabel('Intensity (W)');
legend('\sigma_+','\sigma_-','total field');
yyaxis right
h = plot(t,[phi_linear 270*ones(Nt,1) 90*ones(Nt,1)]);
set(h(1),'linewidth',2);
xlim([-10 0]);
xlabel('time (ps)'); ylabel('the phase difference (deg)');
title('Final fields');
subplot(4,3,10:12);
spectrum_linear = abs(fftshift(ifft(prop_output_linear_circ.fields(:,:,end)),1)).^2;
h = plot(lambda(lambda>0),spectrum_linear(lambda>0,:)*(Nt*dt)^2*c/1e6./lambda(lambda>0).^2);
set(h,'linewidth',2);
xlim([950,1200]);
ylabel('Intensity (nJ/nm)');
xlabel('Wavelength (nm)');
title('Spectrum');

% circular polarization
fig = figure('Name','circular');
fp = get(fig,'position');
screen_size = get(0,'ScreenSize');
original_top = screen_size(4)-fp(2)-fp(4);
set(fig,'position',[fp(1) screen_size(4)-original_top-fp(4)*7/4 fp(3)*5/4 fp(4)*7/4]);
subplot(4,2,1);
h = plot(prop_output_circular.z,energy_circular(:,1));
set(h,'linewidth',2);
xlabel('Propagation distance (m)'); ylabel('Energy (nJ)');
title('\sigma_+ polarization');
subplot(4,2,3);
h = plot(prop_output_circular.z,energy_circular(:,2));
set(h,'linewidth',2);
xlabel('Propagation distance (m)'); ylabel('Energy (nJ)');
title('\sigma_- polarization');
subplot(4,2,[2 4]);
h = plot(prop_output_linear.z,ratio_circular);
set(h,'linewidth',2);
xlabel('Propagation distance (m)'); ylabel('Ratio');
title('Energy ratio (\sigma_+/\sigma_-)');
subplot(4,2,5:6);
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
subplot(4,2,7:8);
spectrum_circular = abs(fftshift(ifft(prop_output_circular.fields(:,:,end)),1)).^2;
h = plot(lambda(lambda>0),spectrum_circular(lambda>0,:)*(Nt*dt)^2*c/1e6./lambda(lambda>0).^2);
set(h,'linewidth',2);
xlim([950,1200]);
ylabel('Intensity (nJ/nm)');
xlabel('Wavelength (nm)');
title('Spectrum');

% elliptical polarization
for i = 1:9
    fig = figure('Name',sprintf('ellipticity: %2.1f',sim_e(i).ellipticity));
    fp = get(fig,'position');
    screen_size = get(0,'ScreenSize');
    original_top = screen_size(4)-fp(2)-fp(4);
    set(fig,'position',[fp(1) screen_size(4)-original_top-fp(4)*7/4 fp(3)*5/4 fp(4)*7/4]);
    subplot(4,3,1);
    h = plot(prop_output_linear.z,energy_elliptical{i}(:,1));
    set(h,'linewidth',2);
    xlabel('Propagation distance (m)'); ylabel('Energy (nJ)');
    title('the 1st polarization');
    subplot(4,3,4);
    h = plot(prop_output_linear.z,energy_elliptical{i}(:,2));
    set(h,'linewidth',2);
    xlabel('Propagation distance (m)'); ylabel('Energy (nJ)');
    title('the 2nd polarization');
    subplot(4,3,2);
    h = plot(prop_output_linear.z,energy_elliptical_circ{i}(:,1));
    set(h,'linewidth',2);
    xlabel('Propagation distance (m)'); ylabel('Energy (nJ)');
    title('\sigma_+ polarization');
    subplot(4,3,5);
    h = plot(prop_output_linear.z,energy_elliptical_circ{i}(:,2));
    set(h,'linewidth',2);
    xlabel('Propagation distance (m)'); ylabel('Energy (nJ)');
    title('\sigma_- polarization');
    subplot(4,3,[3 6]);
    h = plot(prop_output_linear.z,ratio_elliptical_circ{i});
    set(h,'linewidth',2);
    xlabel('Propagation distance (m)'); ylabel('Ratio');
    title('Energy ratio (\sigma_+/\sigma_-)');
    subplot(4,3,7:9);
    yyaxis left
    plot(t,abs(prop_output_elliptical_circ(i).fields(:,:,end)).^2);
    hold on; plot(t,total_field_elliptical{i},'r'); hold off;
    ylabel('Intensity (W)');
    legend('\sigma_+','\sigma_-','total field');
    yyaxis right
    h = plot(t,[phi_elliptical{i} 270*ones(Nt,1) 90*ones(Nt,1)]);
    set(h(1),'linewidth',2);
    xlim([-10 0]);
    xlabel('time (ps)'); ylabel('the phase difference (deg)');
    title('Final fields');
    subplot(4,3,10:12);
    spectrum_elliptical_circ = abs(fftshift(ifft(prop_output_elliptical_circ(i).fields(:,:,end)),1)).^2;
    h = plot(lambda(lambda>0),spectrum_elliptical_circ(lambda>0,:)*(Nt*dt)^2*c/1e6./lambda(lambda>0).^2);
    set(h,'linewidth',2);
    xlim([950,1200]);
    ylabel('Intensity (nJ/nm)');
    xlabel('Wavelength (nm)');
    title('Spectrum');
end

%% Save the result
save('nonlinear_gain_coupling.mat','t','Nt','dt','lambda','sim','sim_e','prop_output_linear','prop_output_circular','prop_output_elliptical');