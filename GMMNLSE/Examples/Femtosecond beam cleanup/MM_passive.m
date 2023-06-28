clearvars; close all;

%% Add the folders of multimode files and others
addpath('../../GMMNLSE algorithm/','../../user_helpers/'); % add where many GMMNLSE-related functions like  "GMMNLSE_propagate" is
fiber.MM_folder = '../../Fibers/GRIN_168_400_wavelength1030nm/';
fiber.betas_filename = 'betas.mat';
fiber.S_tensors_filename = 'S_tensors_21modes.mat';

%% Field and simulation parameters
time_window = 100; % ps
Nt = 2^10; % the number of time points
dt = time_window/Nt;
t = (-Nt/2:Nt/2-1)'*dt; % ps

fiber.L0 = 1; % m; the length of the gain fiber
save_num = 50;
sim.save_period = fiber.L0/save_num;
sim.lambda0 = 1030e-9; % central wavelength; in "m"
%sim.progress_bar = false;
[fiber,sim] = load_default_GMMNLSE_propagate(fiber,sim,'multimode');

%% Random mode coupling
addpath('../../Random_Mode_Coupling/');

sim.rmc.model = true;

% weaker coupling for spatial modes
sim.rmc.varn = 1e-4; % variance of the refractive index
num_modes = length(sim.midx);
sim.deltaZ = 1e-3; % m
save_points = int32(fiber.L0/sim.deltaZ);
sim.rmc.matrices = create_rmc_matrices(fiber,sim,num_modes,save_points);
effective_matrix = calc_effective_rmc_matrix(fiber,sim,Nt,dt,...
                                             sim.rmc.matrices);
fraction = effective_matrix\[1;zeros(num_modes-1,1)];

%% Initial pulse
total_energy = 1000; % nJ
tfwhm = 10; % ps
input_field = build_MMgaussian(tfwhm, time_window, total_energy, length(sim.midx), Nt, {'ifft',0}, [ones(1,5),ones(1,length(sim.midx)-5)*0.0001]);
rand_phase = rand(1,length(sim.midx));
input_field.fields = input_field.fields.*exp(1i*2*pi*rand_phase);

%% Propagation
output_field = GMMNLSE_propagate(fiber,input_field,sim);
%t_spent = datevec(output_field.seconds/3600/24);
%fprintf('Running time for %s: %2u:%3.1f\n','rate-eqn gain',t_spent(5),t_spent(6));

%% Plot results
energy = permute(sum(trapz(abs(output_field.fields).^2),2)*dt/1e3,[3 2 1]);

% Energy
distance = (0:save_num)*sim.save_period;
figure;
plot(distance,energy);
xlabel('Propagation length (m)');
ylabel('Energy (nJ)');
title('Energy');

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
plot(output_field.z,permute(trapz(abs(output_field.fields).^2)*dt/1e3,[3 2 1]));
xlabel('Propagation length (m)');
ylabel('Energy (nJ)');
title('rate-equation gain');

figure;
plot(output_field.z,each_energy(:,1)./sum(each_energy,2));