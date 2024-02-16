% This code demonstrates the generation of multimode soliton through SSFS.
% Strong birefringence is added to demonstrate the polarized simulations.
% Due to "strong" birefringence, its result is the same as the scalar
% simulation.

close all; clearvars;

%% Load folders
addpath('../../GMMNLSE algorithm/','../../user_helpers/');

% Add multimode fiber informaiton
fiber.MM_folder = '../../Fibers/GRIN_62.5um_wavelength1550nm/';
fiber.betas_filename = 'betas.mat';
fiber.S_tensors_filename = 'S_tensors_10modes.mat';

%% Setup fiber parameters
sim.lambda0 = 1550e-9; % the central wavelength
sim.pulse_centering = false;
sim.midx = 1:2; % use two spatial modes for fast demonstration; users are free to try more
sim.scalar = false;
sim.gpu_yes = true; % make it "false" to test CPU

num_modes = length(sim.midx)*(double(~sim.scalar)+1); % 2, in the polarized scenario, accounts for polarized modes

% Load default parameters like 
%
% loading fiber.betas and fiber.SR based on your multimode folder above
% sim.Raman_model = 1; Use isotropic Raman model
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.gpu_yes = true; Use GPU (default to true)
% ......
%
% Please check this function for details.
[fiber,sim] = load_default_GMMNLSE_propagate(fiber,sim,'multimode'); % load default parameters

% Add birefringence to the fiber
beat_length = 1e-3; % 1 mm
% The betas at the other polarization can be easily calculated by 
% beta_e = n_e*w/c = (n_o + delta_n)*w/c
%                  = (n_o + lambda/beat_length)*w/c
if ~sim.scalar
    betas = zeros(size(fiber.betas,1),num_modes);
    for i = 1:length(sim.midx)
        betas(:,i*2-1) = fiber.betas(:,i);
        betas(1,i*2) = fiber.betas(1,i) + 980e-9/beat_length *2*pi/sim.lambda0;
    end
    fiber.betas = betas;
end

num_save = 100;
fiber.L0 = 10; % m
sim.save_period = fiber.L0/num_save;

%% Setup general parameters
Nt = 2^11; % the number of time points
time_window = 50; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

%% Initial condition
total_energy = 5; % nJ
tfwhm = 0.5; % ps

% Initially excited in the same polarization only
ratio = ones(1,num_modes); 
if num_modes > 1
    ratio(2:2:end) = ratio(2:2:end)/1e3;
end
input_field = build_MMgaussian(tfwhm, time_window, total_energy, num_modes, Nt, {'ifft',0}, ratio, -time_window/2*0.4);
input_field.fields = input_field.fields.*exp(1i*2*pi*rand(1,num_modes));

%% Propagate
prop_output = GMMNLSE_propagate(fiber,input_field,sim);

%% Plot
% Time
figure;
plot(t,abs(prop_output.fields(:,:,end)).^2,'linewidth',2);
xlabel('Time (ps)');
ylabel('Power (W)');
set(gca,'fontsize',14);

% Spectrum
figure;
plot(lambda,abs(fftshift(ifft(prop_output.fields(:,:,end)),1)).^2*c./lambda.^2,'linewidth',2);
xlabel('\lambda (nm)');
ylabel('PSD (a.u.)');
xlim([1500,1700]);
set(gca,'fontsize',14);

% Comparison of time
figure;
[x,y] = meshgrid(t,prop_output.z);
pcolor(x,y,permute(abs(prop_output.fields(:,1,:)).^2,[3 1 2]));
shading interp; colormap(jet);
xlabel('Time (ps)');
ylabel('Propagation distance (m)');
set(gca,'fontsize',14);

% Comparison of spectra
figure;
[x,y] = meshgrid((f-sim.f0),prop_output.z(2:end));
tmp = 10*log10(permute(abs(fftshift(ifft(prop_output.fields(:,1,2:end)),1)).^2,[3 1 2])); tmp = tmp - max(tmp(:));
pcolor(x,y,tmp);
shading interp; colormap(jet); caxis([-20,0]);
xlabel('\Deltaf (THz)');
ylabel('Propagation distance (m)');
set(gca,'fontsize',14);