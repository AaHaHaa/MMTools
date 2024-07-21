% Example propagation in GRIN fiber, at a relatively low power.
% You can use this as a starting point for more specific simulations

close all; clearvars;

addpath('../../../GMMNLSE algorithm/','../../../user_helpers/'); % MATLAB needs to know where the propagate files are located

%% Setup fiber parameters
num_modes = 6; % The number of modes to be simulated. 1-~8 will be "very fast", 10-20 will get slower but it still works. More than 20 modes have not been tested

% Multimode parameters
fiber.MM_folder = '../../../Fibers/OM4_wavelength1030nm/';
fiber.betas_filename = 'betas.mat';
fiber.S_tensors_filename = ['S_tensors_' num2str(num_modes) 'modes.mat'];

fiber.L0 = 0.15; % Fiber length in m

%% Setup simulation parameters
c = 299792458*1e-12; % m/ps
wavelength_range = [0.7,2]*1e-6; % m
Nt = 2^11;
sim.lambda0 = 1030e-9; % central pulse wavelength (m
time_window = 7; % ps
dt = time_window/Nt;
t = (-Nt/2:Nt/2-1)'*dt;

[fiber,sim] = load_default_GMMNLSE_propagate(fiber,sim,'multimode');

f = sim.f0 + (-Nt/2:Nt/2-1)'/time_window; % THz

%% Setup initial conditions
tfwhm = 0.05; % ps, FWHM of the initial pulse.
total_energy = 17.04; % nJ, total energy of the initial pulse. By convension this is the total energy in all modes

% This is a helper function to build an evently distributed gaussian
% initial MM pulse
pump_wavelength = 1030e-9; % m
freq_shift = c/pump_wavelength - sim.f0;
initial_condition = build_MMgaussian(tfwhm, time_window, total_energy, num_modes, Nt, {'ifft',freq_shift});

%% Run the propagation
prop_output = GMMNLSE_propagate(fiber, initial_condition, sim); % This actually does the propagation

%% Plot the results

% Plot the time domain
figure();
I_time = abs(prop_output.fields(:, :, end).^2);
tlim = 1;

subplot(1, 2, 1);
plot(t, I_time),axis tight, grid on
ylabel('Power (W)')
xlabel('Time (ps)')
xlim([-tlim, tlim])

% Plot the frequency domain
I_freq = abs(ifftshift(ifft(prop_output.fields(:, :, end)))).^2;
flim = 60;

subplot(1, 2, 2);
plot(f, I_freq),axis tight, grid on
ylabel('PSD (a.u.)')
xlabel('Frequency (THz)')
xlim([sim.f0-flim, sim.f0+flim])

%% Load the spatial modes and plot the full spatial field

% Load the modes
prefix = '../../../Fibers/OM4_wavelength1030nm';
Nx = 400; % The number of spatial grid points that the modes use
mode_profiles = zeros(Nx, Nx, num_modes);
radius = '25'; % Used for loading the file
lambda0 = '10300'; % Used for loading the file
for ii = 1:num_modes
   name = [prefix, '/mode',int2str(ii),'wavelength', lambda0, '.mat'];
   load(name, 'phi');
   mode_profiles(:, :, ii) = phi; % Save the modes
   disp(['Loaded mode ', int2str(ii)])
end
load(name, 'x');
x = (x-mean(x))*1e-6; % The spatial coordinates along one dimension
dx = x(2)-x(1);
mode_profiles = mode_profiles./sqrt(sum(sum(abs(mode_profiles).^2,1),2))/dx;

% Downsample in space to reduce memory usage
factor = 8;
dx = dx*factor;
mode_profiles_sampled = zeros(Nx/factor, Nx/factor, num_modes);
for ii = 1:num_modes
    mode_profiles_sampled(:, :, ii) = downsample(downsample(mode_profiles(:, :, ii), factor)', factor)';
end
x = downsample(x, factor);
Nx = Nx/factor;
[X, Y] = meshgrid(x, x);

% Build the field from the modes and the spatial profiles
E_txy = recompose_into_space(sim.gpu_yes, mode_profiles_sampled, prop_output.fields(:, :, end), sim.cuda_dir_path);
A0 = permute(sum(abs(E_txy).^2, 1)*prop_output.dt/1e12, [2 3 1]); % Integrate over time to get the average spatial field

% Plot the spatial field
figure();
h = pcolor(X*1e6, Y*1e6, A0);
h.LineStyle = 'none';
colorbar;
axis square;
xlabel('x (um)');
ylabel('y (um)');
xlim([-60, 60]);
ylim([-60, 60]);