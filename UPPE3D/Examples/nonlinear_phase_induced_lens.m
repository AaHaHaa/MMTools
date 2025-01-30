% This code simulates the focusing effect induced by the nonlinear phase

close all; clearvars;

addpath('../UPPE3D algorithm/','../user_helpers/');

%% Setup plate parameters
sim.lambda0 = 1030e-9; % the center wavelength
sim.gpuDevice.Index = 2;

[fiber,sim] = load_default_UPPE3D_propagate([],sim); % load default parameters

num_save = 10;

%% Initial condition
MFD0 = 1000e-6; % m
spatial_window = 4e-3; % m
tfwhm = 0.08; % ps
time_window = 1; % ps
energy = 100e3; % nJ
Nt = 2^7; % the number of time points
Nx = 2^9; % the number of spatial points
initial_condition = build_3Dgaussian(MFD0, spatial_window, tfwhm, time_window, energy, Nt, Nx);

%% Setup general parameters
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm
dx = spatial_window/Nx; % m
x = (-Nx/2:Nx/2-1)*dx;

%% Material properties
% Sellmeier coefficients
[a,b] = Sellmeier_coefficients('sapphire');
% Calculate the index difference using the Sellmeier equation to generate n(lambda)
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));

% Plate material properties
n_plate = n_from_Sellmeier(lambda/1e3).*ones(1,Nx,Nx); % refractive index
n2_plate = 3e-20; % nonlinear refractive index

% Air properties
n_air = ones(1,Nx,Nx); % air's n=1
n2_air = 0; % no nonlinearity

%% Plot the initial field
% Show initial real space
figure;
pcolor(x*1e6,x*1e6,abs(squeeze(initial_condition.field(Nt/2,:,:))).^2); colormap(jet);colorbar;
shading interp;colormap(jet);colorbar;
xlabel('x (\mum)');
xlim([-20,20]*1e2);
ylim([-20,20]*1e2);
title('Initial real space');
% Show initial k space
figure;
pcolor(abs(fftshift(fft(fft(squeeze(initial_condition.field(Nt/2,:,:)),[],1),[],2))).^2); colormap(jet);colorbar;
shading interp;colormap(jet);colorbar;
title('Initial k space');

%% 1. After the lens
focal_length = 0.0291; % m
Ef_out = add_thin_lens_phase(ifft(initial_condition.field,[],1),initial_condition.dx,initial_condition.dy,fftshift(lambda,1)/1e9,focal_length);

% free space after the lens (in air)
initial_condition.field = fft(Ef_out(:,:,:,end),[],1);

%% 2. Propagate over the 1-mm plate
fiber_plate = fiber;
fiber_plate.n = n_plate;
fiber_plate.n2 = n2_plate;
fiber_plate.L0 = 1e-3; % m
sim.save_period = fiber_plate.L0/num_save;
prop_output_plate = UPPE3D_propagate(fiber_plate,initial_condition,sim);

%% 3. Propagate up to the focal point
fiber_air = fiber;
fiber_air.n = n_air;
fiber_air.n2 = n2_air; % assume no nonlinearity in air
fiber_air.L0 = focal_length*1.1 - fiber_plate.L0; % m
sim.save_period = fiber_air.L0/num_save;
prop_output_air = UPPE3D_propagate(fiber_air,prop_output_plate,sim);

%% Results
z = [prop_output_plate.z;prop_output_plate.z(end)+prop_output_air.z(2:end)]; % m
field = cat(4,prop_output_plate.field,prop_output_air.field(:,:,:,2:end));

MFD = squeeze(calcMFD(squeeze(field(Nt/2,:,:,:)),spatial_window))*1e6; % um

energy3D_plate = squeeze(sum(abs(prop_output_plate.field).^2,[1,2,3]));
energy3D_air2 = squeeze(sum(abs(prop_output_air.field).^2,[1,2,3]));

%% Plot
% Show final real space
figure;
pcolor(x*1e6,x*1e6,abs(squeeze(prop_output_air.field(Nt/2,:,:,end))).^2); colormap(jet);colorbar;
shading interp;colormap(jet);colorbar;
xlabel('x (\mum)');
xlim([-1,1]*1e2);
ylim([-1,1]*1e2);
title('Final real space');
% Show final k space
figure;
pcolor(abs(fftshift(fft(fft(squeeze(prop_output_air.field(Nt/2,:,:,end)),[],1),[],2))).^2); colormap(jet);colorbar;
shading interp;colormap(jet);colorbar;
title('Final k space')

% Plot MFD
figure;
plot(z*1e2,MFD,'linewidth',2);
xlabel('Propagation distance (cm)');
ylabel('MFD (\mum)');
set(gca,'fontsize',20);