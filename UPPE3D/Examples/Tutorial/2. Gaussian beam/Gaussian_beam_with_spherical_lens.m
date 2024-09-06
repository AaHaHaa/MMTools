% This code simulates a Gaussian beam passing through a thick lens (in air).
% The simulated MFD is compared to the theoretical values.

close all; clearvars;

addpath('../../../UPPE3D algorithm/','../../../user_helpers/');

%% Setup fiber parameters
sim.lambda0 = 1030e-9; % the center wavelength

% Load default parameters like
%
% load fiber.betas and fiber.SR based on your multimode folder above
% sim.f0 = 3e5/sim.lambda0; THz
% sim.save_period = 0; Save only the field at the beginning and the end fiber
% sim.scalar = true; Use scalar propagation
% sim.adaptive_dz.threshold = 1e-5; the threshold of the adaptive-step method
% sim.gpu_yes = true; Use GPU
% sim.include_Raman = 1; Consider the Raman
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.pulse_centering = true; Always shift the pulse to the center of the time window
% sim.gpuDevice.Index = 1; Use the GPU device 1
% sim.progress_bar = true; Show the progress bar
% sim.progress_bar_name = ''; Empty name for the progress bar
% sim.cuda_dir_path = 'UPPE3D/cuda'; Where the cuda files are
% ......
%
% Please check this function for details.
[fiber,sim] = load_default_UPPE3D_propagate([],sim); % load default parameters

fiber.L0 = 1;
num_save = 50;
sim.save_period = fiber.L0/num_save;

%% Initial condition
MFD0 = 500e-6; % m
spatial_window = 10e-3; % m
tfwhm = 1; % ps
time_window = 10; % ps
energy = 1e-3; % nJ
Nt = 2^5; % the number of time points
Nx = 2^7; % the number of spatial points
initial_condition = build_3Dgaussian(MFD0, spatial_window, tfwhm, time_window, energy, Nt, Nx);

% Show initial real space
figure;
pcolor(abs(squeeze(initial_condition.field(Nt/2,:,:))).^2); colormap(jet);colorbar;
shading interp;colormap(jet);colorbar;
title('initial real space');
% Show initial k space
figure;
pcolor(abs(fftshift(ifft(ifft(squeeze(initial_condition.field(Nt/2,:,:)),[],1),[],2))).^2); colormap(jet);colorbar;
shading interp;colormap(jet);colorbar;
title('initial k space');

fiber.n = ones(1,Nx,Nx); % air
fiber.n2 = 0; % no nonlinearity

%% Setup general parameters
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

%% Propagate
% free space before the lens
prop_output1 = UPPE3D_propagate(fiber,initial_condition,sim);
% lens
focal_length = 0.5; % m
radius_of_curvature = 0.25;
Ef_out = add_spherical_lens_phase(ifft(prop_output1.field,[],1),prop_output1.dx,prop_output1.dy,fftshift(lambda,1)/1e9,radius_of_curvature);
% free space after the lens
initial_condition2 = prop_output1; initial_condition2.field = fft(Ef_out(:,:,:,end),[],1);
prop_output2 = UPPE3D_propagate(fiber,initial_condition2,sim);

%% Results
MFD1 = squeeze(calcMFD(squeeze(prop_output1.field(Nt/2,:,:,:)),spatial_window))*1e3;
MFD2 = squeeze(calcMFD(squeeze(prop_output2.field(Nt/2,:,:,:)),spatial_window))*1e3;

z = [prop_output1.z; prop_output2.z(2:end)+prop_output1.z(end)];
MFD = [MFD1;MFD2(2:end)];

%% Theoretical Gaussian propagation
w0 = MFD0/2;
zR0 = pi*w0^2/sim.lambda0; % initial Raylength length
MFD1_theory = MFD0*sqrt(1+(squeeze(prop_output1.z)/zR0).^2)*1e3; % mm

lens_ABCD = @(q,f) q./(1-q/f);
tran_ABCD = @(q,l) q+l;
q1 = fiber.L0 + 1i*zR0;
q1_lens = lens_ABCD(q1,focal_length);
q2 = tran_ABCD(q1_lens,prop_output2.z);
w0_2 = sqrt(imag(q2)/pi*sim.lambda0); % Raylength length after the lens
MFD0_2 = w0_2*2; % MFD0 at the beam waisst after the lens
MFD2_theory = squeeze(MFD0_2).*sqrt(1+(squeeze(real(q2))./imag(q2)).^2)*1e3; % mm

MFD_theory = [MFD1_theory; MFD2_theory(2:end)];

%% Plot
% Show final real space
figure;
pcolor(abs(squeeze(prop_output2.field(Nt/2,:,:,end))).^2); colormap(jet);colorbar;
shading interp;colormap(jet);colorbar;
title('final real space');
% Show final k space
figure;
pcolor(abs(fftshift(ifft(ifft(squeeze(prop_output2.field(Nt/2,:,:,end)),[],1),[],2))).^2); colormap(jet);colorbar;
shading interp;colormap(jet);colorbar;
title('final k space');

% Plot MFD
figure;
plot(z,[MFD,MFD_theory],'linewidth',2);
hold on;
plot(prop_output1.z(end)*[1,1],[0,2],'linewidth',2,'LineStyle','--','Color','k');
hold off
ylim([0,3]);
xlabel('Propagation distance (m)');
ylabel('MFD (mm)');
l = legend('Simulated','Calculated'); set(l,'location','northwest');
set(gca,'fontsize',20);