% This code simulates a Gaussian beam passing through a Kepler telescope with two convex lenses.
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

num_save = 20;

%% Initial condition
MFD0 = 1e-3; % m
spatial_window = 10e-3; % m
tfwhm = 1; % ps
time_window = 10; % ps
energy = 1e-3; % nJ
Nt = 1; % the number of time points
Nx = 2^10; % the number of spatial points
initial_condition = build_3Dgaussian(MFD0, spatial_window, tfwhm, time_window, energy, Nt, Nx);

% Show initial real space
figure;
pcolor(abs(squeeze(initial_condition.field(ceil(Nt/2),:,:))).^2); colormap(jet);colorbar;
shading interp;colormap(jet);colorbar;
title('initial real space');
% Show initial k space
figure;
pcolor(abs(fftshift(fft(fft(squeeze(initial_condition.field(ceil(Nt/2),:,:)),[],1),[],2))).^2); colormap(jet);colorbar;
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
% free space before the telescope
fiber1 = fiber;
fiber1.L0 = 0.01;
sim1 = sim;
sim1.save_period = fiber1.L0/num_save;
prop_output1 = UPPE3D_propagate(fiber1,initial_condition,sim1);

% -------------------------------------------------------------------------
% Telescope
% -------------------------------------------------------------------------
convex_focal_length1 = 0.9;
convex_focal_length2 = 0.3;

% convex lens 1
convex_radius_of_curvature1 = convex_focal_length1/2;
%Ef_out = add_spherical_lens_phase(ifft(prop_output1.field,[],1),prop_output1.dx,prop_output1.dy,fftshift(lambda,1)/1e9,convex_radius_of_curvature1);
Ef_out = add_thin_lens_phase(ifft(prop_output1.field,[],1),prop_output1.dx,prop_output1.dy,fftshift(lambda,1)/1e9,convex_focal_length1);
% free space after the convex lens 1
fiber2 = fiber;
fiber2.L0 = convex_focal_length1 + convex_focal_length2;
sim2 = sim;
sim2.save_period = fiber2.L0/num_save;
initial_condition2 = prop_output1; initial_condition2.field = fft(Ef_out(:,:,:,end),[],1);
prop_output2 = UPPE3D_propagate(fiber2,initial_condition2,sim2);

% convex lens 2
convex_radius_of_curvature2 = convex_focal_length2/2;
%Ef_out = add_spherical_lens_phase(ifft(prop_output2.field,[],1),prop_output2.dx,prop_output2.dy,fftshift(lambda,1)/1e9,convex_radius_of_curvature2);
Ef_out = add_thin_lens_phase(ifft(prop_output2.field,[],1),prop_output2.dx,prop_output2.dy,fftshift(lambda,1)/1e9,convex_focal_length2);
% free space after the convex lens 2
fiber3 = fiber;
fiber3.L0 = 0.5;
sim3 = sim;
sim3.save_period = fiber3.L0/num_save;
initial_condition3 = prop_output1; initial_condition3.field = fft(Ef_out(:,:,:,end),[],1);
prop_output3 = UPPE3D_propagate(fiber3,initial_condition3,sim3);

%% Results
MFD1 = squeeze(calcMFD(squeeze(sum(prop_output1.field,1)),spatial_window))*1e3;
MFD2 = squeeze(calcMFD(squeeze(sum(prop_output2.field,1)),spatial_window))*1e3;
MFD3 = squeeze(calcMFD(squeeze(sum(prop_output3.field,1)),spatial_window))*1e3;

z = [prop_output1.z; prop_output2.z(2:end)+prop_output1.z(end); prop_output3.z(2:end)+prop_output1.z(end)+prop_output2.z(end)];
MFD = [MFD1;MFD2(2:end);MFD3(2:end)];

%% Theoretical Gaussian propagation
w0 = MFD0/2;
zR0 = pi*w0^2/sim.lambda0; % initial Raylength length
MFD1_theory = MFD0*sqrt(1+(squeeze(prop_output1.z)/zR0).^2)*1e3; % mm

lens_ABCD = @(q,f) q./(1-q/f);
tran_ABCD = @(q,l) q+l;

% convex lens 1
q1 = fiber1.L0 + 1i*zR0;
q1_lens = lens_ABCD(q1,convex_focal_length1);
q2 = tran_ABCD(q1_lens,prop_output2.z);
w0_2 = sqrt(imag(q2)/pi*sim.lambda0); % Raylength length after the lens
MFD0_2 = w0_2*2; % MFD0 at the beam waisst after the lens
MFD2_theory = squeeze(MFD0_2).*sqrt(1+(squeeze(real(q2))./imag(q2)).^2)*1e3; % mm

% convex lens 2
q2_lens = lens_ABCD(q2(end),convex_focal_length2);
q3 = tran_ABCD(q2_lens,prop_output3.z);
w0_3 = sqrt(imag(q3)/pi*sim.lambda0); % Raylength length after the lens
MFD0_3 = w0_3*2; % MFD0 at the beam waisst after the lens
MFD3_theory = squeeze(MFD0_3).*sqrt(1+(squeeze(real(q3))./imag(q3)).^2)*1e3; % mm

MFD_theory = [MFD1_theory; MFD2_theory(2:end); MFD3_theory(2:end)];

%% Plot
% Show final real space
figure;
pcolor(abs(squeeze(prop_output3.field(ceil(Nt/2),:,:,end))).^2); colormap(jet);colorbar;
shading interp;colormap(jet);colorbar;
title('final real space');
% Show final k space
figure;
pcolor(abs(fftshift(fft(fft(squeeze(prop_output3.field(ceil(Nt/2),:,:,end)),[],1),[],2))).^2); colormap(jet);colorbar;
shading interp;colormap(jet);colorbar;
title('final k space');

% Plot MFD
figure;
plot(z,[MFD,MFD_theory],'linewidth',2);
hold on;
plot(prop_output1.z(end)*[1,1],[0,10],'linewidth',2,'LineStyle','--','Color','k');
plot((prop_output1.z(end)+prop_output2.z(end))*[1,1],[0,10],'linewidth',2,'LineStyle','--','Color','k');
hold off
ylim([0,2]);
xlabel('Propagation distance (m)');
ylabel('MFD (mm)');
l = legend('Simulated','Calculated'); set(l,'location','northwest');
set(gca,'fontsize',20);