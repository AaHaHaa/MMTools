% This code simulates a Gaussian beam passing through a thick lens (in air).
% The simulated MFD is compared to the theoretical values.
%
% This script uses the radially-symmetric scheme of the UPPE code, rather
% than a full x-y dimension.

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

%% Information for the Hankel transform
Nr = 2^10; % the number of radial sampling points
r_max = 5e-3; % the maximum radius; half of the spatial window
kr_max = 8e4; % the maximum kr vector

[r,kr,...
 l0,exp_prefactor,...
 Q] = Hankel_info(Nr,r_max,kr_max);

% Arrange required Hankel information into "sim" for radially-symmetric
% UPPE to use later.
sim.Hankel = struct('r',r,'kr',kr,'l0',l0,'exp_prefactor',exp_prefactor,'Q',Q);

%% Initial condition
MFD0 = 500e-6; % m
tfwhm = 1; % ps
time_window = 10; % ps
energy = 1e-3; % nJ
Nt = 2^5; % the number of time points
initial_condition = build_3Dgaussian_r(MFD0, tfwhm, time_window, energy, Nt, r);

% Show initial real space
figure;
plot(r,abs(squeeze(initial_condition.field(ceil(Nt/2),:))).^2);
xlabel('r (m)');
title('initial real space');

A0_H = 2*pi*FHATHA(squeeze(initial_condition.field(ceil(Nt/2),:)),...
                   r_max,kr,...
                   l0,exp_prefactor,...
                   Q);

% Show initial k space
figure;
plot(kr,abs(A0_H).^2);
xlabel('k_r (2\pi/m)');
title('initial k space');

fiber.n = ones(1,Nr); % air
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
Ef_out = add_spherical_lens_phase_r(ifft(prop_output1.field,[],1),prop_output1.r,fftshift(lambda,1)/1e9,radius_of_curvature);
% free space after the lens
initial_condition2 = prop_output1; initial_condition2.field = fft(Ef_out(:,:,:,end),[],1);
prop_output2 = UPPE3D_propagate(fiber,initial_condition2,sim);

%% Results
MFD1 = calcMFD_r(squeeze(prop_output1.field(ceil(Nt/2),:,:,:)),r)'*1e3;
MFD2 = calcMFD_r(squeeze(prop_output2.field(ceil(Nt/2),:,:,:)),r)'*1e3;

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
plot(r,abs(squeeze(prop_output2.field(ceil(Nt/2),:,end))).^2);
xlabel('r (m)');
title('final real space');

A_H = 2*pi*FHATHA(squeeze(prop_output2.field(ceil(Nt/2),:,end)),...
                  r_max,kr,...
                  l0,exp_prefactor,...
                  Q);

% Show final k space
figure;
plot(kr,abs(A_H).^2);
xlabel('k_r (2\pi/m)');
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