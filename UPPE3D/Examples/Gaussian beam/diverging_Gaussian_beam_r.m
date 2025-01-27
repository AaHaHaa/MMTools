% This code simulates a diverging Gaussian beam in air. The simulated MFD
% is compared to the theoretical values.
%
% This script employs the radially-symmetric scheme of the UPPE code, 
% rather than a full x-y dimension.

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
num_save = 100;
sim.save_period = fiber.L0/num_save;

%% Information for the Hankel transform
Nr = 2^8; % the number of radial sampling points
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
Nt = 1; % the number of time points
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
prop_output = UPPE3D_propagate(fiber,initial_condition,sim);

%% Results
MFD = calcMFD_r(squeeze(prop_output.field(ceil(Nt/2),:,:)),r)'*1e3;
energy3D = squeeze(sum(abs(prop_output.field).^2,[1,2]));

%% Theoretical Gaussian propagation
w0 = MFD0/2;
zR = pi*w0^2/sim.lambda0; % Raylength length
MFD_theory = MFD0*sqrt(1+(squeeze(prop_output.z)/zR).^2)*1e3; % mm

%% Plot
% Show final real space
figure;
plot(r,abs(squeeze(prop_output.field(ceil(Nt/2),:,end))).^2);
% Show final k space
A_H = 2*pi*FHATHA(squeeze(prop_output.field(ceil(Nt/2),:,end)),...
                  r_max,kr,...
                  l0,exp_prefactor,...
                  Q);
figure;
plot(kr,abs(A_H).^2);
xlabel('k_r (m)');

% Plot MFD
figure;
plot(prop_output.z,[MFD,MFD_theory],'linewidth',2);
xlabel('Propagation distance (m)');
ylabel('MFD (mm)');
l = legend('Simulated','Calculated'); set(l,'location','northwest');
set(gca,'fontsize',20);