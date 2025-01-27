% This code simulates self-focusing of a Gaussian pulse in a bulk silica.
%
% This script employs the radially-symmetric scheme of the UPPE code, 
% rather than a full x-y dimension.

close all; clearvars;

addpath('../../../../UPPE3D algorithm/','../../../../user_helpers/');

%% Setup fiber parameters
sim.lambda0 = 1030e-9; % the center wavelength
sim.gpuDevice.Index = 1;

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

fiber.L0 = 1e-2;
num_save = 100;
sim.save_period = fiber.L0/num_save;

%% Information for the Hankel transform
Nr = 2^10; % the number of radial sampling points
r_max = 5e-4; % the maximum radius; half of the spatial window
kr_max = 6e5; % the maximum kr vector

[r,kr,...
 l0,exp_prefactor,...
 Q] = Hankel_info(Nr,r_max,kr_max);

% Arrange required Hankel information into "sim" for radially-symmetric
% UPPE to use later.
sim.Hankel = struct('r',r,'kr',kr,'l0',l0,'exp_prefactor',exp_prefactor,'Q',Q);

%% Initial condition
MFD0 = 100e-6; % m
tfwhm = 1; % ps
time_window = 10; % ps
energy = 8e3; % nJ
Nt = 2^8; % the number of time points
initial_condition = build_3Dgaussian_r(MFD0, tfwhm, time_window, energy, Nt, r);

%% Setup general parameters
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

%% Fiber waveguide structure
% Sellmeier coefficients
[a,b] = Sellmeier_coefficients(fiber.material);
% Calculate the index difference using the Sellmeier equation to generate n(lambda)
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));

fiber.n = n_from_Sellmeier(lambda/1e3).*ones(1,Nr);

%% Show initial spaces
% Show initial real space
figure;
plot(r*1e6,abs(squeeze(initial_condition.field(ceil(Nt/2),:))).^2,'linewidth',2,'Color','b');
xlabel('r (\mum)');
xlim([0,100]);
set(gca,'fontsize',20);
title('Initial real space');
% Plot the 2D field with pcolor
% However, the Hankel transform doesn't sample at the origin r=0, so we
% need to find it first. This can be done with Hankel_f_at_0().
A0 = Hankel_f_at_0(initial_condition.field(ceil(Nt/2),:,end),l0);
radialPcolor([0,r]*1e6,cat(2,abs(A0).^2,abs(squeeze(initial_condition.field(ceil(Nt/2),:,end))).^2));
xlabel('x (\mum)');
ylabel('y (\mum)');
xlim([-100,100]);
ylim([-100,100]);
set(gca,'fontsize',20);
daspect([1 1 1]); % make aspect ratio = 1
title('Initial real space');

A0_H = 2*pi*FHATHA(squeeze(initial_condition.field(ceil(Nt/2),:)),...
                   r_max,kr,...
                   l0,exp_prefactor,...
                   Q);

% Show initial k space
figure;
plot(kr/1e6,abs(A0_H).^2,'linewidth',2,'Color','r');
xlabel('k_r (2\pi/\mum)');
set(gca,'fontsize',20);
title('Initial k space');
% Plot the 2D field with pcolor
A0_H0 = Hankel_f_at_0(A0_H,l0);
radialPcolor([0,kr]/1e6,cat(2,abs(A0_H0).^2,abs(A0_H).^2));
xlabel('k_x (2\pi/\mum)');
ylabel('k_y (2\pi/\mum)');
set(gca,'fontsize',20);
daspect([1 1 1]); % make aspect ratio = 1
title('Initial k space');

%% Propagate
prop_output = UPPE3D_propagate(fiber,initial_condition,sim);

%% Results
MFD = calcMFD_r(squeeze(prop_output.field(Nt/2,:,:)),r)'*1e6; % um
energy3D = squeeze(sum(abs(prop_output.field).^2,[1,2]));

%% Plot
% Show final real space
figure;
plot(r*1e6,abs(squeeze(prop_output.field(ceil(Nt/2),:,end))).^2,'linewidth',2,'Color','b');
xlabel('r (\mum)');
xlim([0,100]);
set(gca,'fontsize',20);
title('Final real space');
% Plot the 2D field with pcolor
% However, the Hankel transform doesn't sample at the origin r=0, so we
% need to find it first. This can be done with Hankel_f_at_0().
A0 = Hankel_f_at_0(prop_output.field(ceil(Nt/2),:,end),l0);
radialPcolor([0,r]*1e6,cat(2,abs(A0).^2,abs(squeeze(prop_output.field(ceil(Nt/2),:,end))).^2));
xlabel('x (\mum)');
ylabel('y (\mum)');
xlim([-100,100]);
ylim([-100,100]);
set(gca,'fontsize',20);
daspect([1 1 1]); % make aspect ratio = 1
title('Final real space');
% Show final k space
A_H = 2*pi*FHATHA(squeeze(prop_output.field(ceil(Nt/2),:,end)),...
                  r_max,kr,...
                  l0,exp_prefactor,...
                  Q);
figure;
plot(kr/1e6,abs(A_H).^2,'linewidth',2,'Color','r');
xlabel('k_r (2\pi/\mum)');
set(gca,'fontsize',20);
title('Final k space');
% Plot the 2D field with pcolor
A_H0 = Hankel_f_at_0(A_H,l0);
radialPcolor([0,kr]/1e6,cat(2,abs(A_H0).^2,abs(A_H).^2));
xlabel('k_x (2\pi/\mum)');
ylabel('k_y (2\pi/\mum)');
set(gca,'fontsize',20);
daspect([1 1 1]); % make aspect ratio = 1
title('Final k space');

% Plot MFD
figure;
plot(prop_output.z*1e2,MFD,'linewidth',2,'Color','k');
xlabel('Propagation distance (cm)');
ylabel('MFD (\mum)');
set(gca,'fontsize',20);

% Movie
Frame = animator_r(prop_output.field,l0,....
                   prop_output.z,MFD,...
                   Nt,dt,r,lambda,...
                   fiber.L0);
implay(Frame(:),20);
exportVideo = VideoWriter('self-focusing');
exportVideo.FrameRate = 20;
open(exportVideo);
writeVideo(exportVideo, Frame(:));
close(exportVideo);

save('self_focusing_r.mat');