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
num_save = 10;
sim.save_period = fiber.L0/num_save;

%% Information for the Hankel transform
Nr = 2^11; % the number of radial sampling points
r_max = 5e-4; % the maximum radius; half of the spatial window
kr_max = 2e5; % the maximum kr vector

[r,kr,...
 dr,dkr,...
 l0,exp_prefactor,n2_prefactor,...
 ifftQ] = Hankel_info(Nr,r_max,kr_max);

% Arrange required Hankel information into "sim" for radially-symmetric
% UPPE to use later.
sim.Hankel = struct('r',r,'kr',kr,...
                    'dr',dr,'dkr',dkr,...
                    'l0',l0,'exp_prefactor',exp_prefactor,'n2_prefactor',n2_prefactor,...
                    'ifftQ',ifftQ);

%% Initial condition
MFD0 = 100e-6; % m
tfwhm = 1; % ps
time_window = 10; % ps
energy = 9e3; % nJ
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

fiber.n = n_from_Sellmeier(lambda/1e3);

%% Show initial spaces
% Show initial real space
figure;
plot(r*1e6,abs(squeeze(initial_condition.field(floor(Nt/2)+1,:))).^2,'linewidth',2,'Color','b');
xlabel('r (\mum)');
xlim([0,100]);
set(gca,'fontsize',20);
title('Initial real space');
% Plot the 2D field with pcolor
radialPcolor(r*1e6,abs(squeeze(initial_condition.field(floor(Nt/2)+1,:,end))).^2);
xlabel('x (\mum)');
ylabel('y (\mum)');
xlim([-100,100]);
ylim([-100,100]);
set(gca,'fontsize',20);
daspect([1 1 1]); % make aspect ratio = 1
title('Initial real space');

A0_H = 2*pi*FHATHA(squeeze(initial_condition.field(floor(Nt/2)+1,:)),...
                   r_max,...
                   r,kr,...
                   dr,dkr,...
                   l0,exp_prefactor,n2_prefactor,...
                   ifftQ);

% Show initial k space
figure;
plot(kr/1e6,abs(A0_H).^2,'linewidth',2,'Color','r');
xlabel('k_r (2\pi/\mum)');
set(gca,'fontsize',20);
title('Initial k space');
% Plot the 2D field with pcolor
radialPcolor(kr/1e6,abs(A0_H).^2);
xlabel('k_x (2\pi/\mum)');
ylabel('k_y (2\pi/\mum)');
set(gca,'fontsize',20);
daspect([1 1 1]); % make aspect ratio = 1
title('Initial k space');

%% Propagate
prop_output = UPPE3D_propagate(fiber,initial_condition,sim);

%% Results
MFD = calcMFD_r(squeeze(prop_output.field(floor(Nt/2)+1,:,:)),r)'*1e6; % um
optical_energy = squeeze(sum(2*pi*trapz(r,abs(prop_output.field).^2.*r,2),1)*dt/1e3); % nJ

%% Plot
% Show final real space
figure;
plot(r*1e6,abs(squeeze(prop_output.field(floor(Nt/2)+1,:,end))).^2,'linewidth',2,'Color','b');
xlabel('r (\mum)');
xlim([0,100]);
set(gca,'fontsize',20);
title('Final real space');
% Plot the 2D field with pcolor
radialPcolor(r*1e6,abs(squeeze(prop_output.field(floor(Nt/2)+1,:,end))).^2);
xlabel('x (\mum)');
ylabel('y (\mum)');
xlim([-100,100]);
ylim([-100,100]);
set(gca,'fontsize',20);
daspect([1 1 1]); % make aspect ratio = 1
title('Final real space');
% Show final k space
A_H = 2*pi*FHATHA(squeeze(prop_output.field(floor(Nt/2)+1,:,end)),...
                  r_max,...
                  r,kr,...
                  dr,dkr,...
                  l0,exp_prefactor,n2_prefactor,...
                  ifftQ);
figure;
plot(kr/1e6,abs(A_H).^2,'linewidth',2,'Color','r');
xlabel('k_r (2\pi/\mum)');
set(gca,'fontsize',20);
title('Final k space');
% Plot the 2D field with pcolor
radialPcolor(kr/1e6,abs(A_H).^2);
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

% Energy
figure;
plot(prop_output.z,optical_energy,'linewidth',2,'Color','b');
xlabel('Propagation distance (m)');
ylabel('Power (nJ)');
set(gca,'fontsize',20);

% Movie
Frame = animator_r(prop_output.field,....
                   prop_output.z,MFD,...
                   Nt,dt,r,lambda,...
                   fiber.L0);
implay(Frame(:),20);
exportVideo = VideoWriter('self-focusing');
exportVideo.FrameRate = 20;
open(exportVideo);
writeVideo(exportVideo, Frame(:));
close(exportVideo);

%save('self_focusing_r.mat');