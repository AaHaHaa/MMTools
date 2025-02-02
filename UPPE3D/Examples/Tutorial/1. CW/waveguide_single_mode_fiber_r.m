% This code simulates the waveguide effect in a single-mode fiber.
%
% This script employs the radially-symmetric scheme of the UPPE code, 
% rather than a full x-y dimension.
%
% In addition, this script simulates with continuous-wave (CW) light, which
% focuses only on the evolution of its spatial profile through the
% waveguide and diffraction effect.
%
% To run with CW, make Nt = 1.
%
% This code loads a mode profile generated from the FEM mode solver.
% However, in transforming it into the field to be computed in 3D-UPPE with
% an approach based on effective mode-field area, small deviation happens.
% This is why in MFD evolution, we can see variations of the mode spatial
% profile, which stabilizes to a slightly-different mode profile.

close all; clearvars;

addpath('../../../UPPE3D algorithm/','../../../user_helpers/');

% Add the path for the fiber profile
addpath('../../../Fibers/Single-mode-fiber profile for examples/');

lambda0 = 1030e-9; % m

%% Setup fiber parameters
sim.lambda0 = lambda0; % the center wavelength
%im.gpu_yes = false;

% Load default parameters like
%
% load fiber.betas and fiber.SR based on your multimode folder above
% sim.f0 = 3e5/sim.lambda0; THz
% sim.save_period = 0; Save only the field at the beginning and the end fiber
% sim.ellipticity = 0; Linear polarization
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

fiber.L0 = 0.01;
num_save = 10;
sim.save_period = fiber.L0/num_save;

%% Information for the Hankel transform
Nr = 2^11; % the number of radial sampling points
r_max = 10e-6; % the maximum radius; half of the spatial window
kr_max = 20e5; % the maximum kr vector

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

%% Setup general parameters
Nt = 1; % the number of time points
time_window = 0.8; % ps; it is meaningless in this code since Nt = 1 (CW case)
dt = time_window/Nt; % ps; it is meaningless in this code since Nt = 1 (CW case)
f = sim.f0+(-floor(Nt/2):floor((Nt-1)/2))'/(Nt*dt); % THz
t = (-floor(Nt/2):floor((Nt-1)/2))'*dt; % ps; it is meaningless in this code since Nt = 1 (CW case)
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

%% Initial condition
fiber.n = calc_index_profile_r(fiber,r,f);

loaded_data = load('mode1wavelength10300.mat','phi','x','epsilon');

% take only the central part
phi = interp1(loaded_data.x*1e-6,loaded_data.phi(:,ceil(length(loaded_data.x)/2)),r,'linear',0);
phi = phi/sqrt(2*pi*trapz(r,abs(phi).^2.*r));

Aeff = (2*pi*trapz(r,abs(phi).^2.*r)).^2./(2*pi*trapz(r,abs(phi).^4.*r)); % effective mode-field area (m^2)
MFR = sqrt(Aeff/pi)*1e6; % um

% input pulse
% Pulse duration and energy itself have no meaning here. Only the peak
% power is retained as the CW power in building the initial profile.
% However, still ensure that pulse duration is around 5-10x smaller than
% the time window for correct generation of the CW profile.
tfwhm = 1; % ps
total_energy = 1e-3; % nJ
initial_condition = build_MMgaussian(tfwhm, time_window, total_energy, 1, Nt);
initial_condition.field = initial_condition.fields.*phi; initial_condition = rmfield(initial_condition,'fields');
initial_condition.r = r;

%% Show initial spaces
% Show initial real space
figure;
plot(r*1e6,abs(squeeze(initial_condition.field(floor(Nt/2)+1,:))).^2,'linewidth',2,'Color','b');
xlabel('r (\mum)');
set(gca,'fontsize',20);
title('Initial real space');
% Plot the 2D field with pcolor
radialPcolor(r*1e6,abs(squeeze(initial_condition.field(floor(Nt/2)+1,:,end))).^2);
xlabel('x (\mum)');
ylabel('y (\mum)');
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
MFD = calcMFD_r(squeeze(prop_output.field(floor(Nt/2)+1,:,:)),r)'*1e3; % mm
optical_power = squeeze(2*pi*trapz(r,abs(prop_output.field).^2.*r,2)); % W

%% Plot
% Show final real space
figure;
plot(r*1e6,abs(squeeze(prop_output.field(floor(Nt/2)+1,:,end))).^2,'linewidth',2,'Color','b');
xlabel('r (\mum)');
set(gca,'fontsize',20);
title('Final real space');
% Plot the 2D field with pcolor
radialPcolor(r*1e6,abs(squeeze(prop_output.field(floor(Nt/2)+1,:,end))).^2);
xlabel('x (\mum)');
ylabel('y (\mum)');
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

% MFD evolution
figure;
plot(prop_output.z,MFD,'linewidth',2,'Color','b');
xlabel('Propagation distance (m)');
ylabel('MFD (mm)');
set(gca,'fontsize',20);

% Power
% Check that whether it's conserved
figure;
plot(prop_output.z,optical_power,'linewidth',2,'Color','b');
xlabel('Propagation distance (m)');
ylabel('Power (W)');
set(gca,'fontsize',20);

% Mode profile
figure;
plot(r*1e6,squeeze(abs(prop_output.field).^2),'linewidth',2);
xlabel('r (\mum)');
xlim([0,10]);
ylabel('Intensity (W/m^2)');

%% Save the data
%save('SPM3D.mat');