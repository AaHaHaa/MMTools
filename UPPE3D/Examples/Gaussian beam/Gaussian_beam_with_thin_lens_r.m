% This code simulates a Gaussian beam passing through a thin lens (in air).
% The simulated MFD is compared to the theoretical values.
%
% This script employs the radially-symmetric scheme of the UPPE code, 
% rather than a full x-y dimension.
%
% In addition, this script simulates with continuous-wave (CW) light, which
% focuses only on the evolution of its spatial profile through the
% diffraction effect.
%
% To run with CW, make Nt = 1.

close all; clearvars;

addpath('../../UPPE3D algorithm/','../../user_helpers/');

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
kr_max = 3e4; % the maximum kr vector

[r,kr,...
 l0,exp_prefactor,...
 Q] = Hankel_info(Nr,r_max,kr_max);

dr = diff(r,1,2);
dkr = diff(kr,1,2);
% Arrange required Hankel information into "sim" for radially-symmetric
% UPPE to use later.
sim.Hankel = struct('r',r,'kr',kr,...
                    'dr',dr,'dkr',dkr,...
                    'l0',l0,'exp_prefactor',exp_prefactor,'Q',Q);

%% Initial condition
MFD0 = 500e-6; % m
tfwhm = 1; % ps
time_window = 10; % ps
energy = 1e-3; % nJ
Nt = 1; % the number of time points
initial_condition = build_3Dgaussian_r(MFD0, tfwhm, time_window, energy, Nt, r);

fiber.n = ones(1,Nr); % air
fiber.n2 = 0; % no nonlinearity

%% Show initial spaces
% Show initial real space
figure;
plot(r*1e3,abs(squeeze(initial_condition.field(floor(Nt/2)+1,:))).^2,'linewidth',2,'Color','b');
xlabel('r (mm)');
set(gca,'fontsize',20);
title('Initial real space');
% Plot the 2D field with pcolor
% However, the Hankel transform doesn't sample at the origin r=0, so we
% need to find it first. This can be done with Hankel_f_at_0().
A0 = Hankel_f_at_0(initial_condition.field(floor(Nt/2)+1,:,end),l0);
radialPcolor([0,r]*1e3,cat(2,abs(A0).^2,abs(squeeze(initial_condition.field(floor(Nt/2)+1,:,end))).^2));
xlabel('r (mm)');
ylabel('r (mm)');
set(gca,'fontsize',20);
daspect([1 1 1]); % make aspect ratio = 1
title('Initial real space');

A0_H = 2*pi*FHATHA(squeeze(initial_condition.field(floor(Nt/2)+1,:)),...
                   r_max,...
                   r,kr,...
                   dr,dkr,...
                   l0,exp_prefactor,...
                   Q);

% Show initial k space
figure;
plot(kr/1e3,abs(A0_H).^2,'linewidth',2,'Color','r');
xlabel('k_r (2\pi/mm)');
set(gca,'fontsize',20);
title('Initial k space');
% Plot the 2D field with pcolor
A0_H0 = Hankel_f_at_0(A0_H,l0);
radialPcolor([0,kr]/1e3,cat(2,abs(A0_H0).^2,abs(A0_H).^2));
xlabel('k_r (2\pi/mm)');
ylabel('k_r (2\pi/mm)');
set(gca,'fontsize',20);
daspect([1 1 1]); % make aspect ratio = 1
title('Initial k space');

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
Ef_out = add_thin_lens_phase_r(ifft(prop_output1.field,[],1),prop_output1.r,fftshift(lambda,1)/1e9,focal_length);
% free space after the lens
initial_condition2 = prop_output1; initial_condition2.field = fft(Ef_out(:,:,:,end),[],1);
prop_output2 = UPPE3D_propagate(fiber,initial_condition2,sim);

%% Results
% mode-field diameter
MFD1 = calcMFD_r(squeeze(sum(prop_output1.field,1)),r)'*1e3; % mm
MFD2 = calcMFD_r(squeeze(sum(prop_output2.field,1)),r)'*1e3; % mm

% power
optical_power1 = squeeze(sum(2*pi*trapz(r,abs(prop_output1.field).^2.*r,2),1)); % W
optical_power2 = squeeze(sum(2*pi*trapz(r,abs(prop_output2.field).^2.*r,2),1)); % W

z = [prop_output1.z; prop_output2.z(2:end)+prop_output1.z(end)]; % m
MFD = [MFD1;MFD2(2:end)]; % mm
optical_power = [optical_power1;optical_power2(2:end)]; % W

%% Theoretical Gaussian propagation
w0 = MFD0/2; % m
zR0 = pi*w0^2/sim.lambda0; % initial Raylength length (m)
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
plot(r*1e3,abs(squeeze(prop_output2.field(floor(Nt/2)+1,:,end))).^2,'linewidth',2,'Color','b');
xlabel('r (mm)');
set(gca,'fontsize',20);
title('Final real space');
% Plot the 2D field with pcolor
% However, the Hankel transform doesn't sample at the origin r=0, so we
% need to find it first. This can be done with Hankel_f_at_0().
A0 = Hankel_f_at_0(prop_output2.field(floor(Nt/2)+1,:,end),l0);
radialPcolor([0,r]*1e3,cat(2,abs(A0).^2,abs(squeeze(prop_output2.field(floor(Nt/2)+1,:,end))).^2));
xlabel('r (mm)');
ylabel('r (mm)');
set(gca,'fontsize',20);
daspect([1 1 1]); % make aspect ratio = 1
title('Final real space');
% Show final k space
A_H = 2*pi*FHATHA(squeeze(prop_output2.field(floor(Nt/2)+1,:,end)),...
                  r_max,...
                  r,kr,...
                  dr,dkr,...
                  l0,exp_prefactor,...
                  Q);
figure;
plot(kr/1e3,abs(A_H).^2,'linewidth',2,'Color','r');
xlabel('k_r (2\pi/mm)');
set(gca,'fontsize',20);
title('Final k space');
% Plot the 2D field with pcolor
A_H0 = Hankel_f_at_0(A_H,l0);
radialPcolor([0,kr]/1e3,cat(2,abs(A_H0).^2,abs(A_H).^2));
xlabel('k_r (2\pi/mm)');
ylabel('k_r (2\pi/mm)');
set(gca,'fontsize',20);
daspect([1 1 1]); % make aspect ratio = 1
title('Final k space');

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

% Power
% Check that whether it's conserved
figure;
plot(z,optical_power,'linewidth',2,'Color','b');
hold on;
plot(prop_output1.z(end)*[1,1],[0.9,1],'linewidth',2,'LineStyle','--','Color','k');
hold off
xlabel('Propagation distance (m)');
ylabel('Power (W)');
set(gca,'fontsize',20);