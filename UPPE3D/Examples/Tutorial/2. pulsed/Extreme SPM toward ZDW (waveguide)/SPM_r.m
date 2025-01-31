% This code demonstrates the SPM toward the zero-dispersion wavelength 
% around 1300 nm in silica.
%
% This script employs the radially-symmetric scheme of the UPPE code, 
% rather than a full x-y dimension.

close all; clearvars;

addpath('../../../../UPPE3D algorithm/','../../../../user_helpers/');

% Add the path for the fiber profile
addpath('../../../../Fibers/Single-mode-fiber profile for examples/');

lambda0 = 1030e-9; % m

%% Setup fiber parameters
sim.lambda0 = lambda0; % the center wavelength
%sim.gpu_yes = false;
sim.include_Raman = false;

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
num_save = 10;%100;
sim.save_period = fiber.L0/num_save;

%fiber.n2 = 0;

%% Information for the Hankel transform
Nr = 2^10; % the number of radial sampling points
r_max = 25e-6; % the maximum radius; half of the spatial window
kr_max = 20e5; % the maximum kr vector

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

%% Setup general parameters
Nt = 2^8; % the number of time points
time_window = 0.8; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
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
tfwhm = 0.03; % ps
total_energy = 10; % nJ
initial_condition = build_MMgaussian(tfwhm, time_window, total_energy, 1, Nt);
initial_condition.field = initial_condition.fields.*phi; initial_condition = rmfield(initial_condition,'fields');
initial_condition.r = r;

%% Show initial spaces
% Show initial real space
figure;
plot(r*1e6,abs(squeeze(initial_condition.field(floor(Nt/2)+1,:))).^2,'linewidth',2,'Color','b');
xlabel('r (\mum)');
xlim([0,10]);
set(gca,'fontsize',20);
title('Initial real space');
% Plot the 2D field with pcolor
% However, the Hankel transform doesn't sample at the origin r=0, so we
% need to find it first. This can be done with Hankel_f_at_0().
A0 = Hankel_f_at_0(initial_condition.field(floor(Nt/2)+1,:,end),l0);
radialPcolor([0,r]*1e6,cat(2,abs(A0).^2,abs(squeeze(initial_condition.field(floor(Nt/2)+1,:,end))).^2));
xlabel('x (\mum)');
ylabel('y (\mum)');
xlim([-10,10]);
ylim([-10,10]);
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

output_field = zeros(Nt,1,num_save+1);
for i = 1:num_save+1
    output_field(:,:,i) = trapz(r,prop_output.field(:,:,i).*conj(phi).*r,2);
end

%% Results
MFD = calcMFD_r(squeeze(prop_output.field(floor(Nt/2)+1,:,:)),r)'*1e6; % um
energy = squeeze(sum(2*pi*trapz(r,abs(prop_output.field).^2.*r,2),1)*dt/1e3); % nJ

%% Plot
% Time
figure;
h = plot(t,abs(output_field(:,:,end)).^2);
xlim([-0.2,0.2]);
xlabel('t');
ylabel('Power');
title('Field');
set(h,'linewidth',2);
set(gca,'fontsize',14);

% Spectrum
figure;
h = plot(f-sim.f0,abs(fftshift(ifft(output_field(:,:,end)),1)).^2);
xlim([-100,100]);
xlabel('\nu-\nu_0');
ylabel('PSD');
title('Spectrum');
set(h,'linewidth',2);
set(gca,'fontsize',14);

% Comparison of time
figure;
[x,y] = meshgrid(t,prop_output.z);
pcolor(x,y,permute(abs(output_field(:,1,:)).^2,[3 1 2]));
shading interp; colormap(jet);
xlim([-0.2,0.2]);
xlabel('t');
ylabel('z');
title('Field during propagation');
set(gca,'fontsize',14);

% Comparison of spectra
figure;
[x_f,y_z] = meshgrid(f-sim.f0,prop_output.z(2:end));
pcolor(x_f,y_z,permute(abs(fftshift(ifft(output_field(:,1,2:end)),1)).^2,[3 1 2]));
xlim([-100,100]);
shading interp; colormap(jet);
xlabel('\nu-\nu_0');
ylabel('z');
title('Spectrum during propagation');
set(gca,'fontsize',14);

% MFD evolution
figure;
plot(prop_output.z,MFD,'linewidth',2,'Color','b');
xlabel('Propagation distance (m)');
ylabel('MFD (\mum)');
set(gca,'fontsize',20);

% Energy
figure;
plot(prop_output.z,energy,'linewidth',2,'Color','b');
xlabel('Propagation distance (m)');
ylabel('Energy (nJ)');
set(gca,'fontsize',20);

% Show final real space
figure;
plot(r*1e6,abs(squeeze(prop_output.field(floor(Nt/2)+1,:,end))).^2,'linewidth',2,'Color','b');
xlabel('r (\mum)');
xlim([0,10]);
set(gca,'fontsize',20);
title('Final real space');
% Plot the 2D field with pcolor
% However, the Hankel transform doesn't sample at the origin r=0, so we
% need to find it first. This can be done with Hankel_f_at_0().
A0 = Hankel_f_at_0(prop_output.field(floor(Nt/2)+1,:,end),l0);
radialPcolor([0,r]*1e6,cat(2,abs(A0).^2,abs(squeeze(prop_output.field(floor(Nt/2)+1,:,end))).^2));
xlabel('x (\mum)');
ylabel('y (\mum)');
xlim([-10,10]);
ylim([-10,10]);
set(gca,'fontsize',20);
daspect([1 1 1]); % make aspect ratio = 1
title('Final real space');
% Show final k space
A_H = 2*pi*FHATHA(squeeze(prop_output.field(floor(Nt/2)+1,:,end)),...
                  r_max,...
                  r,kr,...
                  dr,dkr,...
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

figure;
plot(r*1e6,abs(squeeze(prop_output.field(floor(Nt/2)+1,:,:))).^2,'linewidth',2);
xlabel('r (\mum)');
xlim([0,10]);
set(gca,'fontsize',20);

%% Save the data
%save('SPM3D.mat');