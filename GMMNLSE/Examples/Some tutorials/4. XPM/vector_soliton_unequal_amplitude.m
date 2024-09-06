% This code shows the soliton trapping for vector solitons of unequal amplitude.
% The example is from Ch.6.5.2 of Nonlinear Fiber Optics (5ed), Agrawal.
% Check its Fig.6.15.
%
% For some details about the equations of delta_beta0, delta_beta1, etc,
% please refer to "Stability of solitons in birefringent optical fibers" I
% and II by Curtis R. Menyuk.

close all; clearvars;

addpath('../../../user_helpers/','../../../GMMNLSE algorithm/');

c = 2.99792458e-4; % m/ps

theta = 30; % deg

lambda0 = 1550e-9; % m
tfwhm = 5; % ps

T0 = tfwhm/(2*asech(1/sqrt(2)));    % ps; 2*asech(1/sqrt(2))=1.7627

D = 14; % ps/nm/km; dispersion
beta2 = -(D*1e6)*lambda0^2/(2*pi*c);

% birefringence
delta = 0.2;
delta_beta1 = delta/T0*2*abs(beta2);
delta_beta0 = delta_beta1*(2*pi*c/lambda0);

%% Setup fiber parameters
sim.lambda0 = lambda0;
sim.scalar = false; % polarized fields
sim.include_Raman = false; % no Raman
fiber.betas = [0     delta_beta0;...
               0     delta_beta1;...
               beta2 beta2];
sim.betas = [delta_beta0/2; delta_beta1/2]; % set the betas for the SVEA and the moving frame of the simulation
LD = T0^2/abs(beta2); % dispersion length
num_save = 20;
fiber.L0 = 10*LD; % propagate for 10*LD
sim.save_period = fiber.L0/num_save;
sim.gpu_yes = false;

% Load default parameters like 
%
% loading fiber.betas and fiber.SR based on your multimode folder above
% sim.include_Raman = true; Consider Raman
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.gpu_yes = true; Use GPU (default to true)
% ......
%
% Please check this function for details.
[fiber,sim] = load_default_GMMNLSE_propagate(fiber,sim); % load default parameters

%% Setup general parameters
Nt = 2^12; % the number of time points
time_window = 500; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

%% Initial condition
theta = theta*pi/180;

Aeff = (fiber.MFD/2*1e-6)^2*pi; % mode area; m^2

% solitons
initial_pulse = build_MMsoliton(tfwhm, fiber.betas(3), 1/Aeff, lambda0, time_window, 1, Nt);
initial_pulse.fields = initial_pulse.fields.*[sin(theta) cos(theta)];

%% Propagate
prop_output = GMMNLSE_propagate(fiber,initial_pulse,sim);

%% Plot
% Time
figure;
plot(t/T0,abs(prop_output.fields(:,:,end)).^2,'linewidth',2);
xlim([-4 6]);
xlabel('t/T_0');
ylabel('Intensity');

% Spectrum
figure;
plot((f-sim.f0)*T0,abs(fftshift(ifft(prop_output.fields(:,:,end)),1)).^2,'linewidth',2);
xlim([-0.25 0.25]);
xlabel('(\nu-\nu_0)*T_0');
ylabel('Intensity');

% Comparison of time
figure;
subplot(2,1,1);
[x,y] = meshgrid(t/T0,prop_output.z/LD);
pcolor(x,y,permute(abs(prop_output.fields(:,1,:)).^2,[3 1 2]));
shading interp; colormap jet;
xlim([-10 10]);
xlabel('t/T_0');
ylabel('z/LD');

subplot(2,1,2);
[x,y] = meshgrid(t/T0,prop_output.z/LD);
pcolor(x,y,permute(abs(prop_output.fields(:,2,:)).^2,[3 1 2]));
shading interp; colormap jet;
xlim([-10 10]);
xlabel('t/T_0');
ylabel('z/LD');

% Comparison of spectra
figure;
subplot(2,1,1);
[x,y] = meshgrid((f-sim.f0)*T0,prop_output.z/LD);
pcolor(x,y,permute(abs(fftshift(ifft(prop_output.fields(:,1,:)),1)).^2,[3 1 2]));
shading interp; colormap jet;
xlim([-0.5 0.5]);
xlabel('(\nu-\nu_0)*T_0');
ylabel('z/LD');

subplot(2,1,2);
[x,y] = meshgrid((f-sim.f0)*T0,prop_output.z/LD);
pcolor(x,y,permute(abs(fftshift(ifft(prop_output.fields(:,2,:)),1)).^2,[3 1 2]));
shading interp; colormap jet;
xlim([-0.5 0.5]);
xlabel('(\nu-\nu_0)*T_0');
ylabel('z/LD');

%% Save the data
save('vector_soliton_unequal_amplitude.mat');