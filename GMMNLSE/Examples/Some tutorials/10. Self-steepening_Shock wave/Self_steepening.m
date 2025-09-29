% This code demonstrates the self-steepening of a pulse.
%
% Self-steepening results from the intensity dependence of the group
% velocity, leading to asymmetries in both the temporal and spectral
% profiles.
% The nonlinear phase from SPM is proportional to the |A(t)|^2, which
% in general exhibits a second-order term that corresponds to the
% intensity-dependent group velocity of the pulse.
%
% Please see Chapter 4 of Nonlinear Fiber Optics, Agrawal (5th ed.) for
% details.
%
% This code uses adaptive-step RK4IP for the passive-fiber propagation.

close all; clearvars;

addpath('../../../GMMNLSE algorithm/','../../../user_helpers/');

%% Setup fiber parameters
sim.lambda0 = 1030e-9; % the center wavelength
sim.include_Raman = false; % no Raman
%sim.gpu_yes = false; % don't use GPU

% Load default parameters like 
%
% load fiber.betas and fiber.SR based on your multimode folder above
% sim.f0 = 3e5/sim.lambda0; THz
% sim.save_period = 0; Save only the fields at the beginning and the end fiber
% sim.ellipticity = 0; Linear polarization
% sim.scalar = true; Use scalar propagation
% sim.adaptive_dz.threshold = 1e-5; the threshold of the adaptive-step method
% sim.gpu_yes = true; Use GPU
% sim.include_Raman = true; Consider Raman
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.pulse_centering = true; Always shift the pulse to the center of the time window
% sim.gpuDevice.Index = 1; Use the GPU device 1
% sim.progress_bar = true; Show the progress bar
% sim.progress_bar_name = ''; Empty name for the progress bar
% sim.cuda_dir_path = 'GMMNLSE/cuda'; Where the cuda files are
% ......
%
% Please check this function for details.
[fiber,sim] = load_default_GMMNLSE_propagate([],sim); % load default parameters

num_save = 100;
fiber.L0 = 10;
fiber.betas = zeros(size(fiber.betas)); % no dispersion
sim.save_period = fiber.L0/num_save;

%% Setup general parameters
Nt = 2^13; % the number of time points
time_window = 10; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

%% Initial condition
tfwhm = 1; % ps
total_energy = 5; % nJ
initial_pulse = build_MMgaussian(tfwhm, time_window, total_energy, 1, Nt);

%% Propagate
prop_output = GMMNLSE_propagate(fiber,initial_pulse,sim);

%% Plot
% Time
figure;
plot(t,abs(prop_output.fields(:,:,1)).^2,'linewidth',2);
hold on;
plot(t,abs(prop_output.fields(:,:,end)).^2,'linewidth',2);
hold off;
legend('Initial','Final');
xlim([-3,3]);
xlabel('t');
ylabel('Power');
title('Field');
set(gca,'fontsize',20);

% Spectrum
figure;
plot(f-sim.f0,abs(fftshift(ifft(prop_output.fields(:,:,end)),1)).^2,'linewidth',2);
xlim([-100,100]);
xlabel('\nu-\nu_0');
ylabel('PSD');
title('Spectrum');
set(gca,'fontsize',20);

% Comparison of time
figure;
[x,y] = meshgrid(t,prop_output.z);
pcolor(x,y,permute(abs(prop_output.fields(:,1,:)).^2,[3 1 2]));
shading interp; colormap(jet);
xlim([-3,3]);
xlabel('t');
ylabel('z');
title('Field during propagation');
set(gca,'fontsize',20);

% Comparison of spectra
figure;
[x,y] = meshgrid(f-sim.f0,prop_output.z(2:end));
pcolor(x,y,permute(abs(fftshift(ifft(prop_output.fields(:,1,2:end)),1)).^2,[3 1 2]));
xlim([-100,100]);
shading interp; colormap(jet);
xlabel('\nu-\nu_0');
ylabel('z/LD');
title('Spectrum during propagation');
set(gca,'fontsize',20);

pause(1);

% High-resolution spectrogram
log_yes = false; % Use "log_yes = true" to see, under log scale, how spectral interference generates temporal fringes.
save_point = size(prop_output.fields,3);
Frame(save_point) = struct('cdata',[],'colormap',[]);
for i = 1:save_point
    [~,~,~,figs,ax] = calc_spectrogram(t,f,prop_output.fields(:,1,i),true,[-3,3],[750,1250],400,400,true,true,log_yes);
    set(figs,'Color',[1,1,1]);

    Frame(i) = getframe(figs);
    close(figs);
end
% Movie
implay(Frame,20);

%% Save the data
save('Self_steepening.mat');