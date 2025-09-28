% This code demonstrates the self-phase modulation of a highly-chirped 
% pulse (to 100 ps).
%
% This code uses adaptive-step RK4IP for the passive-fiber propagation with
% newly-proposed narrowband transformation based on scaled Fourier transform.
%
% Narrowband transformation is applied to speed up the simulation, which is
% controlled by sim.cs.cs, the "scaled factor." It is a method to simulate
% pulse propagation for pulses with fewer sampling points below the
% numerical Nyquist limit. It only applies to pulses with smooth profiles
% in both spectral and temporal domains, such as a highly-chirped pulse in
% a CPA system whose temporal profile is a smooth long pulse with a smooth
% dominantly-quadratic phase, along with its smooth spectral profile.
%
% For scaled Fourier transform and its narrowband transformation, please
% find details in
% Chen and Wise, "Field-based treatment of transient gain in short-pulse
% optical amplifiers," Optica 12(6), 879-889 (2025).

close all; clearvars;

addpath('../../GMMNLSE algorithm/','../../user_helpers/');

%% Setup fiber parameters
sim.lambda0 = 1030e-9; % the center wavelength
sim.include_Raman = false; % no Raman
sim.gpu_yes = false; % don't use GPU

% Load default parameters like 
%
% load fiber.betas and fiber.SR based on your multimode folder above
% sim.f0 = 3e5/sim.lambda0; THz
% sim.save_period = 0; Save only the fields at the beginning and the end fiber
% sim.ellipticity = 0; Linear polarization
% sim.scalar = true; Use scalar propagation
% sim.adaptive_dz.model = 1; Use adaptive-step method
% sim.adaptive_dz.threshold = 1e-5; the threshold of the adaptive-step method
% sim.gpu_yes = true; Use GPU
% sim.step_method = 'RK4IP'; Use RK4IP algorithm for pulse propagation
% sim.include_Raman = true; Consider Raman
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.rmc.model = 0; Don't include random mode coupling
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
fiber.L0 = 1;
fiber.betas = zeros(size(fiber.betas)); % no dispersion
sim.save_period = fiber.L0/num_save;

%% Setup general parameters
Nt = 2^17; % the number of time points
time_window = 1500; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

%% Initial condition
tfwhm = 0.1; % ps
total_energy = 10; % nJ
initial_pulse = build_MMgaussian(tfwhm, time_window, total_energy, 1, Nt);

% Chirp the pulse to 100 ps
chirped_tfwhm = 100; % ps
func = calc_chirp;
omega = ifftshift(2*pi*f,1); % 2*pi*THz
[~,chirped_pulse] = func.General( chirped_tfwhm,omega,ifft(initial_pulse.fields),1 );
initial_pulse.fields = chirped_pulse;

%% Apply the narrowband transformation
scaledFT_func = narrowband_scaledFT();

cs = 4;

transformed_At = scaledFT_func.convert(initial_pulse.fields(:,:,end),cs);
transformed_t = linspace(-floor(Nt/cs/2), floor((Nt/cs-1)/2), Nt/cs)*(dt*cs);

% See the transformed field, which should have a narrower spectrum but the
% same temporal profile.
figure;
plot(transformed_t,abs(transformed_At).^2,'linewidth',2,'Color','b');
xlabel('Time (ps)');
ylabel('Power (W)');
set(gca,'fontsize',20);

hold on;
plot(t,abs(initial_pulse.fields(:,:,end)).^2,'linewidth',2,'Color','r');
hold off;
legend('Scaled-FT','Original');
title('Before and after transformation','FontSize',12);

transformed_f = sim.f0+(-Nt/cs/2:Nt/cs/2-1)'/(Nt*dt); % THz

abs2_transformed_Af = abs(fftshift(ifft(transformed_At,[],1),1)).^2;
norm_factor = max(abs2_transformed_Af);

figure;
plot(transformed_f,abs2_transformed_Af/norm_factor,'linewidth',2,'Color','b');
xlabel('Frequency (THz)');
ylabel('PSD');
set(gca,'fontsize',20);

hold on;
plot(f,abs(fftshift(ifft(initial_pulse.fields(:,:,end),[],1),1)).^2/norm_factor,'linewidth',2,'Color','r');
hold off;
legend('Scaled-FT','Original');
title('Before and after transformation','FontSize',12);

% ===========================================
% See Sec. 3 in the supplement of our Optica paper, reference above, for
% details regarding how to apply the scaled factor to physical quantities
% to obtain with the correct (physical-useful) units.
scaled_transformed_f = sim.f0+(-Nt/cs/2:Nt/cs/2-1)'/(Nt*dt)*cs; % THz

figure;
plot(scaled_transformed_f,abs2_transformed_Af/norm_factor,'linewidth',2,'Color','b');
xlabel('Frequency (THz)');
ylabel('PSD');
set(gca,'fontsize',20);

hold on;
plot(f,abs(fftshift(ifft(initial_pulse.fields(:,:,end),[],1),1)).^2/norm_factor*cs,'linewidth',2,'Color','r');
hold off;
legend('Scaled-FT','Original');
title({'Before and after transformation','(Plot w.r.t. original f in x-axis + normalization)'},'FontSize',12);
% ===========================================

% -------------------------------------------------------------------------
% Try to convert the transformed field back to the original scale
At = scaledFT_func.recover(transformed_At,cs);

figure;
plot(transformed_t,abs(transformed_At).^2,'linewidth',2,'Color','b');
xlabel('Time (ps)');
ylabel('Power (W)');
set(gca,'fontsize',20);

hold on;
plot(t,abs(At).^2,'linewidth',2,'Color','r');
hold off;
legend('Scaled-FT','Recovered');
title({'After transformation and then recovered','(Recovered one should be the same as that without transformation)'},'FontSize',12);

figure;
plot(transformed_f,abs(fftshift(ifft(transformed_At,[],1),1)).^2/norm_factor,'linewidth',2,'Color','b');
xlabel('Frequency (THz)');
ylabel('PSD');
set(gca,'fontsize',20);

hold on;
plot(f,abs(fftshift(ifft(At,[],1),1)).^2/norm_factor,'linewidth',2,'Color','r');
hold off;
legend('Scaled-FT','Recovered');
title({'After transformation and then recovered','(Recovered one should be the same as that without transformation)'},'FontSize',12);

% Compare with the original field
figure;
plot(t,abs(initial_pulse.fields(:,:,end)).^2,'linewidth',2,'Color','b');
hold on;
plot(t,abs(At).^2,'linewidth',2,'Color','r');
hold off;
legend('Original','Recovered');
xlabel('Time (ps)');
ylabel('Power (W)');
set(gca,'fontsize',20);
title('Original vs. Recovered from narrowband transformation','FontSize',12);

%% Propagate
% Propagation with the original field
prop_output = GMMNLSE_propagate(fiber,initial_pulse,sim);

% Propagation with the narrowband-transformed field
sim.cs.cs = cs;
transformed_prop_output = GMMNLSE_propagate(fiber,initial_pulse,sim);

%% Improvement
fprintf('Improvement of narrowband transformation: %2.1f\n',prop_output.seconds/transformed_prop_output.seconds);

%% Plot
% Time
figure;
plot(t,abs(prop_output.fields(:,:,end)).^2,'linewidth',2);
hold on;
plot(t,abs(transformed_prop_output.fields(:,:,end)).^2,'linewidth',2);
hold off;
xlim([-200,200]);
xlabel('t');
ylabel('Power');
title('Field');
set(gca,'fontsize',14);

% Spectrum
figure;
plot(f-sim.f0,abs(fftshift(ifft(prop_output.fields(:,:,end)),1)).^2,'linewidth',2);
hold on;
plot(f-sim.f0,abs(fftshift(ifft(transformed_prop_output.fields(:,:,end)),1)).^2,'linewidth',2);
hold off;
xlim([-20,20]);
xlabel('\nu-\nu_0');
ylabel('PSD');
title('Spectrum');
set(gca,'fontsize',14);

analyze_field( t,f,prop_output.fields(:,:,end),'Treacy-t',pi/6,1e-6 );
analyze_field( t,f,transformed_prop_output.fields(:,:,end),'Treacy-t',pi/6,1e-6 );

%% Save the data
save('scaledFT_SPM.mat');