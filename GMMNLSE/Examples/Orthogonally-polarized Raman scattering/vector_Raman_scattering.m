% This code shows the orthogonally polarized soliton trapping, plus energy 
% transfer, under Raman gain and soliton self-frequency shift (SSFS).
%
% The input pulse includes a pump soliton at 1550 nm in the slow axis and a
% seed pulse at 1650 nm in the fast axis. The pump soliton will transfer
% energy to the seed pulse. They will trap and propagate with the same
% speed despite potential group velocity mismatch, which is compensated by 
% bhe wavelength differernce. Eventually, all the energy will transfer to
% the other polarization with the longer wavelength due to the Raman
% scattering.
%
% This demonstration follows the following paper. I didn't intend to
% perfect duplication the one in the paper, and only focused on the main
% physical phenomenon, so I didn't make all the parameters the same as in
% the paper. By the way, the paper doesn't show what happens with extended
% propagation. Here, I show that all the energy will be transferred to the
% other polarization with the longer wavelength due to the Raman
% scattering.
%
% Reference:
% Nishizawa and Goto, "Trapped pulse generation by femtosecond soliton
% pulse in birefringent optical fibers," 10, 256-261 (2002)

close all; clearvars;

addpath('../../GMMNLSE algorithm/','../../user_helpers/');

lambda0 = 1550e-9; % m
tfwhm = 0.11; % ps

T0 = tfwhm/(2*asech(1/sqrt(2))); % ps; 2*asech(1/sqrt(2))=1.7627

N = 1; % soliton number

%% Setup fiber parameters
sim.lambda0 = lambda0; % the central wavelength
%sim.pulse_centering = false; % comment this out: because the pulse wavelength is shifting, pulse-centering keeps the pulse at the time window with the time delay saved in prop_output.t_delay
%sim.gpu_yes = false; % comment this out: polarized simulation is faster with GPU
sim.scalar = false; % run with polarized modes

% Load default parameters like 
%
% loading fiber.betas and fiber.SR based on your multimode folder above
% sim.include_Raman = true; Consider Raman
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.gpu_yes = true; Use GPU (default to true)
% ......
%
% Please check this function for details.
[fiber,sim] = load_default_GMMNLSE_propagate([],sim); % load default parameters

%% Create birefringence

fiber.betas(3) = -19e-3;

% These lines are after load_default_GMMNLSE_propagate() because it needs 
% to load the default betas at 1550 nm first.

beat_length = 2.7e-3; % 2.7mm at 980nm for PM980
% The betas at the other polarization can be easily calculated by 
% beta_e = n_e*w/c = (n_o + delta_n)*w/c
%                  = (n_o + lambda/beat_length)*w/c
c = 299792458e-12; % m/ps
fiber.betas = [fiber.betas(1), fiber.betas(1) - 980e-9/beat_length *2*pi/sim.lambda0;...
               fiber.betas(2), fiber.betas(2) - 980e-9/beat_length/c;...
               fiber.betas(3), fiber.betas(3);...
               fiber.betas(4), fiber.betas(4);...
               fiber.betas(5), fiber.betas(5)];

LD = T0^2/abs(fiber.betas(3)); % dispersion length
num_save = 100;
fiber.L0 = 500; % m
sim.save_period = fiber.L0/num_save;

%% Setup general parameters
Nt = 2^12; % the number of time points
time_window = 20; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

%% Initial condition
Aeff = 1/fiber.SR;
initial_pulse1 = build_MMsoliton(tfwhm, fiber.betas(3), 1/Aeff, lambda0, time_window, 1, Nt, {'ifft',0}, N);

pulse_lambda0 = 1650e-9;
f_now = c/sim.lambda0*1e-12;
f_pulse = c/pulse_lambda0*1e-12;
freq_shift = f_pulse - f_now;
initial_pulse2 = build_MMsoliton(tfwhm, fiber.betas(3), 1/Aeff, lambda0, time_window, 1, Nt, {'ifft',freq_shift}, N);

input_pulse = initial_pulse1;
input_pulse.fields = [initial_pulse1.fields,initial_pulse2.fields*sqrt(1e-3)];

%% Propagate
prop_output = GMMNLSE_propagate(fiber,input_pulse,sim);

%% Plot
% Time
figure;
plot(t/T0,abs(prop_output.fields(:,:,end)).^2/1e3,'linewidth',2);
xlabel('t/T_0');
ylabel('Power (kW)');
title('Field');
set(gca,'fontsize',20);

% Spectrum
figure;
plot(lambda,abs(fftshift(ifft(prop_output.fields(:,:,end)),1)).^2*c./lambda.^2,'linewidth',2);
xlabel('\lambda (nm)');
ylabel('Intensity (a.u.)');
xlim([1400,2500]);
title('Spectrum');
set(gca,'fontsize',14);

% Comparison of time
figure;
[x,y] = meshgrid(t,prop_output.z);
pcolor(x,y,sum(permute(abs(prop_output.fields).^2,[3 1 2]),3));
shading interp; colormap(jet);
xlabel('t/T_0');
ylabel('z (m)');
title('Pulse evolution');
set(gca,'fontsize',20);

% Comparison of spectra
figure;
[x,y] = meshgrid((f-sim.f0),prop_output.z);
tmp = 10*log10(sum(permute(abs(fftshift(ifft(prop_output.fields),1)).^2,[3 1 2]),3)); tmp = tmp - max(tmp(:));
pcolor(x,y,tmp);
shading interp; colormap(jet); caxis([-20,0]);
xlabel('\Deltaf (THz)');
ylabel('z (m)');
title('Spectral evolution');
set(gca,'fontsize',20);
print(gcf,'Spectral evolution.jpg','-djpeg');

% Energy transfer
figure;
Ex = squeeze(trapz(t,abs(prop_output.fields(:,1,:)).^2,1))/1e3; % nJ
Ey = squeeze(trapz(t,abs(prop_output.fields(:,2,:)).^2,1))/1e3;
plot(prop_output.z,[Ex,Ey],'linewidth',2);
xlabel('z (m)');
ylabel('Energy (nJ)');
set(gca,'fontsize',20);
xlim([0,fiber.L0]);
legend('x','y');
print(gcf,'Energy transfer.jpg','-djpeg');