% This code simulates the propagation of a soliton with soliton number =
% 1.1. It will evolve osciilatorily until it reaches the new fundamental
% soliton with a different pulse duration.
%
% See Ch. 5.2.5, Soliton Stability, in Nonlinear Fiber Optics by Agrawal (5ed)
%
% Raman scattering is turned off to avoid soliton self-frequency shift.

close all; clearvars;

addpath('../../../GMMNLSE algorithm/','../../../user_helpers/');

lambda0 = 1550e-9; % m
tfwhm = 0.5; % ps

T0 = tfwhm/(2*asech(1/sqrt(2))); % ps; 2*asech(1/sqrt(2))=1.7627

N = 1.1; % soliton number

%% Setup fiber parameters
sim.lambda0 = lambda0; % the central wavelength
%sim.gpu_yes = false;
sim.include_Raman = false; % turn off Raman

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

LD = T0^2/abs(fiber.betas(3)); % dispersion length
Ls = LD*pi/2; % soliton length
num_save = 300;
fiber.L0 = 3000; % m
sim.save_period = fiber.L0/num_save;

fiber.betas(4:end) = 0; % turn off higher-order dispersion (TOD, FOD, etc.)

%% Setup general parameters
Nt = 2^17; % the number of time points
time_window = 2000; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

%% Initial condition
Aeff = 1/fiber.SR;
initial_pulse = build_MMsoliton(tfwhm, fiber.betas(3), 1/Aeff, lambda0, time_window, 1, Nt, {'ifft',0}, N);

%% Propagate
prop_output = GMMNLSE_propagate(fiber,initial_pulse,sim);

%% Plot
% Time
figure;
h = plot(t,abs(prop_output.fields(:,:,end)).^2);
xlabel('t (ps)');
ylabel('Power (W)');
title('Field');
xlim([-10,10]);
set(h,'linewidth',2);
set(gca,'fontsize',14);

% Spectrum
figure;
h = plot(lambda,abs(fftshift(ifft(prop_output.fields(:,:,end)),1)).^2*c./lambda.^2);
xlabel('\lambda (nm)');
ylabel('Intensity (a.u.)');
xlim([1530,1570]);
title('Spectrum');
set(h,'linewidth',2);
set(gca,'fontsize',14);

% Comparison of time
figure;
[x,y] = meshgrid(t,prop_output.z/Ls);
pcolor(x,y,permute(abs(prop_output.fields(:,1,:)).^2,[3 1 2]));
shading interp; colormap(jet);
xlim([-10,10]);
xlabel('t (ps)');
ylabel('z/L_s');
title('Pulse evolution');
set(gca,'fontsize',14);

% Comparison of spectra
figure;
[x,y] = meshgrid((f-sim.f0),prop_output.z(2:end)/Ls);
tmp = 10*log10(permute(abs(fftshift(ifft(prop_output.fields(:,1,2:end)),1)).^2,[3 1 2])); tmp = tmp - max(tmp(:));
pcolor(x,y,tmp);
shading interp; colormap(jet); caxis([-20,0]);
xlim([-2,2]);
xlabel('\Deltaf (THz)');
ylabel('z/L_s');
title('Spectral evolution');
set(gca,'fontsize',14);

pulse_FWHM = zeros(length(prop_output.z),1);
for ii = 1:length(prop_output.z)
    threshold = max(abs(prop_output.fields(:,:,ii)).^2)/1.01;
    [~,~,tmp_pulse_width,~] = findpeaks(abs(prop_output.fields(:,:,ii)).^2,t,'MinPeakHeight',threshold,'WidthReference','halfheight','MinPeakProminence',threshold/2);
    pulse_FWHM(ii) = tmp_pulse_width;
end

spectrum_FWHM = zeros(length(prop_output.z),1);
for ii = 1:length(prop_output.z)
    %threshold = max(abs(fftshift(ifft(prop_output.fields(:,:,ii)))).^2)/1.01;
    %[~,~,tmp_spectrum_width,~] = findpeaks(abs(fftshift(ifft(prop_output.fields(:,:,ii)))).^2,f,'MinPeakHeight',threshold,'WidthReference','halfheight','MinPeakProminence',threshold/2);
    %spectrum_FWHM(ii) = tmp_spectrum_width;
    f0ii = sum(f.*abs(fftshift(ifft(prop_output.fields(:,:,ii)))).^2)/sum(abs(fftshift(ifft(prop_output.fields(:,:,ii)))).^2);
    spectrum_FWHM(ii) = sqrt(sum((f-f0ii).^2.*abs(fftshift(ifft(prop_output.fields(:,:,ii)))).^2)/sum(abs(fftshift(ifft(prop_output.fields(:,:,ii)))).^2));
end
figure;
yyaxis left;
plot(prop_output.z/Ls,pulse_FWHM,'linewidth',2,'Color','b');
set(gca,'YColor','b');
ylabel('Duration (ps)');
yyaxis right;
plot(prop_output.z/Ls,spectrum_FWHM,'linewidth',2);
ylabel('Bandwidth (THz)')
xlabel('z/L_s');
set(gca,'fontsize',20);