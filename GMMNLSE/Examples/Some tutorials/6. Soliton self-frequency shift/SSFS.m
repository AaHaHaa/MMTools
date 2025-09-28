% This code shows the soliton self-frequency shift (SSFS).

close all; clearvars;

addpath('../../../GMMNLSE algorithm/','../../../user_helpers/');

lambda0 = 1550e-9; % m
tfwhm = 0.03; % ps

T0 = tfwhm/(2*asech(1/sqrt(2))); % ps; 2*asech(1/sqrt(2))=1.7627

N = 1.3; % soliton number

%% Setup fiber parameters
sim.lambda0 = lambda0; % the central wavelength
sim.pulse_centering = false;
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
[fiber,sim] = load_default_GMMNLSE_propagate([],sim); % load default parameters

LD = T0^2/abs(fiber.betas(3)); % dispersion length
num_save = 100;
fiber.L0 = 1; % m
sim.save_period = fiber.L0/num_save;

%% Setup general parameters
Nt = 2^11; % the number of time points
time_window = 10; % ps
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
factor_correct_unit = (Nt*dt)^2/1e3; % to make the spectrum of the correct unit "nJ/THz"
                                     % "/1e3" is to make pJ into nJ
factor = c./lambda.^2; % change the spectrum from frequency domain into wavelength domain

% Time
figure;
h = plot(t/T0,abs(prop_output.fields(:,:,end)).^2);
xlabel('t/T_0');
ylabel('Power (W)');
title('Field');
set(h,'linewidth',2);
set(gca,'fontsize',14);

% Spectrum
figure;
h = plot(lambda,abs(fftshift(ifft(prop_output.fields(:,:,end)),1)).^2*factor_correct_unit.*factor);
xlabel('\lambda (nm)');
ylabel('Intensity (nJ/nm)');
xlim([1000,3000]);
title('Spectrum');
set(h,'linewidth',2);
set(gca,'fontsize',14);

% Comparison of time
figure;
[x,y] = meshgrid(t/T0,prop_output.z/LD);
pcolor(x,y,permute(abs(prop_output.fields(:,1,:)).^2,[3 1 2]));
shading interp; colormap(jet);
xlabel('t/T_0');
ylabel('z/L_D');
title('Pulse evolution');
set(gca,'fontsize',14);

% Comparison of spectra
figure;
[x,y] = meshgrid((f-sim.f0),prop_output.z(2:end)/LD);
tmp = 10*log10(permute(abs(fftshift(ifft(prop_output.fields(:,1,2:end)),1)).^2*factor_correct_unit,[3 1 2])); tmp = tmp - max(tmp(:));
pcolor(x,y,tmp);
shading interp; colormap(jet); caxis([-20,0]);
xlabel('\Deltaf (THz)');
ylabel('z/L_D');
title('Spectral evolution');
set(gca,'fontsize',14);

%% Animation
% Temporal and spectral evolutions
max_PowerT = max(abs(prop_output.fields(:)).^2);
max_PSD = max(abs(fftshift(ifft(prop_output.fields),1)).^2*factor_correct_unit.*factor,[],'all');
save_point = size(prop_output.fields,3);
Frame(save_point) = struct('cdata',[],'colormap',[]);
for i = 1:save_point
    figs = figure;
    subplot(1,2,1);
    plot(t/T0,abs(prop_output.fields(:,:,i)).^2/max_PowerT,'linewidth',2,'Color','b');
    xlabel('Time');
    ylabel('Power (norm.)');
    set(gca,'XTick',[0,106.4],'XTickLabel',{'0','\Deltat'});
    set(gca,'YTick',[0,1]);
    xlim([-10,120]);
    ylim([0,1]);
    set(gca,'fontsize',25);
    subplot(1,2,2);
    plot(lambda,abs(fftshift(ifft(prop_output.fields(:,:,i)),1)).^2*factor_correct_unit.*factor/max_PSD,'linewidth',2,'Color','r');
    xlabel('Wavelength');
    ylabel('PSD (norm.)');
    set(gca,'XTick',[1550,1845.8],'XTickLabel',{'\lambda_0','\lambda_0+\Delta\lambda'});
    set(gca,'YTick',[0,1]);
    xlim([1400,2000]);
    ylim([0,1]);
    set(gca,'fontsize',25);
    cpos = get(figs,'Position');
    set(figs,'Position',[400,290,840,420]);
    drawnow;

    set(figs,'Color',[1,1,1]);

    Frame(i) = getframe(figs);
    close(figs);
end
% Movie
implay(Frame,20);

exportVideo = VideoWriter('SSFS');
exportVideo.FrameRate = 20;
open(exportVideo);
writeVideo(exportVideo, Frame);
close(exportVideo);

%% Save the data
%save('SSFS.mat');