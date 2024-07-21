clearvars; close all;

filename = 'MM_YDFA_check_single_mode.mat';

load(filename);

%% Plot results
energy_Gaussianain  = permute(sum(trapz(abs(output_field{1}.fields).^2),2)*dt/1e3,[3 2 1]);
energy_SMrategain = permute(sum(trapz(abs(output_field{2}.fields).^2),2)*dt/1e3,[3 2 1]);
energy_MMrategain = permute(sum(trapz(abs(output_field{3}.fields).^2),2)*dt/1e3,[3 2 1]);

% Energy
figure;
plot(output_field{1}.z,[energy_GaussianGain energy_SMrategain energy_MMrategain]);
legend('Gaussian gain','SMrategain','MMrategain');
xlabel('Propagation length (m)');
ylabel('Energy (nJ)');
title('Energy');

c = 299792458e-12; % m/ps
f = (-N/2:N/2-1)'/N/dt+c/sim{1}.lambda0;
lambda = c./f*1e9;

c = 299792.458; % nm/ps
factor = c./lambda.^2; % change the spectrum from frequency domain into wavelength domain

% -------------------------------------------------------------------------
% Gaussian gain
% -------------------------------------------------------------------------
% Field
figure;
subplot(2,1,1);
plot(t,abs(output_field{1}.fields(:,:,end)).^2);
legend('mode 1','mode 2','mode 3');
xlabel('Time (ps)');
ylabel('Power (W)');
title('The final output field of YDFA (Gaussian gain)');

% Spectrum
subplot(2,1,2);
plot(lambda,abs(ifftshift(ifft(output_field{1}.fields(:,:,end)),1)).^2);
legend('mode 1','mode 2','mode 3');
xlabel('Wavelength (nm)');
ylabel('PSD (a.u.)');
title('The final output spectrum of YDFA (Gaussian gain)');
xlim([1000 1080]);

% -------------------------------------------------------------------------
% SM rategain
% -------------------------------------------------------------------------
% Field
figure;
subplot(2,1,1);
plot(t,abs(output_field{2}.fields(:,:,end)).^2);
xlabel('Time (ps)');
ylabel('Power (W)');
title('The final output field of YDFA (SM rate-equation gain)');

% Spectrum
subplot(2,1,2);
plot(lambda,abs(ifftshift(ifft(output_field{2}.fields(:,:,end)),1)).^2);
xlabel('Wavelength (nm)');
ylabel('PSD (a.u.)');
title('The final output spectrum of YDFA (SM rate-equation gain)');
xlim([1000 1080]);

% -------------------------------------------------------------------------
% MM rategain
% -------------------------------------------------------------------------
% Field
figure;
subplot(2,1,1);
plot(t,abs(output_field{3}.fields(:,:,end)).^2);
legend('mode 1','mode 2','mode 3');
xlabel('Time (ps)');
ylabel('Power (W)');
title('The final output field of YDFA (MM rate-equation gain)');

% Spectrum
subplot(2,1,2);
plot(lambda,abs(ifftshift(ifft(output_field{3}.fields(:,:,end)),1)).^2);
legend('mode 1','mode 2','mode 3');
xlabel('Wavelength (nm)');
ylabel('PSD (a.u.)');
title('The final output spectrum of YDFA (MM rate-equation gain)');
xlim([1000 1080]);

% =========================================================================
% Pump
figure;
plot(output_field{1}.z,[output_field{2}.Power.pump.forward(:) output_field{3}.Power.pump.forward(:)]);
legend('SM rate-equation gain','MM rate-equation gain');
xlabel('Propagation length (m)');
ylabel('Power (W)');
title('Pump power for SM and MM rate-equation gain models');