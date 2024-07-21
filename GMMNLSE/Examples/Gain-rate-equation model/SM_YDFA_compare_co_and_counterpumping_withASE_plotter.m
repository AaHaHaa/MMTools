clearvars; close all;

addpath('../../');

filename = 'SM_YDFA_compare_co_and_counterpumping.mat';

load(filename);

% nonlinear phase
distance = 0:sim.save_period:fiber.L0;
nonlinear_phase_copumping = accumulated_nonlinear_phase(fiber.L0,1/fiber.SR,sim.f0,output_field{1}.fields,distance,output_field{1}.dt);
nonlinear_phase_counterpumping = accumulated_nonlinear_phase(fiber.L0,1/fiber.SR,sim.f0,output_field{2}.fields,distance,output_field{2}.dt);
fprintf('nonlinear phase (copumping): %6.4f\n',nonlinear_phase_copumping);
fprintf('nonlinear phase (counterpumping): %6.4f\n',nonlinear_phase_counterpumping);

energy_copumping   = permute(sum(trapz(abs(output_field{1}.fields).^2),2)*dt/1e3,[3 2 1]);
energy_counterpumping = permute(sum(trapz(abs(output_field{2}.fields).^2),2)*dt/1e3,[3 2 1]);

% Energy
figure;
h = plot(distance,[energy_copumping energy_counterpumping]);
legend('copumping','counterpumping');
xlabel('Propagation length (m)');
ylabel('Energy (nJ)');
title('Energy');
set(h,'linewidth',2);
set(gca,'fontsize',14);

% Pump
figure;
h = plot(distance,[output_field{1}.Power.pump.forward(:) output_field{2}.Power.pump.backward(:)]);
legend('copumping','counterpumping');
xlabel('Propagation length (m)');
ylabel('Power (W)');
title('Pump power for both pumping directions');
set(h,'linewidth',2);
set(gca,'fontsize',14);

% ASE
figure;
h1 = plot(distance,permute(sum([output_field{1}.Power.ASE.forward output_field{1}.Power.ASE.backward]*1e3)/time_window,[3 2 1]));
hold on;
h2 = plot(distance,permute(sum([output_field{2}.Power.ASE.forward output_field{2}.Power.ASE.backward]*1e3)/time_window,[3 2 1]));
hold off;
legend('copumping (forward)','copumping (backward)','counterpumping (forward)','counterpumping (backward)');
xlabel('Propagation length (m)');
ylabel('Power (mW)');
title('ASE power for both pumping directions');
set(h1,'linewidth',2);
set(h2,'linestyle','--',{'color'},{[0 0 1];[1 0 0]},'linewidth',2);
set(gca,'fontsize',14);

% N2
figure;
h = plot(distance,1-[output_field{1}.population(:) output_field{2}.population(:)]);
legend('copumping','counterpumping');
xlabel('Propagation length (m)');
ylabel('N_1');
title('N_1 for both pumping directions');
set(h,'linewidth',2);
set(gca,'fontsize',14);

c = 299792458e-12; % m/ps
f = (-N/2:N/2-1)'/N/dt+c/sim.lambda0;
lambda = c./f*1e9; % nm

c = 299792.458; % nm/ps
factor = c./lambda.^2; % change the spectrum from frequency domain into wavelength domain

% -------------------------------------------------------------------------
% copumping
% -------------------------------------------------------------------------
% Field
figure;
subplot(2,1,1);
h = plot(t,abs(output_field{1}.fields(:,:,end)).^2);
xlabel('Time (ps)');
ylabel('Power (W)');
title('The final output field of YDFA (copumping)');
set(h,'linewidth',2);
set(gca,'fontsize',14);

% Spectrum
subplot(2,1,2);
h = plot(lambda,abs(fftshift(ifft(output_field{1}.fields(:,:,end)),1)).^2.*factor);
xlabel('Wavelength (nm)');
ylabel('PSD (a.u.)');
title('The final output spectrum of YDFA (copumping)');
xlim([1010 1050]);
set(h,'linewidth',2);
set(gca,'fontsize',14);

% -------------------------------------------------------------------------
% counterpumping
% -------------------------------------------------------------------------
% Field
figure;
subplot(2,1,1);
h = plot(t,abs(output_field{2}.fields(:,:,end)).^2);
xlabel('Time (ps)');
ylabel('Power (W)');
title('The final output field of YDFA (counterpumping)');
set(h,'linewidth',2);
set(gca,'fontsize',14);

% Spectrum
subplot(2,1,2);
h = plot(lambda,abs(fftshift(ifft(output_field{2}.fields(:,:,end)),1)).^2.*factor);
xlabel('Wavelength (nm)');
ylabel('PSD (a.u.)');
title('The final output spectrum of YDFA (counterpumping)');
xlim([1010 1050]);
set(h,'linewidth',2);
set(gca,'fontsize',14);