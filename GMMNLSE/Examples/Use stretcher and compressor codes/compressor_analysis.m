% This code analyzes the compressor's GDD and TOD behaviors.

close all; clearvars;

addpath('../../GMMNLSE algorithm/','../../user_helpers/');

%% Setup the pulse
Nt = 2^12;
time_window = 1; % ps
tfwhm = 0.05; % ps
pulse = build_MMgaussian(tfwhm, time_window, 1, 1, Nt);
dt = time_window/Nt; % ps
t = (-Nt/2:Nt/2-1)'*dt; % ps
wavelength0 = 1060; % nm
f = (-Nt/2:Nt/2-1)'/time_window + 299792.458/wavelength0; % THz

separation = 1e-4;

%% Dechirped by a prism compressor
alpha = pi/3;
dechirped_field_prism = pulse_compressor_single( 'prism',separation,[],wavelength0,t,pulse.fields,alpha,'N-SF11' );
[quardratic_phase_prism,cubic_phase_prism,~,quintic_phase_prism] = characterize_spectral_phase( f,fftshift(ifft(ifftshift(dechirped_field_prism,1)),1),7 );
GDD_prism = quardratic_phase_prism/separation*1e-3;
TOD_prism = cubic_phase_prism/separation*1e-3;
FOD_prism = quintic_phase_prism/separation*1e-3;

fprintf('    GDD (prism) = %6.4f(fs^2/mm)\n',GDD_prism);
fprintf('TOD/GDD (prism) = %6.4f(fs)\n',TOD_prism/GDD_prism);
fprintf('FOD/GDD (prism) = %6.4f(fs^2)\n',FOD_prism/GDD_prism);
fprintf('-----------------------------\n');

%% Dechirped by a grism compressor 1
alpha1 = 45*pi/180;
incident_angle = 35; % deg
grating_spacing = 1e-3/1000; % m
dechirped_field_grism1 = pulse_compressor_single( 'grism1',separation,incident_angle*pi/180,wavelength0,t,pulse.fields,grating_spacing,alpha1,'N-SF11' );
[quardratic_phase_grism1,cubic_phase_grism1,~,quintic_phase_grism1] = characterize_spectral_phase( f,fftshift(ifft(ifftshift(dechirped_field_grism1,1)),1),7 );
GDD_grism1 = quardratic_phase_grism1/separation*1e-3;
TOD_grism1 = cubic_phase_grism1/separation*1e-3;
FOD_grism1 = quintic_phase_grism1/separation*1e-3;

fprintf('    GDD (grism1) = %6.4f(fs^2/mm)\n',GDD_grism1);
fprintf('TOD/GDD (grism1) = %6.4f(fs)\n',TOD_grism1/GDD_grism1);
fprintf('FOD/GDD (grism1) = %6.4f(fs^2)\n',FOD_grism1/GDD_grism1);
fprintf('------------------------------\n');

%% Dechirped by a grism compressor 2
alpha2 = 69.1*pi/180;
incident_angle = 30; % deg
grating_spacing = 1e-3/1000; % m
dechirped_field_grism2 = pulse_compressor_single( 'grism2',separation,incident_angle*pi/180,wavelength0,t,pulse.fields,grating_spacing,alpha2,'fused silica' );
[quardratic_phase_grism2,cubic_phase_grism2,~,quintic_phase_grism2] = characterize_spectral_phase( f,fftshift(ifft(ifftshift(dechirped_field_grism2,1)),1),7 );
GDD_grism2 = quardratic_phase_grism2/separation*1e-3;
TOD_grism2 = cubic_phase_grism2/separation*1e-3;
FOD_grism2 = quintic_phase_grism2/separation*1e-3;

fprintf('    GDD (grism2) = %6.4f(fs^2/mm)\n',GDD_grism2);
fprintf('TOD/GDD (grism2) = %6.4f(fs)\n',TOD_grism2/GDD_grism2);
fprintf('FOD/GDD (grism2) = %6.4f(fs^2)\n',FOD_grism2/GDD_grism2);
fprintf('------------------------------\n');

%% Dechirped by a Treacy compressor
incident_angle = 60; % deg
grating_spacing = 1e-3/1000; % m
dechirped_field_Treacy = pulse_compressor_single( 'Treacy-t',separation,incident_angle*pi/180,wavelength0,t,pulse.fields,grating_spacing );
[quardratic_phase_Treacy,cubic_phase_Treacy,~,quintic_phase_Treacy] = characterize_spectral_phase( f,fftshift(ifft(ifftshift(dechirped_field_Treacy,1)),1),7 );
GDD_Treacy = quardratic_phase_Treacy/separation*1e-3;
TOD_Treacy = cubic_phase_Treacy/separation*1e-3;
FOD_Treacy = quintic_phase_Treacy/separation*1e-3;

fprintf('    GDD (Treacy-t) = %6.4f(fs^2/mm)\n',GDD_Treacy);
fprintf('TOD/GDD (Treacy-t) = %6.4f(fs)\n',TOD_Treacy/GDD_Treacy);
fprintf('FOD/GDD (Treacy-t) = %6.4f(fs^2)\n',FOD_Treacy/GDD_Treacy);
fprintf('--------------------------------\n');