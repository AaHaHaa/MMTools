% This code stretches the pulse by adding anomalous dispersion.
%
% It has two parts.
% 1.
% The first part stretches and compresses back the 30-fs pulse @1030 nm
% with different types of compressors. It's stretched with a Martinez 
% stretcher which introduces aberration.
% Its aberration can't be corrected with any types of compressors used
% here.
%
% 2.
% The second part stretches and compresses back the 30-fs pulse @1030 nm
% with different types of compressors. It's stretched with a Treacy
% stretcher which has no aberration due to no lens.
% You can see only the aberration-free double-grating Offner compressor is
% able to compress it back.
%
% This simple code demonstrates how to use the stretcher and the compressor code.
% Please go check their functions for details.

clearvars; close all;

addpath('../');

%% Martinez stretcher
Nt = 2^18;
time_window = 2000; % ps
tfwhm = 0.03; % ps
input = build_MMgaussian(tfwhm, time_window, 300e3, 1, Nt);
dt = time_window/Nt; % ps
t = (-Nt/2:Nt/2-1)'*dt; % ps

target_duration = 500; % ps
incident_angle = 30; % deg
wavelength0 = 1030; % nm
grating_spacing = 1e-3/1000; % m
f = 0.75; % m; focal length of the Martinez stretcher

% stretched by a Martinez stretcher
fprintf('Start stretching...... \n');
[stretched_l,stretched_field] = pulse_stretcher_addAnomalousDispersion( 'Martinez',target_duration,incident_angle*pi/180,wavelength0,t,input.fields,grating_spacing,f,true,true,-1 );
title('Stretched pulse (Martinez stretcher)');
fprintf('Finish stretching.\n\n');

% compressed by an Offner-type single-grating compressor
R = 2; % m; radius of curvature of the mirror in an Offner stretcher/compressor
fprintf('Start compressing with an Offner single-grating compressor...... \n');
[dechirped_offcenter,optimal_FWHM,dechirped_field] = pulse_compressor( 'Offner1',incident_angle*pi/180,wavelength0,t,stretched_field,grating_spacing,R,true,true,-1 );
title('Dechirped pulse (single-grating Offner dechirper)');
fprintf('\n');
% compressed by an Offner-type double-grating compressor
%R = 2; % m; radius of curvature of the mirror in an Offner stretcher/compressor
%compressed_offcenter = 0.8;
%fprintf('Start compressing with an Offner double-grating compressor...... \n');
%[dechirped_separation,optimal_FWHM,dechirped_field] = pulse_compressor( 'Offner2',incident_angle*pi/180,wavelength0,t,stretched_field,grating_spacing,R,compressed_offcenter,true,true,-1 );
%title('Dechirped pulse (double-grating Offner dechirper)');
%fprintf('\n');
% compressed by an aberration-free Offner-type double-grating compressor
R = 2; % m; radius of curvature of the mirror in an Offner stretcher/compressor
fprintf('Start compressing with an aberration-free Offner double-grating compressor...... \n');
[dechirped_separation,optimal_FWHM,dechirped_field] = pulse_compressor( 'Offner3',incident_angle*pi/180,wavelength0,t,stretched_field,grating_spacing,R,true,true,-1 );
title('Dechirped pulse (aberration-free Offner dechirper)');
fprintf('\n');

fprintf('====================\n\n');

%% Treacy stretcher
Nt = 2^18;
time_window = 2000; % ps
tfwhm = 0.03; % ps
input = build_MMgaussian(tfwhm, time_window, 300e3, 1, Nt);
dt = time_window/Nt; % ps
t = (-Nt/2:Nt/2-1)'*dt; % ps

target_duration = 500; % ps
incident_angle = 30; % deg
wavelength0 = 1030; % nm
grating_spacing = 1e-3/1000; % m

% stretched by a Treacy stretcher
fprintf('Start stretching...... \n');
[stretched_separation,stretched_field] = pulse_stretcher_addAnomalousDispersion( 'Treacy-t',target_duration,incident_angle*pi/180,wavelength0,t,input.fields,grating_spacing,true,true,-1 );
title('Stretched pulse (Treacy stretcher)');
fprintf('Finish stretching.\n\n');

% compressed by an Offner-type single-grating compressor
R = 2; % m; radius of curvature of the mirror in an Offner stretcher/compressor
fprintf('Start compressing with an Offner single-grating compressor...... \n');
[dechirped_offcenter,optimal_FWHM,dechirped_field] = pulse_compressor( 'Offner1',incident_angle*pi/180,wavelength0,t,stretched_field,grating_spacing,R,true,true,-1 );
title('Dechirped pulse (single-grating Offner dechirper)');
fprintf('\n');
% compressed by an Offner-type double-grating compressor
%R = 2; % m; radius of curvature of the mirror in an Offner stretcher/compressor
%compressed_offcenter = 0.8;
%fprintf('Start compressing with an Offner double-grating compressor...... \n');
%[dechirped_separation,optimal_FWHM,dechirped_field] = pulse_compressor( 'Offner2',incident_angle*pi/180,wavelength0,t,stretched_field,grating_spacing,R,compressed_offcenter,true,true,-1 );
%title('Dechirped pulse (double-grating Offner dechirper)');
%fprintf('\n');
% compressed by an aberration-free Offner-type double-grating compressor
R = 2; % m; radius of curvature of the mirror in an Offner stretcher/compressor
fprintf('Start compressing with an aberration-free Offner double-grating compressor...... \n');
[dechirped_separation,optimal_FWHM,dechirped_field] = pulse_compressor( 'Offner3',incident_angle*pi/180,wavelength0,t,stretched_field,grating_spacing,R,true,true,-1 );
title('Dechirped pulse (aberration-free Offner dechirper)');
fprintf('\n');