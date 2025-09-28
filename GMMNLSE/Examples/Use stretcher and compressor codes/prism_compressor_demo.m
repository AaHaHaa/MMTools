% This code stretches the pulse by adding normal dispersion, followed by
% compression with a prism pair.
%
% This simple code demonstrates how to use the stretcher and the compressor code.
% Please go check their functions for details.

clearvars; close all;

addpath('../../user_helpers/','../../GMMNLSE algorithm/');

%% Single-grating Offner stretcher
% 500 fs is more forgiving with a single-grating Offner stretcher with aberration.
Nt = 2^11;
time_window = 10; % ps
tfwhm = 0.5; % ps
input = build_MMgaussian(tfwhm, time_window, 300e3, 1, Nt);
dt = time_window/Nt; % ps
t = (-Nt/2:Nt/2-1)'*dt; % ps

target_duration = 1; % ps
incident_angle = 30; % deg
wavelength0 = 1030; % nm
grating_spacing = 1e-3/1000; % m
R = 1; % m; radius of curvature of the mirror in an Offner stretcher/compressor

% stretched by an Offner single-grating stretcher
fprintf('Start stretching...... \n');
[stretched_offcenter,stretched_field,~,~,~,~,recover_info] = pulse_stretcher_addNormalDispersion( 'single-Offner',target_duration,incident_angle*pi/180,wavelength0,t,input.fields,grating_spacing,R,true,true,-1 );
title('Stretched pulse (single-grating Offner stretcher)');
fprintf('Finish stretching.\n\n');

% Recover the field directly by applying the phase reversely
fprintf('Start compressing by applying a sign-reversed phase...... \n');
recovered_field = pulse_recover_stretching_dechirping(stretched_field,recover_info,true,t);
fprintf('\n');

% compressed by a prism-type compressor
fprintf('Start compressing with a prism compressor...... \n');
alpha = pi/3;
[dechirped_separation,optimal_FWHM,dechirped_field] = pulse_compressor( 'prism',[],wavelength0,t,stretched_field,alpha,'N-SF10',true,true,false );
title('Dechirped pulse (prism dechirper)');
fprintf('\n');