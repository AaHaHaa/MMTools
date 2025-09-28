% This code stretches the pulse by adding normal dispersion.
%
% It has three parts.
% 1.
% The first part stretches and compresses back the 500-fs pulse @1030 nm
% with different types of compressors. It's stretched with a sinlge-grating
% Offner stretcher which introduces aberration for easy setup.
% 2.
% The second part stretches and compresses back the 30-fs pulse @1030 nm
% with different types compressors. It's stretched with an aberration-free
% double-grating Offner stretcher.
% 3.
% The third part stretches and compresses back the 30-fs pulse @1030 nm
% with different types of compressors. It's stretched with a Martinez 
% stretcher which introduces aberration.
% It's interesting to see that the single-grating Offner compressor can
% compress the pulse back better than the double-grating one. The
% aberration introduced by the compressor compensates the one from the
% Martinez stretcher.
%
% You can play with it. You'll see that for the 30-fs case, it can't be 
% compressed back to 30 fs well if it's stretched with the single-grating 
% Offner stretcher. This results from the broadband nature of the pulse.
% In this case, an aberration-free Offner stretcher is preferred.
% It needs a double-grating Offner stretcher which features aberration-free.
%
% This simple code demonstrates how to use the stretcher and the compressor code.
% Please go check their functions for details.

clearvars; close all;

addpath('../../user_helpers/','../../GMMNLSE algorithm/');

%% Single-grating Offner stretcher
% 500 fs is more forgiving with a single-grating Offner stretcher with aberration.
Nt = 2^18;
time_window = 2000; % ps
tfwhm = 0.5; % ps
input = build_MMgaussian(tfwhm, time_window, 300e3, 1, Nt);
dt = time_window/Nt; % ps
t = (-Nt/2:Nt/2-1)'*dt; % ps

target_duration = 100; % ps
incident_angle = 30; % deg
wavelength0 = 1030; % nm
grating_spacing = 1e-3/1000; % m
R = 4; % m; radius of curvature of the mirror in an Offner stretcher/compressor

% stretched by an Offner single-grating stretcher
fprintf('Start stretching...... \n');
[stretched_offcenter,stretched_field,~,~,~,~,recover_info] = pulse_stretcher_addNormalDispersion( 'single-Offner',target_duration,incident_angle*pi/180,wavelength0,t,input.fields,grating_spacing,R,true,true,-1 );
title('Stretched pulse (single-grating Offner stretcher)');
fprintf('Finish stretching.\n\n');

% Recover the field directly by applying the phase reversely
fprintf('Start compressing by applying a sign-reversed phase...... \n');
recovered_field = pulse_recover_stretching_dechirping(stretched_field,recover_info,true,t);
fprintf('\n');

% compressed by a Treacy-type transmissive grating compressor
fprintf('Start compressing with a Treacy grating compressor...... \n');
[dechirped_separation,optimal_FWHM,dechirped_field] = pulse_compressor( 'Treacy-t',incident_angle*pi/180,wavelength0,t,stretched_field,grating_spacing,true,true,false,-1 );
title('Dechirped pulse (Treacy dechirper)');
fprintf('\n');
% compressed by an Offner-type single-grating compressor
fprintf('Start compressing with an Offner single-grating compressor...... \n');
[dechirped_offcenter,optimal_FWHM,dechirped_field] = pulse_compressor( 'Offner1',incident_angle*pi/180,wavelength0,t,stretched_field,grating_spacing,R,true,true,false,-1 );
title('Dechirped pulse (single-grating Offner dechirper)');
fprintf('\n');
% compressed by an Offner-type double-grating compressor
fprintf('Start compressing with an Offner double-grating compressor...... \n');
[dechirped_separation,optimal_FWHM,dechirped_field] = pulse_compressor( 'Offner2',incident_angle*pi/180,wavelength0,t,stretched_field,grating_spacing,R,-stretched_offcenter,true,true,false,-1 );
title('Dechirped pulse (double-grating Offner dechirper)');
fprintf('\n');

fprintf('====================\n\n');

%% Aberration-free double-grating Offner stretcher
% 30 fs isn't forgiving with a single-grating Offner stretcher with aberration.
Nt = 2^20;
time_window = 8000; % ps
tfwhm = 0.03; % ps
input = build_MMgaussian(tfwhm, time_window, 300e3, 1, Nt);
dt = time_window/Nt; % ps
t = (-Nt/2:Nt/2-1)'*dt; % ps

target_duration = 1000; % ps
incident_angle = 30; % deg
wavelength0 = 1030; % nm
grating_spacing = 1e-3/1000; % m
R = 2; % m; radius of curvature of the mirror in an Offner stretcher/compressor

f = (-Nt/2:Nt/2-1)'/time_window + 299792.458/wavelength0; % THz

% stretched by an Offner double-grating stretcher
fprintf('Start stretching with a double-grating Offner stretcher...... \n');
[stretched_offcenter,stretched_field] = pulse_stretcher_addNormalDispersion( 'double-Offner',target_duration,incident_angle*pi/180,wavelength0,t,input.fields,grating_spacing,R,true,true,-1 );
title('Stretched pulse (double-grating Offner stretcher)');
fprintf('Finish stretching.\n\n');

fprintf('----------\n');
characterize_spectral_phase( f,stretched_field,5,true );
fprintf('----------\n');

% compressed by a Treacy-type transmissive grating compressor
fprintf('Start compressing with a Treacy grating compressor...... \n');
[dechirped_separation,optimal_FWHM,dechirped_field] = pulse_compressor( 'Treacy-t',incident_angle*pi/180,wavelength0,t,stretched_field,grating_spacing,true,true,false,-1 );
title('Dechirped pulse (double-grating Offner stretcher)');
fprintf('\n');
fprintf('----------\n');
characterize_spectral_phase( f,dechirped_field,5,true );
fprintf('----------\n');
% compressed by an Offner-type single-grating compressor
%fprintf('Start compressing with an Offner single-grating compressor...... \n');
%[dechirped_offcenter,optimal_FWHM,dechirped_field] = pulse_compressor( 'Offner1',incident_angle*pi/180,wavelength0,t,stretched_field,grating_spacing,R,true,true,false,-1 );
%title('Dechirped pulse (double-grating Offner stretcher)');
%fprintf('\n');
% compressed by an Offner-type double-grating compressor
%fprintf('Start compressing with an Offner double-grating compressor...... \n');
%[dechirped_separation,optimal_FWHM,dechirped_field] = pulse_compressor( 'Offner2',incident_angle*pi/180,wavelength0,t,stretched_field,grating_spacing,R,-stretched_offcenter,true,true,false,-1 );
%title('Dechirped pulse (double-grating Offner stretcher)');
%fprintf('\n');

fprintf('====================\n');

% =========================================================================
% stretched by an Offner single-grating stretcher
fprintf('Start stretching with a single-grating Offner stretcher...... \n');
[stretched_offcenter,stretched_field] = pulse_stretcher_addNormalDispersion( 'single-Offner',target_duration,incident_angle*pi/180,wavelength0,t,input.fields,grating_spacing,R,true,true,-1 );
title('Stretched pulse (single-grating Offner stretcher)');
fprintf('Finish stretching.\n\n');

fprintf('----------\n');
characterize_spectral_phase( f,stretched_field,5,true );
fprintf('----------\n');

% compressed by a Treacy-type transmissive grating compressor
fprintf('Start compressing with a Treacy grating compressor...... \n');
[dechirped_separation,optimal_FWHM,dechirped_field] = pulse_compressor( 'Treacy-t',incident_angle*pi/180,wavelength0,t,stretched_field,grating_spacing,true,true,false,-1 );
title('Dechirped pulse (single-grating Offner stretcher)');
fprintf('\n');
fprintf('----------\n');
characterize_spectral_phase( f,dechirped_field,5,true );
fprintf('----------\n');
% compressed by an Offner-type single-grating compressor
%fprintf('Start compressing with an Offner single-grating compressor...... \n');
%[dechirped_offcenter,optimal_FWHM,dechirped_field] = pulse_compressor( 'Offner1',incident_angle*pi/180,wavelength0,t,stretched_field,grating_spacing,R,true,true,false,-1 );
%title('Dechirped pulse (single-grating Offner stretcher)');
%fprintf('\n');
% compressed by an Offner-type double-grating compressor
%fprintf('Start compressing with an Offner double-grating compressor...... \n');
%[dechirped_separation,optimal_FWHM,dechirped_field] = pulse_compressor( 'Offner2',incident_angle*pi/180,wavelength0,t,stretched_field,grating_spacing,R,-stretched_offcenter,true,true,false,-1 );
%title('Dechirped pulse (single-grating Offner stretcher)');
%fprintf('\n');

%% Martinez stretcher
Nt = 2^20;
time_window = 8000; % ps
tfwhm = 0.03; % ps
input = build_MMgaussian(tfwhm, time_window, 300e3, 1, Nt);
dt = time_window/Nt; % ps
t = (-Nt/2:Nt/2-1)'*dt; % ps

target_duration = 200; % ps
incident_angle = 30; % deg
wavelength0 = 1030; % nm
grating_spacing = 1e-3/1000; % m
f = 0.3; % m; focal length of the Martinez stretcher

% stretched by a Martinez stretcher
fprintf('Start stretching...... \n');
[stretched_l,stretched_field] = pulse_stretcher_addNormalDispersion( 'Martinez',target_duration,incident_angle*pi/180,wavelength0,t,input.fields,grating_spacing,f,true,true,-1 );
title('Stretched pulse (Martinez stretcher)');
fprintf('Finish stretching.\n\n');

% compressed by a Treacy-type transmissive grating compressor
fprintf('Start compressing with a Treacy grating compressor...... \n');
[dechirped_separation,optimal_FWHM,dechirped_field] = pulse_compressor( 'Treacy-t',incident_angle*pi/180,wavelength0,t,stretched_field,grating_spacing,true,true,false,-1 );
title('Dechirped pulse (Treacy dechirper)');
fprintf('\n');
% compressed by a Martinez stretcher when grating-to-lens separation is
% larger than the focal length of the lens
%fprintf('Start compressing with a Martinez stretcher...... \n');
%[dechirped_separation,optimal_FWHM,dechirped_field] = pulse_compressor( 'Martinez',incident_angle*pi/180,wavelength0,t,stretched_field,grating_spacing,f,true,true,false,-1 );
%title('Dechirped pulse (Martinez dechirper)');
%fprintf('\n');