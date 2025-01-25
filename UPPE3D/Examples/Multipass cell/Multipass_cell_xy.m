% This code simulates a multipass cell and aims to duplicate Fig.4 of the
% reference paper:
%
% Hanna et al., "Nonlinear temporal compression in multipass cells:
% theory," J. Opt. Soc. Am. B 34(7), 1340-1347 (2017)
%
% This script uses the 3D-UPPE that uses full x-y dimension. For
% more-efficient modeling, pelase see its radially-symmetric version.

clearvars; close all;

addpath('../../UPPE3D algorithm/','../../user_helpers/'); % add package paths
sim.gpuDevice.Index = 1; % choose which GPU to use if you have multiple GPUs: 1,2,3...

%% Setup simulation parameters
sim.lambda0 = 1030e-9; % m; the center wavelength
%sim.include_Raman = false; % no Raman

fiber.material = 'silica'; % for finding the refractive index in UPPE3D_propagate()

% Please check this function for details.
[fiber,sim] = load_default_UPPE3D_propagate(fiber,sim); % load default parameters

%% Setup MPC parameters
num_roundtrip = 9;%18;
MPC_length = 0.529; % m
focal_R = 0.35; % m; radius of curvature of the mirror
focal_length = focal_R/2; % m
silica_thickness = 13e-3; % m
mirror_dispersion = -400e-6; % ps^2
MFD0 = 450e-6; % m; mode-field diameter; calculated from Eq.(5) in the reference paper

plate_z = [(MPC_length-silica_thickness*2)/2;...
           (MPC_length-silica_thickness*2)/2 + silica_thickness];
for i = 1:num_roundtrip*6
    switch mod(i,6)
        case {0,1,3,4}
            zi = silica_thickness;
        case {2,5}
            zi = MPC_length - silica_thickness*2;
    end
    plate_z = [plate_z;plate_z(end)+zi];
end

%% Initial condition
spatial_window = 8e-3; % m
tfwhm = 0.85; % ps
time_window = 10; % ps
energy = 30e3; % nJ
Nt = 2^8; % the number of time points
Nx = 2^8; % the number of spatial points
initial_condition = build_3Dgaussian_xy(MFD0, spatial_window, tfwhm, time_window, energy, Nt, Nx);

%% Setup general parameters
dt = time_window/Nt; % ps
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

%% Material properties
% Sellmeier coefficients
[a,b] = Sellmeier_coefficients(fiber.material);
% Calculate the index difference using the Sellmeier equation to generate n(lambda)
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));

% silica material properties
n_silica = n_from_Sellmeier(lambda/1e3).*ones(1,Nx,Nx); % refractive index
n2_silica = 2.3e-20; % m^2/W; nonlinear refractive index

% Air properties
n_air = ones(1,Nx,Nx); % air's n=1
n2_air = 0; % no nonlinearity

% =========================================================================
% Start the simulation of MPC
% =========================================================================
num_save = 5;
z_all = zeros(num_save,6*num_roundtrip+2);
MFD_all = zeros(num_save,6*num_roundtrip+2);

Frame(num_save,6*num_roundtrip+2) = struct('cdata',[],'colormap',[]);

%% Before hitting the first mirror
% 1. air
fiber.L0 = (MPC_length-silica_thickness)/2; % m
fiber.n =  n_air;
fiber.n2 = n2_air;

sim.save_period = fiber.L0/num_save;

% Simulate the propagation
prop_output = UPPE3D_propagate(fiber,initial_condition,sim);

% Calculate MFD during propagation
z_all(:,1) = prop_output.z(2:end); % the first one is the input z=0, so ignore it
remove_noise_model = 2; % see calcMFD for details
MFD = squeeze(calcMFD_xy(squeeze(sqrt(sum(abs(prop_output.field).^2,1))),spatial_window,remove_noise_model))*1e3; % mm
MFD_all(:,1) = MFD(2:end); % the first one is the input, so ignore it
fprintf('air %u\n',0);
num_to_plot = num_save;
fig = plotter_xy([],...
                 prop_output.field,...
                 z_all(1:num_to_plot),MFD_all(1:num_to_plot),...
                 Nt,Nx,spatial_window/Nx,lambda);

% Animation
Frame(:,1) = animator_xy(Frame(:,1),...
                         prop_output.field,...
                         z_all(1:num_to_plot),MFD_all(1:num_to_plot),0,...
                         Nt,dt,Nx,spatial_window/Nx,lambda,...
                         plate_z);

% 2. silica
fiber.L0 = silica_thickness; % m
fiber.n =  n_silica;
fiber.n2 = n2_silica;

sim.save_period = fiber.L0/num_save;

% Simulate the propagation
prop_output = UPPE3D_propagate(fiber,prop_output,sim);

% Calculate MFD during propagation
z_all(:,2) = prop_output.z(2:end) + z_all(end,1); % the first one is the input z=0, so ignore it
remove_noise_model = 2; % see calcMFD for details
MFD = squeeze(calcMFD_xy(squeeze(sqrt(sum(abs(prop_output.field).^2,1))),spatial_window,remove_noise_model))*1e3; % mm
MFD_all(:,2) = MFD(2:end); % the first one is the input, so ignore it
fprintf('silica %u\n',0);
num_to_plot = 2*num_save;
fig = plotter_xy(fig,...
                 prop_output.field,...
                 z_all(1:num_to_plot),MFD_all(1:num_to_plot),...
                 Nt,Nx,spatial_window/Nx,lambda);

% Animation
Frame(:,2) = animator_xy(Frame(:,2),...
                         prop_output.field,...
                         z_all(1:num_to_plot),MFD_all(1:num_to_plot),num_save,...
                         Nt,dt,Nx,spatial_window/Nx,lambda,...
                         plate_z);

%% PLKM: plate, air, plate, air,......
for i = 2+(1:num_roundtrip*6)
    initial_condition.field = prop_output.field(:,:,:,end);

    if ismember(mod(i-2,6),[1,4])
        Ef_reflected = add_thin_lens_phase_xy(ifft(initial_condition.field,[],1),initial_condition.dx,initial_condition.dy,fftshift(lambda,1)/1e9,focal_length);
        Ef_reflected = Ef_reflected.*exp(1i*mirror_dispersion/2*(2*pi*ifftshift(f-sim.f0,1)).^2);
        initial_condition.field = fft(Ef_reflected,[],1);
    end
    
    switch mod(i-2,6)
        case {0,1,3,4}
            fiber.L0 = silica_thickness; % m
            fiber.n =  n_silica;
            fiber.n2 = n2_silica;
        case {2,5}
            fiber.L0 = MPC_length - silica_thickness*2;
            fiber.n = n_air; % air
            fiber.n2 = n2_air; % no nonlinearity
    end
    
    sim.save_period = fiber.L0/num_save;
    
    % Simulate the propagation
    prop_output = UPPE3D_propagate(fiber,initial_condition,sim);
    
    % Calculate MFD during propagation
    z_all(:,i) = prop_output.z(2:end) + z_all(end,i-1); % the first one is the input z=0, so ignore it
    MFD = squeeze(calcMFD_xy(squeeze(sqrt(sum(abs(prop_output.field).^2,1))),spatial_window,remove_noise_model))*1e3; % mm
    MFD_all(:,i) = MFD(2:end); % the first one is the input, so ignore it
    switch mod(i-2,6)
        case {0,1,3,4}
            fprintf('silica %u\n',floor((i-2)/6));
        case {2,5}
            fprintf('air %u\n',floor((i-2)/6));
    end
    num_to_plot = i*num_save;
    fig = plotter_xy(fig,... ...
                     prop_output.field,...
                     z_all(1:num_to_plot),MFD_all(1:num_to_plot),...
                     Nt,Nx,spatial_window/Nx,lambda);

    % Animation
    Frame(:,i) = animator_xy(Frame(:,i),...
                             prop_output.field,...
                             z_all(1:num_to_plot),MFD_all(1:num_to_plot),(i-1)*num_save,...
                             Nt,dt,Nx,spatial_window/Nx,lambda,...
                             plate_z);
end

% Movie
implay(Frame(:),30);
exportVideo = VideoWriter('MPC');
exportVideo.FrameRate = 30;
open(exportVideo);
writeVideo(exportVideo, Frame(:));
close(exportVideo);

%% Dechirping
% Use the center field as the whole temporal profile for dechirping. This
% ignores spatiotemporal coupling.
output_mode_field = prop_output.field(:,1,1,end)*sqrt(pi*(MFD(end)/2*1e-3)^2);

theta_in = pi/6;
wavelength0 = sim.lambda0*1e6;
grating_spacing = 1e-6;
pulse_compressor( 'Treacy-t',theta_in,wavelength0,t,output_mode_field,grating_spacing,true );