% This code simulates a 12-plate PKLM compressor of a Gaussian pulse.
% The medium is N-SF11.
% A PLKM-based compressor is used to establish a discrete waveguide for
% nonlinear pulse compressor, which, here, compresses a 310-fs pulse to <80
% fs.
%
% Please check our paper for details:
%   Wang et al., "Efficient temporal compression of 10-μJ pulses in
%   periodic layered Kerr media," Opt. Lett. 49(20), 5787-5790 (2024)
%
% This script employs the 3D-UPPE that uses full x-y dimension. For
% more-efficient modeling, pelase see its radially-symmetric version.

clearvars; close all;

addpath('../../UPPE3D algorithm/','../../user_helpers/'); % add package paths
sim.gpuDevice.Index = 1; % choose which GPU to use if you have multiple GPUs: 1,2,3...

%% Setup simulation parameters
sim.lambda0 = 1030e-9; % m; the center wavelength
sim.include_Raman = false; % no Raman

fiber.material = 'N-SF11'; % for finding the refractive index in UPPE3D_propagate()

% Please check this function for details.
[fiber,sim] = load_default_UPPE3D_propagate(fiber,sim); % load default parameters

%% Setup PLKM parameters
num_plates = 12;
Plate_thickness = 0.5e-3; % m
Plate_spacing = 10.3e-3; % m
D = 12e-3; % m; distance between the focal point and the first plate
MFD0 = 120e-6; % m; mode-field diameter

plate_z = D;
for i = 2:num_plates*2+1
    if mod(i,2) == 0
        zi = Plate_thickness;
    else
        zi = Plate_spacing;
    end
    plate_z = [plate_z;plate_z(end)+zi];
end

%% Initial condition
spatial_window = 2e-3; % m
tfwhm = 0.31; % ps
time_window = 2; % ps
energy = ((20)*0.85)*1e3; % nJ
Nt = 2^8; % the number of time points
Nx = 2^7; % the number of spatial points
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

% Plate material properties
n_plate = n_from_Sellmeier(lambda/1e3).*ones(1,Nx,Nx); % refractive index
n2_plate = (10)*1e-20; % m^2/W; nonlinear refractive index

% Air properties
n_air = ones(1,Nx,Nx); % air's n=1
n2_air = 0; % no nonlinearity

% =========================================================================
% Start the simulation of PLKM
% =========================================================================
num_save = 50;
z_all = zeros(num_save,2*num_plates+1);
MFD_all = zeros(num_save,2*num_plates+1);

Frame(num_save,2*num_plates+1) = struct('cdata',[],'colormap',[]);

%% air0; before the beam hits the first plate
fiber.L0 = D; % m
fiber.n =  n_air;
fiber.n2 = n2_air;

sim.save_period = fiber.L0/num_save;

% Simulate the propagation
prop_output = UPPE3D_propagate(fiber,initial_condition,sim);

% Calculate MFD during propagation
z_all(:,1) = prop_output.z(2:end); % the first one is the input z=0, so ignore it
remove_noise_model = 2; % see calcMFD_xy for details
MFD = squeeze(calcMFD_xy(squeeze(prop_output.field(Nt/2,:,:,:)),spatial_window,remove_noise_model))*1e3; % mm
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

%% PLKM: plate, air, plate, air,......
for i = 1+(1:num_plates*2)
    initial_condition.field = prop_output.field(:,:,:,end);
    
    if mod(i,2) == 0 % plate
        fiber.L0 = Plate_thickness; % m
        fiber.n =  n_plate;
        fiber.n2 = n2_plate;
    else % mod(i,2) == 1; air
        fiber.L0 = Plate_spacing;
        fiber.n = n_air; % air
        fiber.n2 = n2_air; % no nonlinearity
    end
    
    sim.save_period = fiber.L0/num_save;
    
    % Simulate the propagation
    prop_output = UPPE3D_propagate(fiber,initial_condition,sim);
    
    % Calculate MFD during propagation
    z_all(:,i) = prop_output.z(2:end) + z_all(end,i-1); % the first one is the input z=0, so ignore it
    MFD = squeeze(calcMFD_xy(squeeze(prop_output.field(Nt/2,:,:,:)),spatial_window,remove_noise_model))*1e3; % mm
    MFD_all(:,i) = MFD(2:end); % the first one is the input, so ignore it
    if mod(i,2) == 0 % plate
        fprintf('plate %u\n',i/2);
    else
        fprintf('air %u\n',floor(i/2));
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
implay(Frame(:),20);
exportVideo = VideoWriter('PLKM_xy');
exportVideo.FrameRate = 20;
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