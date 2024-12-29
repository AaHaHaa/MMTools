% This code simulates a 12-plate PKLM of a Gaussian pulse.

clearvars; close all;

addpath('../../UPPE3D algorithm/','../../user_helpers/'); % add package paths
sim.gpuDevice.Index = 1; % choose which GPU to use if you have multiple GPUs: 1,2,3...

%% Setup simulation parameters
sim.lambda0 = 1030e-9; % m; the center wavelength
sim.include_Raman = false; % no Raman

fiber.material = 'N-SF11'; % for finding the refractive index in GNLSE3D_propagate()

% Please check this function for details.
[fiber,sim] = load_default_UPPE3D_propagate(fiber,sim); % load default parameters

%% Setup PLKM parameters
Plate_Thick = 0.5e-3; % m
Plate_spacing = 9e-3; % m
D = 9e-3; % m; distance between the focal point and the first plate
MFD0 = 120e-6; % m; mode-field diameter

%% Initial condition
spatial_window = 1e-3; % m
tfwhm = 0.31; % ps
time_window = 2; % ps
energy = ((20)*0.85)*1e3; % nJ
Nt = 2^8; % the number of time points
Nx = 2^6; % the number of spatial points
initial_condition = build_3Dgaussian(MFD0, spatial_window, tfwhm, time_window, energy, Nt, Nx);

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
n2_plate = (10)*1e-20; % nonlinear refractive index

n_freespace = ones(1,Nx,Nx); % air's n=1
n2_freespace = 0; % no nonlinearity

% =========================================================================
% Start the simulation of PLKM
% =========================================================================
num_save = 100;
num_plates = 12;
z_all = zeros(num_save,2*num_plates+1);
MFD_all = zeros(num_save,2*num_plates+1);

Frame(num_save,2*num_plates+1) = struct('cdata',[],'colormap',[]);

%% air0; before the beam hits the first plate
fiber.L0 = D; % m
fiber.n =  n_freespace;
fiber.n2 = n2_freespace;

sim.save_period = fiber.L0/num_save;

% Simulate the propagation
prop_output = UPPE3D_propagate(fiber,initial_condition,sim);

% Calculate MFD during propagation
z_all(:,1) = prop_output.z(2:end); % the first one is the input z=0, so ignore it
remove_noise_model = 2; % see calcMFD for details
MFD = squeeze(calcMFD(squeeze(prop_output.field(Nt/2,:,:,:)),spatial_window,remove_noise_model))*1e3; % mm
MFD_all(:,1) = MFD(2:end); % the first one is the input, so ignore it
fprintf('air %u\n',0);
num_to_plot = num_save;
fig = plotter([],...
              prop_output.field,...
              z_all(1:num_to_plot),MFD_all(1:num_to_plot),...
              Nt,Nx,spatial_window/Nx,lambda);

% Animation
Frame(:,1) = animator(Frame(:,1),...
                      prop_output.field,...
                      z_all(1:num_to_plot),MFD_all(1:num_to_plot),0,...
                      Nt,dt,Nx,spatial_window/Nx,lambda);

%% PLKM: plate, air, plate, air,......
for i = 1+(1:num_plates*2)
    initial_condition.field = prop_output.field(:,:,:,end);
    
    if mod(i,2) == 0 % plate
        fiber.L0 = Plate_Thick; % m
        fiber.n =  n_plate;
        fiber.n2 = n2_plate;
    else % mod(i,2) == 1; air
        fiber.L0 = Plate_spacing;
        fiber.n = n_freespace; % air
        fiber.n2 = n2_freespace; % no nonlinearity
    end
    
    sim.save_period = fiber.L0/num_save;
    
    % Simulate the propagation
    prop_output = UPPE3D_propagate(fiber,initial_condition,sim);
    
    % Calculate MFD during propagation
    z_all(:,i) = prop_output.z(2:end) + z_all(end,i-1); % the first one is the input z=0, so ignore it
    MFD = squeeze(calcMFD(squeeze(prop_output.field(Nt/2,:,:,:)),spatial_window,remove_noise_model))*1e3; % mm
    MFD_all(:,i) = MFD(2:end); % the first one is the input, so ignore it
    if mod(i,2) == 0 % plate
        fprintf('plate %u\n',i/2);
    else
        fprintf('air %u\n',floor(i/2));
    end
    num_to_plot = i*num_save;
    fig = plotter(fig,... ...
                  prop_output.field,...
                  z_all(1:num_to_plot),MFD_all(1:num_to_plot),...
                  Nt,Nx,spatial_window/Nx,lambda);

    % Animation
    Frame(:,i) = animator(Frame(:,i),...
                          prop_output.field,...
                          z_all(1:num_to_plot),MFD_all(1:num_to_plot),(i-1)*num_save,...
                          Nt,dt,Nx,spatial_window/Nx,lambda);
end

% Movie
implay(Frame(:),20);