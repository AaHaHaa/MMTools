% This code runs the multimode all-normal dispersion oscillator (ANDi) with
% an OM4 multimode fiber (with 50-um-diameter core) and a passive/Yb-doped 
% 10-um LMA fiber. The gain fiber is modeled with the rate equations. I put
% the OM4 fiber after the few-mode 10-um fiber to enhance the multimode
% nonlinear effect for this demonstration. Spatial filter and offset
% splicing are included to artificially couple various modes for multimode
% linear interactions, which play a crucial role in multimode modelocking.
%
% In addition, to demonstrate how to use the multimode functions, the
% system is assumed to be
%
%                                                                           output pulse
%                                                                                ^
%                                                                                |
%                                                                                |
%   --> 10-um passive fiber --> 10-um gain fiber --> OM4 passive fiber --> output coupler --
%   |                                                                            |
%   -------- spectral filter <-- spatial filter <-- saturable absorber <----------
%
% When the pulse couples back into the 10-um fiber, we assume that there is
% a perfect imaging system such that the end-facet of the OM4 fiber is
% imaged to the 10-um fiber. With this imaging system, the coupling matrix
% is based on the spatial profiles of two fiber modes.
% With this assumption, we use two equivalent methods for coupling:
%    (1) recompose_into_space() and decompose_into_modes()
%    (2) calc_filter_matrix() and spatial_filter_moderesolved()
%
% One important point I notice is that it's necessary to reduce the 
% adaptive-step threshold beyond the default 1e-3, such as 1e-4 or 1e-5, to
% guarantee convergence.
%
% The aim of this code is to demonstrate how to use several functions
% useful in multimode modelocking. For more details about multimode
% modelocking, please read the following papers:
% 1. Wright et al., "Spatiotemporal mode-locking in multimode fiber
%    lasers," Science 358, 94-97 (2017)
% 2. Wright et al., "Mechanisms of spatiotemporal mode-locking," Nat. Phys.
%    16, 565-570 (2020)

clearvars; close all;

%% Add the folders of multimode files and others
addpath('../../../GMMNLSE algorithm/','../../../user_helpers/'); % add where many GMMNLSE-related functions like  "GMMNLSE_propagate" is

% For 50-um OM4 highly-multimode fiber
% This code is only for demonstration purpose, so ony 6 modes are considered.
num_modes_OM4 = 6;
fiber_OM4.MM_folder = '../../../Fibers/OM4_wavelength1030nm/';
fiber_OM4.betas_filename = 'betas.mat';
fiber_OM4.S_tensors_filename = sprintf('S_tensors_%umodes.mat',num_modes_OM4);

% For 10-um 3-mode/few-mode gain fiber
num_modes_10um = 3;
fiber_10um.MM_folder = '../../../Fibers/YB1200-10_125DC-PM_wavelength1030nm/';
fiber_10um.betas_filename = 'betas.mat';
fiber_10um.S_tensors_filename = sprintf('S_tensors_%umodes.mat',num_modes_10um);

%% Gain info
gain_rate_eqn.gain_medium = 'Yb'; % specify the gain medium
gain_rate_eqn.reuse_data = false; % For a ring or linear cavity, the pulse will enter a steady state eventually.
                                  % If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
gain_rate_eqn.linear_oscillator = false; % For a linear oscillator, there are pulses from both directions simultaneously, which will deplete the gain;
                                         % therefore, the backward-propagating pulses need to be taken into account.
gain_rate_eqn.core_diameter = 10; % um
gain_rate_eqn.cladding_diameter = 125; % um
gain_rate_eqn.core_NA = 0.08; % in fact, this is only used in single-mode
gain_rate_eqn.absorption_wavelength_to_get_N_total = 920; % nm
gain_rate_eqn.absorption_to_get_N_total = 1.7; % dB/m
gain_rate_eqn.pump_wavelength = 976; % nm
gain_rate_eqn.copump_power = 0.7; % W
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.downsampling_factor = 1; % an integer; downsample the eigenmode profiles to run faster
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 10; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the iteration
gain_rate_eqn.verbose = true; % show the information(final pulse energy) during iterations of computing the gain

%% Setup fiber parameters
% General parameters
sim.lambda0 = 1030e-9; % m; center wavelength
sim.f0 = 2.99792458e-4/sim.lambda0; % THz; center frequency
%sim.progress_bar = false;
sim.adaptive_dz.threshold = 1e-4; % reduce the adaptive-step threshold (default: 1e-3)

% -------------------------------------------------------------------------
% Passive fiber (10-um fiber)
% -------------------------------------------------------------------------
sim_10um = sim; % this loads the above general parameters
sim_10um.progress_bar_name = 'Fiber before gain (10um)';
sim_10um.save_period = 0.1; % m
fiber_10um.L0 = 20; % m; fiber length
[fiber_10um,sim_10um] = load_default_GMMNLSE_propagate(fiber_10um,sim_10um,'multimode'); % for the passive fiber before the gain fiber

% -------------------------------------------------------------------------
% Gain fiber (10-um Yb-doped double-cladding fiber)
% -------------------------------------------------------------------------
sim_Gain = sim; % this loads the above general parameters
sim_Gain.gain_model = 2; % use the rate-eqn-gain model
sim_Gain.progress_bar_name = 'Gain (10um)';
sim_Gain.save_period = 0.1; % m
fiber_Gain.MM_folder = fiber_10um.MM_folder;
fiber_Gain.betas_filename = fiber_10um.betas_filename;
fiber_Gain.S_tensors_filename = fiber_10um.S_tensors_filename;
fiber_Gain.L0 = 2; % m; fiber length
[fiber_Gain,sim_Gain] = load_default_GMMNLSE_propagate(fiber_Gain,sim_Gain,'multimode'); % for the gain fiber

% -------------------------------------------------------------------------
% Passive OM4 fiber
% -------------------------------------------------------------------------
sim_OM4 = sim; % this loads the above general parameters
sim_OM4.progress_bar_name = 'MMF (OM4)';
sim_OM4.save_period = 0.1; % m

fiber_OM4.L0 = 1; % m; fiber length
[fiber_OM4,sim_OM4] = load_default_GMMNLSE_propagate(fiber_OM4,sim_OM4,'multimode');

%% Setup general cavity parameters
max_rt = 500; % maximum roundtrips (in case it doesn't converge)
Nt = 2^10; % the number of points
time_window = 50; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm
OC = 0.8; % output coupling
loss = 0.5; % the total loss of the cavity
saturation_power = 2500; % the saturation power of the saturable absorber; W
moddepth = 0.99; % the modulation depth of the saturable absorber
Aeff = 1/fiber_10um.SR(1); % assume 10-um fundamental mode for calculating "saturation_'intensity'=saturation_power/Aeff"
tol_convergence = 1e-3; % the tolerance of the convergence of the output ANDi pulse energy

%% Spectral filter parameters
gaussexpo = 1;
plot_filter_yes = true;
spectral_filter = struct('bw',3, ... % bandwidth (nm)
                         'cw',1030); % center wavelength (nm)

%% Setup initial conditions
tfwhm = 3; % ps
total_energy = 3; % nJ

pulse_lambda0 = 1030e-9;
f_now = c/sim.lambda0*1e-12;
f_pulse = c/pulse_lambda0*1e-12;
freq_shift = f_pulse - f_now;

% Initial condition is a Gaussian pulse with three modes with equal energies
prop_output = build_MMgaussian(tfwhm, time_window, total_energy, num_modes_OM4, Nt, {'ifft',freq_shift});

%% Saved field information
L0 = fiber_10um.L0 + fiber_Gain.L0 + fiber_OM4.L0; % total fiber length
save_num = int64(fiber_10um.L0/sim_10um.save_period + fiber_Gain.L0/sim_Gain.save_period + fiber_OM4.L0/sim_OM4.save_period - 1); % the total number of fields to save
save_num = double(save_num);
saved_z = zeros(1,save_num); % the propagation distance to save
field = cell(3,max_rt);
pump = cell(1,max_rt); % the pump power to save
output_field = zeros(Nt,num_modes_OM4,max_rt); % the output pulse

%% Load gain parameters
L_air = 1; % 1 is the free-space length
c = 299792458; % m/s
v = 1/fiber_OM4.betas(2)*1e12; % velocity in the fiber

% Approximate roundtrip time for the cavity
% The actual value will be continuously updated during iterations
t_rep = L0/v + L_air/c; % s; the time required to finish a roundtrip (the inverse repetition rate of the pulse)
                        % This gain model solves the gain of the fiber under the steady-state condition; therefore, the repetition rate must be high compared to the lifetime of the doped ions.

gain_rate_eqn.t_rep = t_rep;
gain_rate_eqn = gain_info( fiber_Gain,sim_Gain,gain_rate_eqn,ifftshift(lambda,1) );

%% Spatial mode profiles
% Note: it makes computation easier if two mode profiles have the same
% spatial window and the same spatial sampling rate.

% Downsampling is required; otherwise, memory usage is too high
% Downsampling function handle
downsample2 = @(x,r) x(1:r:end,1:r:end,:);

% Load OM4 mode profiles
Nx_OM4 = 400; % this is known in advance (=length(data.x) below)
mode_profiles_OM4 = zeros(Nx_OM4,Nx_OM4,num_modes_OM4);
for midx = 1:num_modes_OM4
    data = load(sprintf('%smode%uwavelength10300.mat',fiber_OM4.MM_folder,midx));
    mode_profiles_OM4(:,:,midx) = data.phi; % 1/sqrt(um)
end
downsample_ratio = 4;
% Normalize the mode fields for later computations
norm_mode_profiles_OM4 = sqrt(sum(sum(mode_profiles_OM4.^2,1),2))*mean(diff(data.x*1e-6));
mode_profiles_OM4 = downsample2(mode_profiles_OM4./norm_mode_profiles_OM4,downsample_ratio); % 1/sqrt(m)
mode_profiles_OM4_x = downsample(data.x,downsample_ratio); % um

% Load 10-um-fiber mode profiles
Nx_10um = 800; % this is known in advance (=length(data.x) below)
mode_profiles_10um = zeros(Nx_10um,Nx_10um,num_modes_10um);
for midx = 1:num_modes_10um
    data = load(sprintf('%smode%uwavelength10300.mat',fiber_10um.MM_folder,midx));
    mode_profiles_10um(:,:,midx) = data.phi; % 1/sqrt(um)
end
% Normalize the mode fields for later computations
norm_mode_profiles_10um = sqrt(sum(sum(mode_profiles_10um.^2,1),2))*mean(diff(data.x*1e-6));
mode_profiles_10um = downsample2(mode_profiles_10um./norm_mode_profiles_10um,Nx_10um/Nx_OM4*downsample_ratio); % 1/sqrt(m)
mode_profiles_10um_x = downsample(data.x,Nx_10um/Nx_OM4*downsample_ratio); % um

% Coupling matrix from OM4 to the 10-um fiber
fiber_coupling_matrix_OM42Gain = calc_mode_coupling_matrix(mode_profiles_OM4_x,  mode_profiles_OM4,...
                                                           mode_profiles_10um_x, mode_profiles_10um);
% Coupling matrix from the 10-um fiber to OM4
fiber_coupling_matrix_Gain2OM4 = calc_mode_coupling_matrix(mode_profiles_10um_x, mode_profiles_10um,...
                                                           mode_profiles_OM4_x,  mode_profiles_OM4);

% Offset in offset splicing between 10-um and OM4 fibers
offset = 5; % um

% Transmission matrix of passing the spatial filter
spatial_hwhm = 5e-6; % um
filter_type = 'pinhole';
offcenter = [0,0]*1e-6; % m
% Because calc_filter_matrix() relies on the unit "m", we scale our inputs accordingly.
spatial_filter_matrix = calc_filter_matrix(mode_profiles_OM4, mode_profiles_OM4_x*1e-6, spatial_hwhm, offcenter, filter_type);

% Make the two mode profiles have the same spatial window size
ratio = abs(mode_profiles_10um_x(1)/mode_profiles_OM4_x(1));
resized_mode_profiles_10um = spatial_field_resizing(mode_profiles_10um, ratio);
resized_mode_profiles_10um_x = mode_profiles_10um_x/ratio;

%% Run the cavity simulation

func = analyze_sim; % the container of several analyzing functions

% Initialize some parameters
output_energy = zeros(max_rt,1);
rt_num = 0;
pulse_survives = true;
while rt_num < max_rt
    time_delay = 0;
    current_z = 0;
    zn = 1;
    rt_num = rt_num + 1;
    
    t_iteration_start = tic;
    cprintf('*[1 0.5 0.31]','Iteration %d', rt_num);

    % ---------------------------------------------------------------------
    % 10-um passive fiber
    % ---------------------------------------------------------------------
    % I assume we have an imaging system such that the end-facet of OM4 is
    % imaged to 10-um fiber.
    prop_output.fields = (fiber_coupling_matrix_OM42Gain*prop_output.fields(:,:,end).').';

    % Fiber propagation
    prop_output = GMMNLSE_propagate(fiber_10um, prop_output, sim_10um);
    
    time_delay = time_delay + prop_output.t_delay(end);
    field{1,rt_num} = prop_output.fields;
    saved_z(zn:zn+length(prop_output.z)-1) = prop_output.z;
    zn = zn + length(prop_output.z)-1;

    % ---------------------------------------------------------------------
    % 10-um gain fiber
    % ---------------------------------------------------------------------

    % Fiber propagation
    prop_output = GMMNLSE_propagate(fiber_Gain, prop_output, sim_Gain, gain_rate_eqn);
    
    time_delay = time_delay + prop_output.t_delay(end);
    field{2,rt_num} = prop_output.fields;
    saved_z(zn:zn+length(prop_output.z)-1) = saved_z(zn) + prop_output.z;
    zn = zn + length(prop_output.z)-1;
    
    % ---------------------------------------------------------------------
    % Offset splicing
    % ---------------------------------------------------------------------
    % Uncomment each and comment the other one to see their effects
    % (1) offset splicing
    % The 10-um field is converted into the spatial domain and then converted into the mode domain of the OM4 fiber
    full_field_txy = recompose_into_space(sim_Gain.gpu_yes, resized_mode_profiles_10um, prop_output.fields(:,:,end), sim_Gain.cuda_dir_path);
    % Offset the field of offset splicing
    full_field_txy = circshift(full_field_txy,ceil(offset/mean(diff(resized_mode_profiles_10um_x))),2);
    prop_output.fields = decompose_into_modes(sim_Gain.gpu_yes, mode_profiles_OM4, full_field_txy, mean(diff(mode_profiles_OM4_x*1e-6)), sim_Gain.cuda_dir_path);
    % (2) no offset splicing if directly calculating their coupling matrix based on mode profiles
    % Apply the fiber coupling matrix for coupling 10-um-fiber fields into OM4 fields.
    %prop_output = spatial_filter_moderesolved(prop_output,fiber_coupling_matrix_Gain2OM4);

    % ---------------------------------------------------------------------
    % Passive OM4 fiber
    % ---------------------------------------------------------------------
    % Fiber propagation
    prop_output = GMMNLSE_propagate(fiber_OM4, prop_output, sim_OM4);
    
    time_delay = time_delay + prop_output.t_delay(end);
    field{3,rt_num} = prop_output.fields;
    saved_z(zn:zn+length(prop_output.z)-1) = saved_z(zn) + prop_output.z;
    zn = zn + length(prop_output.z)-1;
    
    % -----------------------------------------------------------------
    % Output coupler
    output_field(:,:,rt_num) = sqrt(OC)*prop_output.fields(:,:,end);
    prop_output.fields = sqrt(1-OC)*prop_output.fields(:,:,end);
    
    % -----------------------------------------------------------------
    % Saturable absorber
    prop_output = saturable_absorber_action_3d(prop_output, saturation_power, moddepth, mode_profiles_OM4, mean(diff(mode_profiles_OM4_x*1e-6)), Aeff, sim_Gain.gpu_yes*false, sim_Gain.cuda_dir_path);
    %prop_output = saturable_absorber_action_simple(prop_output, saturation_power, moddepth);

    % -----------------------------------------------------------------
    % Spatial filter
    prop_output = spatial_filter_moderesolved(prop_output,spatial_filter_matrix);
    
    % -----------------------------------------------------------------
    % Spectral filter
    if rt_num ~= 1
        close(fig_filter); % close the previous figure and plot a new one
    end
    [prop_output,fig_filter] = gaussian_spectral_filter(prop_output, sim.f0, spectral_filter.cw, spectral_filter.bw, gaussexpo ,plot_filter_yes); % Filter
    
    prop_output.fields = sqrt(1-loss)*prop_output.fields(:,:,end);
    
    % -----------------------------------------------------------------
    % Energy of the output field
    output_energy(rt_num) = sum(trapz(abs(output_field(:,:,rt_num)).^2,1),2)*prop_output.dt/1e3; % energy in nJ

    % If the energies stop changing, then we're done!
    if rt_num ~= 1
        close(fig); % close the previous figure and plot a new one
    end
    warning('off')
    output_field(:,:,rt_num) = pulse_tracker(output_field(:,:,rt_num));
    warning('on');
    [converged_yes,fig] = check_convergence( output_energy,output_field(:,:,rt_num),f,t,tol_convergence,true );
    
    % ---------------------------------------------------------------------
    % Display running time of each roundtrip
    t_iteration_end = toc(t_iteration_start);
    t_iteration_spent = datevec(t_iteration_end/3600/24);
    fprintf(': %1u:%2u:%3.1f\n',t_iteration_spent(4),t_iteration_spent(5),t_iteration_spent(6));
    
    % ---------------------------------------------------------------------
    % Update the repetition rate based on "time_delay"
    % The pulse will shift in time with respect to the moving frame.
    % Since I implement the pulse centering function in the code,
    % I can use the information of "time_delay" in each roundtrip to calibrate to get the actual repetition rate.
    gain_rate_eqn.t_rep = t_rep + time_delay*1e-12;
    
    % ---------------------------------------------------------------------
    % Break if converged
    if converged_yes
        cprintf('blue','The field has converged!\n');
        break;
    end
    % Break if pulse dies
    if output_energy(rt_num) < 0.0001 % smaller than 0.01 nJ
        disp('The pulse dies.');
        pulse_survives = false;
        break;
    end
end

%% Finish the simulation and save the data
% Clear reducdant parts of the data
field = field(:,1:rt_num);
output_field = output_field(:,:,1:rt_num);
energy = output_energy(arrayfun(@any,output_energy)); % clear zero

% -------------------------------------------------------------------------
% Save the final output field
if pulse_survives
    save('MM_ANDi.mat', 't','f','lambda','output_field','energy',...
                        'saved_z','field',...
                        '-v7.3'); % saved mat file version
end
% -------------------------------------------------------------------------

close(fig,fig_filter);

%% Compress the pulse in the fundamental mode only
[Strehl_ratio,dechirped_FWHM,transform_limited_FWHM,peak_power] = analyze_field( t,f,output_field(:,1,end),'Treacy-t',pi/6,1e-3/1000,true );

%% Plot the output spatial profile
full_field_txy = recompose_into_space(sim_Gain.gpu_yes, mode_profiles_10um, output_field(:,:,end), sim_Gain.cuda_dir_path);

figure;
pcolor(mode_profiles_10um_x,mode_profiles_10um_x,squeeze(sum(abs(full_field_txy).^2,1)));
shading interp;
colormap(jet);