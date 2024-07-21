% This code computes the 4-um Nd linear Mamyshev oscillator with only one
% gain fiber without any passive one.
% Note that the final pump power and the inversion, N1, are the same at
% each point of the gain fiber for co- and counter-pumping because the gain
% can't respond to the pulse in time and sees only the average effect
% including both the forward and backward propagating pulses.
close all; clearvars;

addpath('../../../GMMNLSE algorithm/','../../../user_helpers/');

%% Gain info
% Please find details of all the parameters in "gain_info.m" if not specified here.
% Note that the use of single spatial mode is different from multi-spatial modes.
% "reuse_data" and "linear_oscillator_model" are activated and some related parameters are set.
gain_rate_eqn.gain_medium = 'Nd'; % specify the gain medium
gain_rate_eqn.base_medium = 'silica'; % specify the base medium
gain_rate_eqn.core_diameter = 4.7; % um
gain_rate_eqn.cladding_diameter = 77; % um
gain_rate_eqn.core_NA = 0.16;
gain_rate_eqn.absorption_wavelength_to_get_N_total = 800; % nm
gain_rate_eqn.absorption_to_get_N_total = 0.6; % dB/m
gain_rate_eqn.pump_wavelength = 808; % nm
gain_rate_eqn.copump_power = 0; % W; it's set below
gain_rate_eqn.counterpump_power = 0; % W; it's set below
gain_rate_eqn.reuse_data = true; % For a ring or linear cavity, the pulse will enter a steady state eventually.
                                 % If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
gain_rate_eqn.linear_oscillator = true; % For a linear oscillator, there are pulses from both directions simultaneously, which will deplete the gain;
                                        % therefore, the backward-propagating pulses need to be taken into account.
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 50; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-3; % the tolerance for the above iterations
gain_rate_eqn.verbose = true; % show the information(final pulse energy) during iterations of computing the gain

%% Setup fiber parameters
% load_default_GMMNLSE_propagate loads MFD=6.2um.

% General parameters
sim.lambda0 = 920e-9; % m
sim.f0 = 2.99792458e-4/sim.lambda0; % THz
sim.dz = 1000e-6; % um
sim.save_period = 0.05; % m
sim.gpu_yes = false; % For only the fundamental mode, running with CPU is faster if the number of points is lower than 2^(~18).

% -------------------------------------------------------------------------
% -------------------------------- Arm (6um) ------------------------------
% -------------------------------------------------------------------------
% Gain fiber
sim_Gain = sim;
sim_Gain.gain_model = 2;
fiber_Gain.L0 = 1.5; % m; the length of the gain fiber
fiber_Gain.MFD = 4; % um; mode-field diameter

% Load default parameters like 
%
% loading fiber.betas and fiber.SR based on your multimode folder above
% sim.Raman_model = 1; Use isotropic Raman model
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.gpu_yes = true; Use GPU (default to true)
% ......
%
% Please check this function for details.
[fiber_Gain,sim_Gain] = load_default_GMMNLSE_propagate(fiber_Gain,sim_Gain,'single_mode');

fiber_cavity = [fiber_Gain fiber_Gain];

%% Setup general cavity parameters
max_rt = 100; % maximum number of roundtrips
Nt = 2^13; % the number of points
time_window = 50; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % 
OC = 0.8; % output coupler
loss= 0.5;
tol_convergence = 1e-3; % the tolerance of pulse energy convergence for the oscillator

%% calculate fiber betas from silica refractive index
% Sellmeier coefficients
material = 'fused silica';
[a,b] = Sellmeier_coefficients(material);
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
n_silica = n_from_Sellmeier(lambda/1e3);

fiber.betas = n_silica*2*pi./(lambda*1e-9);

%% Filter parameters
gaussexpo = 1;
plot_filter_yes = true;
spectral_filter = struct('bw',1.7, ... % bandwidth (nm)
                         'cw',{920,912}); % center wavelength (nm)

%% Setup initial conditions
tfwhm = 1; % ps
total_energy = 10; % nJ
pedestal_energy = 0.1; % nJ
prop_output = build_noisy_MMgaussian(tfwhm, inf, time_window, total_energy,pedestal_energy,1,Nt,0.01);
            
%% Saved field information
field = cell(1,max_rt);
N0 = cell(1,max_rt);
pump = cell(1,max_rt);
splice_z = cumsum([fiber_cavity.L0]);
filter_displacement = sim.save_period/25;
splice_z = [splice_z(1) splice_z(1)+filter_displacement splice_z(2)+filter_displacement];
output_field = zeros(Nt,1,max_rt);
max_save_per_fiber = 100;

%% Load gain parameters
L_air = 1; % 1 is the free-space length
c = 299792458; % m/s
v = 1/fiber_Gain.betas(2)*1e12; % velocity in the fiber
gain_rate_eqn.t_rep = sum([fiber_cavity.L0])/v + L_air/c; % s

gain_rate_eqn.saved_data = [];

gain_rate_eqn_copumping = gain_rate_eqn;
gain_rate_eqn_counterpumping = gain_rate_eqn;
gain_rate_eqn_copumping.copump_power = 10; % W
gain_rate_eqn_counterpumping.counterpump_power = gain_rate_eqn_copumping.copump_power; % W

gain_rate_eqn_copumping = gain_info( fiber_Gain,sim_Gain,gain_rate_eqn_copumping,ifftshift(lambda,1) );
gain_rate_eqn_counterpumping = gain_info( fiber_Gain,sim_Gain,gain_rate_eqn_counterpumping,ifftshift(lambda,1) );

%% Run the cavity simulation
func = analyze_sim;

% Initialize some parameters
output_energy = zeros(max_rt,1);
rt_num = 0;
pulse_survives = true;
while rt_num < max_rt
    time_delay = 0;
    current_z = 0;
    zn = 1;
    rt_num = rt_num +1;
    
    t_iteration_start = tic;
    cprintf('*[1 0.5 0.31]','Iteration %d', rt_num);
    
    % ---------------------------------------------------------------------
    % Counterpumping gain fiber (fiber arm 1)
    % ---------------------------------------------------------------------
    prop_output = GMMNLSE_propagate(fiber_Gain, prop_output, sim_Gain, gain_rate_eqn_counterpumping);
    
    gain_rate_eqn_copumping.saved_data = prop_output.saved_data;
    time_delay = time_delay + prop_output.t_delay(end);
    
    % Save the information
    if isstruct(prop_output.fields)
        prop_output.fields = prop_output.fields.forward;
    end
    [saved_field,saved_z_this_fiber] = func.extract_saved_field(prop_output.fields,max_save_per_fiber,current_z,prop_output.z);
    field{rt_num}(:,:,zn:zn+size(saved_field,3)-1) = saved_field;
    saved_z(zn:zn+size(saved_field,3)-1) = saved_z_this_fiber;
    % Save the gain info
    saved_N0 = func.extract_saved_field( prop_output.population,max_save_per_fiber,current_z,sim.save_period );
    N0{rt_num}(:,:,zn:zn+size(saved_field,3)-1) = saved_N0;
    saved_pump_backward = func.extract_saved_field( prop_output.Power.pump.backward,max_save_per_fiber,current_z,sim.save_period );
    pump{rt_num}(:,:,zn:zn+size(saved_field,3)-1) = saved_pump_backward;
    
    current_z = saved_z_this_fiber(end);
    zn = zn + size(saved_field,3)-1;
    
    % ---------------------------------------------------------------------
    % Loss and spectral filter after the fiber arm 1
    % ---------------------------------------------------------------------
    % Spectral filter
    if rt_num ~= 1
        close(fig_filter);
    end
    [prop_output,fig_filter] = gaussian_spectral_filter(prop_output, sim.f0, spectral_filter(1).cw, spectral_filter(1).bw, gaussexpo ,plot_filter_yes); % Filter
    % Loss
    prop_output.fields(:,:,end) = prop_output.fields(:,:,end)*sqrt(1-loss);

    % Save the field after the filter
    current_z = current_z + filter_displacement;
    zn = zn + 1;

    % ---------------------------------------------------------------------
    % Copumping gain fiber (fiber arm 2)
    % ---------------------------------------------------------------------
    prop_output = GMMNLSE_propagate(fiber_Gain, prop_output, sim_Gain, gain_rate_eqn_copumping);
    
    gain_rate_eqn_counterpumping.saved_data = prop_output.saved_data;
    time_delay = time_delay + prop_output.t_delay(end);
    
    % Save the information
    if isstruct(prop_output.fields)
        prop_output.fields = prop_output.fields.forward;
    end
    [saved_field,saved_z_this_fiber] = func.extract_saved_field(prop_output.fields,max_save_per_fiber,current_z,prop_output.z);
    field{rt_num}(:,:,zn:zn+size(saved_field,3)-1) = saved_field;
    saved_z(zn:zn+size(saved_field,3)-1) = saved_z_this_fiber;
    % Save the gain info
    saved_N0 = func.extract_saved_field( prop_output.population,max_save_per_fiber,current_z,sim.save_period );
    N0{rt_num}(:,:,zn:zn+size(saved_field,3)-1) = saved_N0;
    saved_pump_forward = func.extract_saved_field( prop_output.Power.pump.forward,max_save_per_fiber,current_z,sim.save_period );
    pump{rt_num}(:,:,zn:zn+size(saved_field,3)-1) = saved_pump_forward;
    
    current_z = saved_z_this_fiber(end);
    zn = zn + size(saved_field,3)-1;
    
    % ---------------------------------------------------------------------
    % Finish propagation
    % ---------------------------------------------------------------------
    saved_z = saved_z(1:zn);
    
    % ---------------------------------------------------------------------
    % Output coupler
    output_field(:,:,rt_num) = sqrt(OC)*prop_output.fields(:,:,end);
    prop_output.fields = sqrt(1-OC)*sqrt(1-loss)*prop_output.fields(:,:,end);
    
    % ---------------------------------------------------------------------
    % Spectral filter
    close(fig_filter);
    [prop_output,fig_filter] = gaussian_spectral_filter(prop_output, sim.f0, spectral_filter(2).cw, spectral_filter(2).bw, gaussexpo ,plot_filter_yes); % Filter
    
    % ---------------------------------------------------------------------
    % Energy of the output field
    output_energy(rt_num) = sum(trapz(abs(output_field(:,:,rt_num)).^2))*prop_output.dt/1e3; % energy in nJ

    % If the energies stop changing, then we're done!
    if rt_num ~= 1
        close(fig);
    end
    warning('off')
    output_field(:,:,rt_num) = pulse_tracker(output_field(:,:,rt_num));
    warning('on');
    [converged_yes,fig] = check_convergence( output_energy,output_field(:,:,rt_num),f,t,tol_convergence,true );
    
    % ---------------------------------------------------------------------
    % Display running time
    t_iteration_end = toc(t_iteration_start);
    t_iteration_spent = datevec(t_iteration_end/3600/24);
    fprintf(': %1u:%2u:%3.1f\n',t_iteration_spent(4),t_iteration_spent(5),t_iteration_spent(6));
    
    % ---------------------------------------------------------------------
    % Update the repetition rate based on "time_delay"
    gain_rate_eqn_counterpumping.t_rep = gain_rate_eqn.t_rep + time_delay*1e-12;
    gain_rate_eqn_copumping.t_rep = gain_rate_eqn_counterpumping.t_rep;
    
    % ---------------------------------------------------------------------
    % Plot
    if rt_num ~= 1
        close(fig_evolution);
    end
    fig_evolution = func.analyze_fields(t,f,field{rt_num},saved_z,splice_z);
    
    if rt_num ~= 1
        close(fig_gain);
    end
    pump_plot.forward = cat(3,zeros(1,1,length(saved_z)/2),pump{rt_num}(1,1,length(saved_z)/2+1:end));
    pump_plot.backward = cat(3,pump{rt_num}(1,1,1:length(saved_z)/2),zeros(1,1,length(saved_z)/2));
    fig_gain = func.analyze_gain(saved_z,splice_z,pump_plot,N0{rt_num});
    
    % ---------------------------------------------------------------------
    % Break if converged
    if converged_yes
        cprintf('blue','The field has converged!\n');
        break;
    end
    % Break if pulse dies
    if output_energy(rt_num) < 0.01
        disp('The pulse dies.');
        pulse_survives = false;
        break;
    end
end

%% Finish the simulation and save the data
% Clear reducdant parts of the data
field = field(1:rt_num);
N0 = N0(1:rt_num);
pump = pump(1:rt_num);
output_field = output_field(:,:,1:rt_num);
energy = output_energy(arrayfun(@any,output_energy)); % clear zero

% -------------------------------------------------------------------------
% Save the final output field
if pulse_survives
    save('Nd_linear_Mamyshev_oscillator_noASE.mat', 't','f','output_field','time_delay','energy',...
                                                    'saved_z','splice_z','field',...
                                                    'pump','N0','gain_rate_eqn_copumping','gain_rate_eqn_counterpumping',...
                                                    '-v7.3'); % saved mat file version
end
% -------------------------------------------------------------------------

close(fig,fig_filter,fig_evolution,fig_gain);