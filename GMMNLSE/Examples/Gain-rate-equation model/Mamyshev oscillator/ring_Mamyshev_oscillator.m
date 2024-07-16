% This code runs a ring Mamyshev oscillator with a rate-equation-gain model.

close all; clearvars;

addpath('../../../GMMNLSE algorithm/','../../../user_helpers/');

%% Gain info
% Please find details of all the parameters in "gain_info.m" if not specified here.
% Note that the use of single spatial mode is different from multi-spatial modes.
% "reuse_data" and "linear_oscillator_model" are activated and some related parameters are set.
gain_rate_eqn.cross_section_filename = 'Liekki Yb_AV_20160530.txt';
gain_rate_eqn.core_diameter = 6; % um
gain_rate_eqn.cladding_diameter = 125; % um
gain_rate_eqn.core_NA = 0.12;
gain_rate_eqn.absorption_wavelength_to_get_N_total = 920; % nm
gain_rate_eqn.absorption_to_get_N_total = 0.55; % dB/m
gain_rate_eqn.pump_wavelength = 976; % nm
gain_rate_eqn.copump_power = 2; % W
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.reuse_data = false; % For a ring or linear cavity, the pulse will enter a steady state eventually.
                                 % If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
gain_rate_eqn.linear_oscillator = false; % For a linear oscillator, there are pulses from both directions simultaneously, which will deplete the gain;
                                         % therefore, the backward-propagating pulses need to be taken into account.
gain_rate_eqn.tau = 840e-6; % lifetime of Yb in F_(5/2) state (Paschotta et al., "Lifetme quenching in Yb-doped fibers"); in "s"
gain_rate_eqn.export_N2 = true; % whether to export N2, the ion density in the upper state or not
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 50; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the above iterations
gain_rate_eqn.verbose = true; % show the information(final pulse energy) during iterations of computing the gain

%% Setup fiber parameters
% load_default_GMMNLSE_propagate loads default parameters like 
%
% loading fiber.betas and fiber.SR based on your multimode folder above
% fiber.MFD = 6.2; % um
% sim.Raman_model = 1; Use isotropic Raman model
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.gpu_yes = true; Use GPU (default to true)
% ......

% General parameters
sim.lambda0 = 1030e-9;
%sim.progress_bar = false;
sim.save_period = 0.1;
sim.gpu_yes = false; % For only the fundamental mode, running with CPU is faster if the number of points is lower than 2^(~18).
sim.gain_model = 2;
sim.progress_bar_name = 'Gain';
fiber.L0 = 3;
[fiber,sim] = load_default_GMMNLSE_propagate(fiber,sim,'single_mode');

%% Setup general cavity parameters
max_rt = 100;
N = 2^13; % the number of time/freq points
time_window = 50; % ps
dt = time_window/N;
f = sim.f0+(-N/2:N/2-1)'/(N*dt); % THz
t = (-N/2:N/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm
OC = 0.9; % output coupler
loss = 0.1; % coupling loss
tol_convergence = 1e-5;

%% Filter parameters
spectral_filter = struct('bw',4, ...    % bandwidth (nm)
                         'cw',{1025,1045}); % center wavelength (nm)

%% Setup initial conditions
tfwhm = 1; % ps
total_energy = 1; % nJ
pedestal_energy = 0.01; % nJ

prop_output = build_noisy_MMgaussian(tfwhm, inf, time_window, total_energy,pedestal_energy,1,N,0.01);

%% Saved field information
save_num = int64(fiber.L0*2/sim.save_period + 1);
save_num = double(save_num);
saved_z = zeros(1,save_num);
splice_z = cumsum([fiber.L0,fiber.L0]);
filter_displacement = sim.save_period/25;
splice_z = [splice_z(1) splice_z(1)+filter_displacement splice_z(2)+filter_displacement];
output_field = zeros(N,1,max_rt);
output_field2 = zeros(N,1,max_rt);

%% Load gain parameters
L_air = 1; % 1 is the free-space length
c = 299792458; % m/s
v = 1/fiber.betas(2)*1e12; % velocity in the fiber
gain_rate_eqn.t_rep = fiber.L0*2/v + L_air/c; % s

% arm 1: SPM
gain_rate_eqn.copump_power = gain_rate_eqn.copump_power/3;
gain_rate_eqn_SPM = gain_info( fiber,sim,gain_rate_eqn,ifftshift(lambda,1) );
% arm 2: GMNA
gain_rate_eqn_GMNA = gain_info( fiber,sim,gain_rate_eqn,ifftshift(lambda,1) );

%% Run the cavity simulation
func = analyze_sim;

% Initialize some parameters
output_energy = zeros(max_rt,1);
rt_num = 0;
pulse_survives = true;
while rt_num < max_rt
    current_z = 0;
    time_delay = 0;
    zn = 1;
    rt_num = rt_num +1;
    
    t_iteration_start = tic;
    cprintf('*[1 0.5 0.31]','Iteration %d', rt_num);
    
    % ---------------------------------------------------------------------
    % fiber arm 1: SPM
    % ---------------------------------------------------------------------
    prop_output = GMMNLSE_propagate(fiber, prop_output, sim, gain_rate_eqn_SPM);
    time_delay = time_delay + prop_output.t_delay(end);

    % Save the field information
    field(:,:,zn:zn+length(prop_output.z)-1) = prop_output.fields;
    saved_z(zn:zn+length(prop_output.z)-1) = current_z + prop_output.z;

    % Save the gain info
    N2(:,:,zn:zn+length(prop_output.z)-1) = prop_output.N2;
    pump(:,:,zn:zn+length(prop_output.z)-1) = prop_output.Power.pump.forward;

    current_z = prop_output.z(end);
    zn = zn + length(prop_output.z)-1;
    
    % ---------------------------------------------------------------------
    % Loss and spectral filter after the fiber arm 1
    % ---------------------------------------------------------------------
    output_field2(:,:,rt_num) = prop_output.fields(:,:,end);
    if rt_num ~= 1
        close(fig_filter);
    end
    [prop_output,fig_filter] = gaussian_spectral_filter(prop_output, sim.f0, spectral_filter(1).cw, spectral_filter(1).bw,3,true); % Filter

    prop_output.fields = prop_output.fields(:,:,end)*sqrt(1-loss);

    % Save the field after the filter
    current_z = current_z + filter_displacement;
    zn = zn + 1;
    
    % ---------------------------------------------------------------------
    % fiber arm 2: GMNA
    % ---------------------------------------------------------------------
    prop_output = GMMNLSE_propagate(fiber, prop_output, sim, gain_rate_eqn_GMNA);
    time_delay = time_delay + prop_output.t_delay(end);

    % Save the field information
    field(:,:,zn:zn+length(prop_output.z)-1) = prop_output.fields;
    saved_z(zn:zn+length(prop_output.z)-1) = current_z + prop_output.z;

    % Save the gain info
    N2(:,:,zn:zn+length(prop_output.z)-1) = prop_output.N2;
    pump(:,:,zn:zn+length(prop_output.z)-1) = prop_output.Power.pump.forward;

    current_z = prop_output.z(end);
    zn = zn + length(prop_output.z)-1;

    % ---------------------------------------------------------------------
    % Finish propagation
    % ---------------------------------------------------------------------
    saved_z = saved_z(1:zn);
    
    % Output couplier
    output_field(:,:,rt_num) = sqrt(OC)*prop_output.fields(:,:,end);
    prop_output.fields = sqrt(1-OC)*sqrt(1-loss)*prop_output.fields(:,:,end);
    
    % -----------------------------------------------------------------
    % Spectral filter (after arm 2)
    close(fig_filter);
    [prop_output,fig_filter] = gaussian_spectral_filter(prop_output, sim.f0, spectral_filter(2).cw, spectral_filter(2).bw,3,true); % Filter
    
    % -----------------------------------------------------------------
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
    gain_rate_eqn_SPM.t_rep = gain_rate_eqn.t_rep + time_delay*1e-12;
    gain_rate_eqn_GMNA.t_rep = gain_rate_eqn_SPM.t_rep;
    
    % ---------------------------------------------------------------------
    % Plot
    if rt_num ~= 1
        close(fig_evolution);
    end
    fig_evolution = func.analyze_fields(t,f,field,saved_z,splice_z);
    
    if rt_num ~= 1
        close(fig_gain);
    end
    pump_plot.forward = pump;
    pump_plot.backward = zeros(1,1,length(saved_z));
    fig_gain = func.analyze_gain(saved_z,splice_z,pump_plot,N2);
    
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
output_field = output_field(:,:,1:rt_num);
output_field2 = output_field2(:,:,1:rt_num);
energy = output_energy(arrayfun(@any,output_energy)); % clear zero

% -------------------------------------------------------------------------
% Save the final output field
save('ring_Mamyshev_oscillator.mat', 't','f','output_field','output_field2','time_delay','energy',...
                         'saved_z','splice_z','field',...
                         'N2','pump','gain_rate_eqn_SPM','gain_rate_eqn_GMNA',...
                         'fiber','sim',... % cavity parameters
                         '-v7.3'); % saved mat file version
% -------------------------------------------------------------------------

close(fig,fig_filter,fig_evolution,fig_gain);