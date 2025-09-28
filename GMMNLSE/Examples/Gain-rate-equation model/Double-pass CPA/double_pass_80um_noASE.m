% This code computes the double-pass amplifier.

close all; clearvars;

addpath('../../../GMMNLSE algorithm/','../../../user_helpers/');

pump_power_forward  = 70; % W
pump_power_backward = 0; % W

%% Gain info
% Please find details of all the parameters in "gain_info.m" if not specified here.
% Note that the use of single spatial mode is different from multi-spatial modes.
% "reuse_data" and "linear_oscillator_model" are activated and some related parameters are set.
gain_rate_eqn.gain_medium = 'Yb'; % specify the gain medium
gain_rate_eqn.base_medium = 'silica'; % specify the base medium
gain_rate_eqn.core_diameter = 80; % um
gain_rate_eqn.cladding_diameter = 200; % um
gain_rate_eqn.core_NA = 0.025;
gain_rate_eqn.absorption_wavelength_to_get_N_total = 976; % nm
gain_rate_eqn.absorption_to_get_N_total = 40; % dB/m
gain_rate_eqn.pump_wavelength = 976; % nm
gain_rate_eqn.copump_power = 0; % W; it's set below
gain_rate_eqn.counterpump_power = 0; % W; it's set below
gain_rate_eqn.reuse_data = true; % For a ring or linear cavity, the pulse will enter a steady state eventually.
                                 % If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
gain_rate_eqn.linear_oscillator = true; % For a linear oscillator, there are pulses from both directions simultaneously, which will deplete the gain;
                                        % therefore, the backward-propagating pulses need to be taken into account.
gain_rate_eqn.t_rep = 1/1e6; % Assume 1 MHz here; s; the time required to finish a roundtrip (the inverse repetition rate of the pulse)
                               % This gain model solves the gain of the fiber under the steady-state condition; therefore, the repetition rate must be high compared to the lifetime of the doped ions.
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 100; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-3; % the tolerance for the above iterations
gain_rate_eqn.verbose = true; % show the information(final pulse energy) during iterations of computing the gain

%% Setup fiber parameters
% load_default_GMMNLSE_propagate loads MFD=6.2um.

% General parameters
sim.lambda0 = 1030e-9; % m
sim.f0 = 2.99792458e-4/sim.lambda0; % THz
sim.dz = 500e-6; % um
sim.save_period = 0.05; % m
sim.gpuDevice.Index = 1; % choose which GPU to use if you have multiple GPUs: 1,2,3...
%sim.gpu_yes = false;
% -------------------------------------------------------------------------
% -------------------------------- Arm (6um) ------------------------------
% -------------------------------------------------------------------------
% Gain fiber
sim.gain_model = 2;
fiber.L0 = 0.95; % m; the length of the gain fiber
fiber.MFD = 64; % um; mode-field diameter

% Load default parameters like 
%
% loading fiber.betas and fiber.SR based on your multimode folder above
% sim.include_Raman = true; Consider Raman
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.gpu_yes = true; Use GPU (default to true)
% ......
%
% Please check this function for details.
[fiber,sim] = load_default_GMMNLSE_propagate(fiber,sim,'single_mode');

%% Setup general cavity parameters
max_iterations = 100; % maximum number of iterations
Nt = 2^17; % the number of points
time_window = 6e3; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

%% calculate fiber betas from silica refractive index
% Sellmeier coefficients
material = 'silica';
[a,b] = Sellmeier_coefficients(material);
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
n_silica = n_from_Sellmeier(lambda/1e3);

fiber.betas = n_silica*2*pi./(lambda*1e-9);

%% Setup initial conditions
tfwhm = 0.2; % ps
total_energy = 3; % nJ

pulse_lambda0 = 1028e-9;
f_now = c/sim.lambda0*1e-12;
f_pulse = c/pulse_lambda0*1e-12;
freq_shift = f_pulse - f_now;
seed_pulse = build_MMgaussian(tfwhm, time_window, total_energy, 1, Nt, {'ifft',freq_shift});

stretched_duration = 500; % ps
[~,stretched_field] = pulse_stretcher_addNormalDispersion( 'double-Offner',stretched_duration,30*pi/180,sim.lambda0*1e9,t,seed_pulse.fields,1e-3/1000,5,true );
input_pulse = seed_pulse;
input_pulse.fields = stretched_field;
input_pulse.Power.ASE.forward = zeros(Nt,1);
input_pulse.Power.ASE.backward = zeros(Nt,1);

input_pulse = {input_pulse,input_pulse};

%% Load gain parameters
gain_rate_eqn = gain_info( fiber,sim,gain_rate_eqn,ifftshift(lambda,1) );

gain_rate_eqn = {gain_rate_eqn,gain_rate_eqn}; % forward, backward

gain_rate_eqn{1}.saved_data = [];
gain_rate_eqn{2}.saved_data = [];

%% Run the simulation
% container for analysis functions
func = analyze_sim;

% Initialize some parameters
output_energy = zeros(max_iterations,1);
iterations_num = 0;

gain_rate_eqn{1}.copump_power = pump_power_forward; % W
gain_rate_eqn{1}.counterpump_power = pump_power_backward; % W
gain_rate_eqn{2}.copump_power = pump_power_backward; % W
gain_rate_eqn{2}.counterpump_power = pump_power_forward; % W

while iterations_num < max_iterations
    iterations_num = iterations_num +1;

    if iterations_num > 1
        previous_saved_data = prop_output.saved_data;
    end

    t_iteration_start = tic;
    cprintf('*[1 0.5 0.31]','Iteration %d', iterations_num);
    % -----------------------------------------------------------------
    % Propagation inside fibers
    for i = 1:2
        prop_output = GMMNLSE_propagate(fiber, input_pulse{i}, sim, gain_rate_eqn{i});

        if i == 1
            input_pulse{2}.fields = prop_output.fields.forward(:,:,end);
        end

        gain_rate_eqn{mod(i,2)+1}.saved_data = prop_output.saved_data;
    end

    % ---------------------------------------------------------------------
    % Energy of the output field
    output_energy(iterations_num) = sum(trapz(abs(prop_output.fields.forward(:,:,end)).^2))*prop_output.dt/1e3; % energy in nJ

    % If the energies stop changing, then we're done!
    if iterations_num ~= 1
        close(fig);
    end
    [converged_yes,fig] = check_convergence( output_energy,prop_output.fields.forward(:,:,end),f,t,1e-2 );

    % ---------------------------------------------------------------------
    % Average results from the previous iteration to avoid slow convergence
    if iterations_num > 5 && rand(1) > 0.5
        for k = 1:length(gain_rate_eqn{1}.saved_data.signal_fields)
            gain_rate_eqn{1}.saved_data.signal_fields_backward{k} = (abs(gain_rate_eqn{1}.saved_data.signal_fields_backward{k}) + abs(previous_saved_data.signal_fields_backward{k}))/2;
            gain_rate_eqn{1}.saved_data.Power_ASE_backward{k} = (gain_rate_eqn{1}.saved_data.Power_ASE_backward{k} + previous_saved_data.Power_ASE_backward{k})/2;
        end
    end
    
    % ---------------------------------------------------------------------
    % Display running time
    t_iteration_end = toc(t_iteration_start);
    t_iteration_spent = datevec(t_iteration_end/3600/24);
    fprintf(': %1u:%2u:%3.1f\n',t_iteration_spent(4),t_iteration_spent(5),t_iteration_spent(6));

    % ---------------------------------------------------------------------
    % Plot
    if iterations_num ~= 1
        close(fig_evolution);
    end
    fig_evolution = figure;
    if size(prop_output.fields.backward,1)~=1
        plot(prop_output.z,squeeze(sum(abs([prop_output.fields.forward,prop_output.fields.backward]).^2,1)*dt/1e3),'linewidth',2);
        xlabel('Propagation length (m)');
        ylabel('Pulse energy (nJ)');
        legend('Forward seed','Backward seed');
        set(gca,'fontsize',20);
    end

    if iterations_num ~= 1
        close(fig_gain);
    end
    fig_gain = func.analyze_gain(prop_output.z,[],prop_output.Power.pump,prop_output.population);

    % ---------------------------------------------------------------------
    % Break if converged
    if converged_yes
        cprintf('blue','The field has converged!\n');
        break;
    end
end

%% Finish the simulation
% Clear reducdant parts of the data
energy = output_energy(arrayfun(@any,output_energy)); % clear zeros

[Strehl_ratio,dechirped_FWHM,transform_limited_FWHM,peak_power] = analyze_field( t,f,prop_output.fields.forward(:,:,end),'Treacy-t',pi/6,1e-6 );

%% save
close all;
prop_output = rmfield(prop_output,'saved_data');
gain_rate_eqn{1} = rmfield(gain_rate_eqn{1},'saved_data');
gain_rate_eqn{2} = rmfield(gain_rate_eqn{2},'saved_data');
save('double_pass_noASE.mat','-v7.3');