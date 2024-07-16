% This code computes the 6-um bidirectionally pumped GMNA with a seed
% entering at each end.
% I assume they are the same seed, so the output will be two GMN pulses.
% My code can run bidirectional pumping with only one seed pulse, so I use
% the linear-oscillator scheme for this type of simulation.
% For a linear oscillator, a gain fiber is affected by pulses from both
% directions due to slow response of Yb ions and sees only the average 
% effect.

close all; clearvars;

addpath('../../../../GMMNLSE algorithm/','../../../../user_helpers/');

pump_power_forward  = 1.4; % W
pump_power_backward = 1.4; % W

%% Gain info
% Please find details of all the parameters in "gain_info.m" if not specified here.
% Note that the use of single spatial mode is different from multi-spatial modes.
% "reuse_data" and "linear_oscillator_model" are activated and some related parameters are set.
gain_rate_eqn.gain_medium = 'Yb'; % specify the gain medium
gain_rate_eqn.core_diameter = 6; % um
gain_rate_eqn.cladding_diameter = 125; % um
gain_rate_eqn.core_NA = 0.12;
gain_rate_eqn.absorption_wavelength_to_get_N_total = 920; % nm
gain_rate_eqn.absorption_to_get_N_total = 0.55; % dB/m
gain_rate_eqn.pump_wavelength = 976; % nm
gain_rate_eqn.copump_power = 0; % W; it's set below
gain_rate_eqn.counterpump_power = 0; % W; it's set below
gain_rate_eqn.reuse_data = true; % For a ring or linear cavity, the pulse will enter a steady state eventually.
                                 % If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
gain_rate_eqn.linear_oscillator = true; % For a linear oscillator, there are pulses from both directions simultaneously, which will deplete the gain;
                                        % therefore, the backward-propagating pulses need to be taken into account.
gain_rate_eqn.t_rep = 1/15e6; % Assume 15 MHz here; s; the time required to finish a roundtrip (the inverse repetition rate of the pulse)
                              % This gain model solves the gain of the fiber under the steady-state condition; therefore, the repetition rate must be high compared to the lifetime of the doped ions.
gain_rate_eqn.ignore_ASE = false;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 50; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-3; % the tolerance for the above iterations
gain_rate_eqn.verbose = true; % show the information(final pulse energy) during iterations of computing the gain

%% Setup fiber parameters
% load_default_GMMNLSE_propagate loads MFD=6.2um.

% General parameters
sim.lambda0 = 1030e-9; % m
sim.f0 = 2.99792458e-4/sim.lambda0; % THz
sim.dz = 1000e-6; % um
sim.save_period = 0.05;
sim.gpu_yes = false; % For only the fundamental mode, running with CPU is faster if the number of points is lower than 2^(~18).

% -------------------------------------------------------------------------
% -------------------------------- Arm (6um) ------------------------------
% -------------------------------------------------------------------------
% Gain fiber
sim.gain_model = 2;
fiber.L0 = 3; % m; fiber length
%fiber.MFD = 6; % um; mode-field diameter; this is the default settings loaded below

% Load default parameters like 
%
% loading fiber.betas and fiber.SR based on your multimode folder above
% sim.Raman_model = 1; Use isotropic Raman model
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.gpu_yes = true; Use GPU (default to true)
% ......
%
% Please check this function for details.
[fiber,sim] = load_default_GMMNLSE_propagate(fiber,sim,'single_mode');

%% Setup general cavity parameters
max_iterations = 100; % maximum number of iterations
N = 2^14; % the number of points
time_window = 50; % ps
dt = time_window/N;
f = sim.f0+(-N/2:N/2-1)'/(N*dt); % THz
t = (-N/2:N/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

%% calculate fiber betas from silica refractive index
% Sellmeier coefficients
material = 'fused silica';
[a,b] = Sellmeier_coefficients(material);
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
n_silica = n_from_Sellmeier(lambda/1e3);

fiber.betas = n_silica*2*pi./(lambda*1e-9);

%% Setup initial conditions
tfwhm = 0.5; % ps
total_energy = 1; % nJ
seed_wavelength = 1025e-9; % m
freq_shift = c/seed_wavelength*1e-12 - sim.f0;
Gaussian_pulse = build_MMgaussian(tfwhm, time_window, total_energy,1,N,{'ifft',freq_shift});
Gaussian_pulse.Power.ASE.forward = zeros(N,1);
Gaussian_pulse.Power.ASE.backward = zeros(N,1);

seed_forward  = Gaussian_pulse;
seed_backward = Gaussian_pulse;
input_pulse = {seed_forward,seed_backward};

%% Load gain parameters
gain_rate_eqn_forward  = gain_rate_eqn;
gain_rate_eqn_backward = gain_rate_eqn;

gain_rate_eqn_forward = gain_info( fiber,sim,gain_rate_eqn_forward,ifftshift(lambda,1) );
gain_rate_eqn_backward = gain_info( fiber,sim,gain_rate_eqn_backward,ifftshift(lambda,1) );

gain_rate_eqn = {gain_rate_eqn_forward,gain_rate_eqn_backward};

%% Run the simulation
% container for analysis functions
func = analyze_sim;

% Initialize some parameters
output_energy = zeros(max_iterations,1);
iterations_num = 0;
gain_rate_eqn{1}.saved_data = [];
gain_rate_eqn{2}.saved_data = [];

% I run the simulations from small pump power to the set values; otherwise,
% the pulse may blow up because of bidirectional pumping without considering
% another seed pulse. In practice, two seed pulses deplete the pump, but my
% code considers only a single pulse.
power_num = 1;
pump_power_forward_j  = linspace(pump_power_forward/power_num,pump_power_forward,power_num);
pump_power_backward_j = linspace(pump_power_backward/power_num,pump_power_backward,power_num);

for j = 1:power_num
    gain_rate_eqn{1}.copump_power = pump_power_forward_j(j); % W
    gain_rate_eqn{1}.counterpump_power = pump_power_backward_j(j); % W
    gain_rate_eqn{2}.copump_power = pump_power_backward_j(j); % W
    gain_rate_eqn{2}.counterpump_power = pump_power_forward_j(j); % W

    while iterations_num < max_iterations
        iterations_num = iterations_num +1;

        t_iteration_start = tic;
        cprintf('*[1 0.5 0.31]','Iteration %d', iterations_num);
        % -----------------------------------------------------------------
        % Propagation inside fibers
        for i = 1:2
            prop_output = GMMNLSE_propagate(fiber, input_pulse{i}, sim, gain_rate_eqn{i});

            gain_rate_eqn{mod(i,2)+1}.saved_data = prop_output.saved_data;
        end

        % ---------------------------------------------------------------------
        % Energy of the output field
        output_energy(iterations_num) = sum(trapz(abs(prop_output.fields.forward(:,:,end)).^2))*prop_output.dt/1e3; % energy in nJ

        % If the energies stop changing, then we're done!
        if iterations_num ~= 1
            close(fig);
        end
        [converged_yes,fig] = check_convergence( output_energy,prop_output.fields.forward(:,:,end),f,t,1e-3 );

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
            close(fig_ASE);
        end
        fig_ASE = func.analyze_ASE(f,prop_output.Power.ASE,prop_output.z);
        
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
end

close(fig,fig_ASE(2));
fig_ASE = fig_ASE(1);

%% Finish the simulation
% Clear reducdant parts of the data
energy = output_energy(arrayfun(@any,output_energy)); % clear zeros

[Strehl_ratio,dechirped_FWHM,transform_limited_FWHM,peak_power] = analyze_field( t,f,prop_output.fields.forward(:,:,end),'Treacy-t',pi/6,1e-6 );