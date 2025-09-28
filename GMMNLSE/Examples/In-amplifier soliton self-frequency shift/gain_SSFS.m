% This code shows the soliton self-frequency shift (SSFS) in an amplifier
% (CorActive's DCF_TM_10-128).

close all; clearvars;

addpath('../../GMMNLSE algorithm/','../../user_helpers/');

lambda0 = 2000e-9; % m
tfwhm = 0.5; % ps

T0 = tfwhm/(2*asech(1/sqrt(2))); % ps; 2*asech(1/sqrt(2))=1.7627

N = 1; % soliton number

%% Setup fiber parameters
sim.lambda0 = lambda0; % the central wavelength
sim.gpu_yes = false;
sim.gain_model = 2;

% Load default parameters like 
%
% loading fiber.betas and fiber.SR based on your multimode folder above
% sim.include_Raman = true; Consider Raman
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.gpu_yes = true; Use GPU (default to true)
% ......
%
% Please check this function for details.
[fiber,sim] = load_default_GMMNLSE_propagate([],sim); % load default parameters

fiber.MFD = 10;
num_save = 100;
fiber.L0 = 3; % m
sim.save_period = fiber.L0/num_save;

%% Gain info
% Please find details of all the parameters in "gain_info.m" if not specified here.
% Note that the use of single spatial mode is different from multi-spatial modes.
% Activating "reuse_data" or "linear_oscillator_model" requires setting other parameters.
% Check the example or "gain_info.m".
gain_rate_eqn.gain_medium = 'Tm'; % specify the gain medium
gain_rate_eqn.base_medium = 'silica'; % specify the base medium
gain_rate_eqn.core_diameter = 10; % um
gain_rate_eqn.cladding_diameter = 128; % um
gain_rate_eqn.core_NA = 0.22;
gain_rate_eqn.absorption_wavelength_to_get_N_total = 790; % nm
gain_rate_eqn.absorption_to_get_N_total = 5; % dB/m
gain_rate_eqn.pump_wavelength = 1550; % nm
gain_rate_eqn.copump_power = 1.5; % W
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.reuse_data = false; % For a ring or linear cavity, the pulse will enter a steady state eventually.
                                  % If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
gain_rate_eqn.linear_oscillator = false; % For a linear oscillator, there are pulses from both directions simultaneously, which will deplete the gain;
                                         % therefore, the backward-propagating pulses need to be taken into account.
gain_rate_eqn.t_rep = 1/20e6; % Assume 20 MHz here; s; the time required to finish a roundtrip (the inverse repetition rate of the pulse)
                              % This gain model solves the gain of the fiber under the steady-state condition; therefore, the repetition rate must be high compared to the lifetime of the doped ions.
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.sponASE_spatial_modes = []; % In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE. If empty like [], it's length(sim.midx).
gain_rate_eqn.max_iterations = 50; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the above iterations
gain_rate_eqn.verbose = false; % show the information(final pulse energy) during iterations of computing the gain

%% Setup general parameters
Nt = 2^13; % the number of time points
time_window = 100; % ps
dt = time_window/Nt;
f = sim.f0+(-Nt/2:Nt/2-1)'/(Nt*dt); % THz
t = (-Nt/2:Nt/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm

% Precompute some parameters related to the gain to save the computational time
% Check "gain_info.m" for details.
gain_rate_eqn = gain_info( fiber,sim,gain_rate_eqn,ifftshift(lambda,1) );

%% calculate fiber betas from silica refractive index
% This is important to correctly simulate the broadband situations.
% Taylor-series coefficients is only good in narrowband situations.

% Sellmeier coefficients
material = 'silica';
[a,b] = Sellmeier_coefficients(material);
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
n_silica = n_from_Sellmeier(lambda/1e3);

fiber.betas = n_silica*2*pi./(lambda*1e-9);

%% Initial condition
idx = find(lambda < lambda0*1e9,1);
p = polyfit(2*pi*f(idx-floor(Nt/100):idx+floor(Nt/100)),fiber.betas(idx-floor(Nt/100):idx+floor(Nt/100)),2);
beta2 = p(end-2)*2;

Aeff = 1/fiber.SR;
initial_pulse = build_MMsoliton(tfwhm, beta2, 1/Aeff, lambda0, time_window, 1, Nt, {'ifft',0}, N);

all_duration = 0.5:0.01:0.8;
all_lambda1 = zeros(length(duration),1);
all_energy = zeros(length(duration),1);
for i = 1:length(all_duration)
    func = calc_chirp;
    duration = all_duration(i); % ps
    spectrum_amplitude = ifft(initial_pulse.fields);
    [chirp,chirped_pulse] = func.General( duration,2*pi*ifftshift(f,1),spectrum_amplitude,1 );
    initial_pulse.fields = chirped_pulse;
    
    %% Propagate
    prop_output = GMMNLSE_propagate(fiber,initial_pulse,sim,gain_rate_eqn);
    
    %% Raman wavelength
    filtered_pulse = prop_output.fields(:,:,end);
    filtered_pulse(t<-1 | t>1) = 0; % Remove residual pump
    [~,lambda1] = calc_RMS(lambda,abs(fftshift(ifft(filtered_pulse),1)).^2./lambda.^2);

    %% Energy
    energy = trapz(t,abs(prop_output.fields(:,:,end)).^2)/1e3; % nJ
    
    %% Save the data
    all_lambda1(i) = lambda1;
    all_energy(i) = energy;
end

%% Save the data
save('SSFS.mat','all_lambda1','all_energy','all_duration');