function foutput = GMMNLSE_propagate_with_adaptive(fiber, initial_condition, sim, gain_rate_eqn)
%GMMNLSE_PROPAGATE_WITH_ADAPTIVE Propagate an initial multimode pulse 
%through an arbitrary distance of an optical fiber with an adaptive-step 
%method.
%
% -------------------------------------------------------------------------
%
%   "fiber" is a structure with the fields:
%
%       Basic properties -->
%
%           betas - a (?,nm) matrix; "nm" = num_spatial_modes if under scalar fields;
%                                    otherwise, "nm" can be both num_spatial_modes or 2*num_spatial_modes depending on whether there's birefringence.
%                   betas(i, :) = (i-1)th order dispersion coefficient for each mode, in ps^n/m
%
%           n2 - the nonlinear coefficient (default to 2.3e-20 if not set)
%
%           SR - SR tensor, in m^-2
%           L0 - length of fiber, in m
%           fiber_type - 'silica', 'chalcogenide', or 'ZBLAN' (default: 'silica')
%
%       Gain properties (for Gaussian-gain model; for rate-eqn gain model, see "gain_info.m") -->
%
%           dB_gain - the small-signal gain amplification of the pulse energy in dB;
%                     This is used to calculate the gain_coeff  (default to 30)
%           gain_coeff - small signal gain coefficient in m^-1, defined by "g" in A(z)=exp(gz/2)A(0)
%           gain_fwhm - FWHM of the gain spectrum, in m
%           gain_doped_diameter - the diameter of the doped core to compute the overlap integral between mode fields and the doped core and accurately find their gain
%
%           saturation ==>
%                   The following three parameters are used to calculate the saturation intensity,
%                   or you can set the saturation intensity directly.
%
%                   gain_tau - the lifetime of the upper energy level
%                   t_rep - the repetition rate of the pulse
%                   gain_cross_section - the absorption+emission cross sections
%
%               saturation_intensity - for multimode Gaussian gain model (J/m^2)
%                                      This is defined by h*f/(sigma*tau)*t_rep, where f is the center frequency,
%                                                                                sigma is the sum of the emission and absorption cross sections,
%                                                                                tau is the lifetime of the higher energy level for population inversion,
%                                                                                t_rep is the repetition rate of the pulse or laser
%                   OR
%               saturation_energy - for single-mode Gaussian gain model (nJ)
%                                   This is defined by "saturation_intensity*Aeff"
%
% -------------------------------------------------------------------------
%
%   "initial_condition" is a structure with the fields:
%
%       dt - time step
%       fields - initial field, in W^1/2, (N-by-num_modes).
%                If the size is (N-by-num_modes-by-S), then it will take the last S.
%
%                num_modes = num_spatial_modes if "sim.scalar = true"
%                num_modes = num_spatial_modes*2 (spatial modes + polarization modes-x,y) otherwise
%
% -------------------------------------------------------------------------
%
%   "sim" is a structure with the fields:
%
%       Basic settings -->
%
%           betas - the betas for the slowly varying approximation and the moving frame, 
%                   that is to say, fiber.betas([1 2],:) = fiber.betas([1 2],:) - sim.betas;
%                   (2,1) column vector;
%                   if not set, no "sim.betas", the simulation will be run relative to the first mode
%           midx - mode index; an integer array
%           f0 - center frequency, in THz
%           save_period - spatial period between saves, in m
%                         0 = only save input and output (save_period = fiber.L0)
%
%       MPA -->
%
%           MPA.M - parallel extent for MPA;
%                   1 is no parallelization,
%                   5-20 is recommended; there are strongly diminishing returns after 5-10.
%           MPA.n_tot_max - maximum number of iterations for MPA
%                           This doesn't really matter because if the step size is too large, the algorithm will diverge after a few iterations.
%           MPA.n_tot_min - minimum number of iterations for MPA
%           MPA.tol - tolerance for convergence for MPA
%                     Value of the average NRMSE between consecutive itertaions in MPA at which the step is considered converged.
%
%       Polarization included -->
%
%           scalar - 0(false) = consider polarization mode coupling
%                    1(true) = don't consider polarization mode coupling
%
%           *If under scalar field, the input field takes only the scalar fields, e.g., [mode1, mode2, mode3......].
%           *Otherwise, the input field of each polarization needs to be specified in the order of [mode1_+ mode1_- mode2_+ mode2_-......], where (+,-) can be (x,y) or any orthogonal modes.
%           *SRSK is always loaded in a dimension of num_spatial_modes^4. It's automatically calculated to its polarized version in the code.
%
%           ellipticity - the ellipticity of the polarization modes; Please refer to "Nonlinear Fiber Optics, eq (6.1.18) Agrawal" for the equations.
%                         0: linear polarization   -> (+,-)=(x,y)
%                         1: circular polarization -> (+,-)=(right,left)
%
%       Adaptive method -->
%
%           adaptive_dz.threshold - a scalar;
%                                       the accuracy used to determined whether to increase or decrease the step size.
%           adaptive_dz.max_dz - a scalar; the maximum adaptive step size
%
%       Algorithms to use -->
%
%           gpu_yes - 1(true) = GPU
%                     0(false) = CPU
%
%                     Whether or not to use the GPU. Using the GPU is HIGHLY recommended, as a speedup of 50-100x should be possible.
%
%           include_Raman - 0(false) = ignore Raman effect
%                           1(true) = Either (a) Raman model approximated analytically by a single vibrational frequency of silica molecules
%                                                (Ch. 2.3, p.42, Nonlinear Fiber Optics (5th), Agrawal)
%                                                Extensions to other materials are also included. Please check Raman_model().
%                                     or     (b) Raman model including the anisotropic contribution
%                                                ("Ch. 2.3, p.43" and "Ch. 8.5, p.340", Nonlinear Fiber Optics (5th), Agrawal)
%                                                For more details, please read "Raman response function for silica fibers", by Q. Lin and Govind P. Agrawal (2006)
%
%           gain_model - 0 = no gain
%                        1 = Gaussian-gain model
%                        2 = Gain-rate-equation model: see "gain_info.m" for details
%
%       Others -->
%
%           pulse_centering - 1(true) = center the pulse according to the time window, 0(false) = do not
%                             The time delay will be stored in time_delay after running GMMNLSE_propagate().
%           gpuDevice.Index - a scalar; the GPU to use
%           gpuDevice.Device - the output of MATLAB "gpuDevice(gpu_index)"
%           cuda_dir_path - path to the cuda directory into which ptx files will be compiled and stored
%           progress_bar - 1(true) = show progress bar, 0(false) = do not
%                          It'll slow down the code slightly. Turn it off for performance.
%           progress_bar_name - the name of the GMMNLSE propagation shown on the progress bar.
%                               If not set (no "sim.progress_bar_name"), it uses a default empty string, ''.
%
% =========================================================================
%
% All the cases of num_modes for input fields and betas:
%
%                 | field   betas   SRSK
%   --------------------------------------------
%    s            |   m       m      m
%    p            |   2m     m,2m    m
%
%   m: num_spatial_modes
%   s: scalar fields, p: polarized fields
%
%   If num_modes(betas) = m, it's assumed that polarization modes are
%   degenerate in betas and is expanded into 2m modes automatically in 
%   the code, that is, (m,2m)->2m.
%
%   SRSK is always in a dimension of num_spatial_modes^4.
%
% =========================================================================
% Output:
%   foutput.fields - (N, num_modes, num_save_points) matrix with the multimode field at each save point
%   foutput.dt - time grid point spacing, to fully identify the field
%   foutput.z - the propagation length of the saved points
%   foutput.dz - the (small) step size at each saved points
%   foutput.betas - the [betas0,betas1] used for the moving frame
%   foutput.t_delay - the time delay of each pulse which is centered in the time window during propagation
%   foutput.seconds - time spent in the main loop
%
% For gain-rate-equation model:
%   foutput.Power.pump.forward  - (1,1,num_save_points); the forward pump power along the fiber
%   foutput.Power.pump.backward - (1,1,num_save_points); the backward pump power along the fiber
%   foutput.population - (Nx,Nx,num_save_points); the doped ion density of the various states
%                        For single-mode, Nx=1.

%% Consider only the last fields of the initial condition
initial_condition.fields = initial_condition.fields(:,:,end);

%% Pick the stepping algorithm to use depending on the number of modes
if isscalar(sim.midx) % single mode
    sim.step_method = 'RK4IP';
else % multimode
    sim.step_method = 'MPA';
    sim.MPA.M = ceil(abs(sim.MPA.M)/2)*2; % force M to be even; if odd, make it a larger even number
    sim.MPA.coeff = permute(MPA_AM_coeff(sim.MPA.M),[1,3,2]);
end

%% Check the validity of input parameters
if sim.save_period == 0
    sim.save_period = fiber.L0;
end

num_saves_total = fiber.L0/sim.save_period;
if rem(num_saves_total,1) && rem(num_saves_total+eps(num_saves_total),1) && rem(num_saves_total-eps(num_saves_total),1)
    error('GMMNLSE_propagate:SizeIncommemsurateError',...
          'The save period is %f m and the fiber length is %f m, which are not commensurate', sim.save_period, fiber.L0)
else
    num_saves_total = round(num_saves_total);
end

% Error check on the dimensions (num_modes) of matrices
check_nummodes(sim, fiber, initial_condition.fields);

%% Some Parameters
% Get the numerical parameters from the initial condition.
[Nt, num_modes,~] = size(initial_condition.fields);
num_spatial_modes = size(fiber.SR,1);

% For polarized fields, the dimension of the input betas is allowed to be
% "num_spatial_modes", so it needs to be expanded into
% "2*num_spatial_modes" in the rest of the computation.
fiber = betas_expansion_including_polarization_modes(sim,fiber,num_modes);

%% Set up the GPU details
if sim.gpu_yes
    [sim.gpuDevice.Device,...
     sim.cuda_SRSK,sim.cuda_num_operations_SRSK,...
     sim.cuda_MPA_psi_update] = setup_stepping_kernel(sim,Nt,num_modes);
end

%% Pre-calculate the dispersion term
% The "omegas" here is actually (omega - omega0), omega: true angular frequency
%                                                 omega0: central angular frequency (=2*pi*f0)
if sim.gpu_yes
    dt = gpuArray(initial_condition.dt);
else
    dt = initial_condition.dt;
end
omegas = 2*pi*ifftshift(linspace(-floor(Nt/2), floor((Nt-1)/2), Nt))'/(Nt*dt); % in rad/ps, in the order that the ifft gives

% The dispersion term in the GMMNLSE, in frequency space
[D_op,sim] = calc_D_op(fiber,sim,Nt,dt,omegas,initial_condition.fields);

%% Pre-calculate the factor used in GMMNLSE (the nonlinear constant)
c = 2.99792458e-4; % speed of ligth m/ps
if ~isfield(fiber,'n2') || isempty(fiber.n2)
    fiber.n2 = 2.3e-20; % m^2/W
end
n2_prefactor = 1i*fiber.n2*(omegas+2*pi*sim.f0)/c; % m/W

%% Pre-compute the Raman response in frequency space
if ~isfield(fiber,'material')
    fiber.material = 'silica';
end
[fiber,haw,hbw] = Raman_model(fiber,sim,Nt,dt);

%% Pre-calculate the Gaussian gain term if necessary
[G,saturation_parameter] = Gaussian_gain(fiber,sim,omegas);

%% Work out the overlap tensor details
[SK_info, SRa_info, SRb_info] = calc_SRSK(fiber,sim,num_spatial_modes);

%% Create a damped frequency window to prevent aliasing
sim.damped_freq_window = create_damped_freq_window(Nt);

%% Setup the exact save points
% We will always save the initial condition as well
save_points = int64(num_saves_total + 1);
save_z = double(0:save_points-1)'*sim.save_period;

%% Modified shot-noise for noise modeling
At_noise = shot_noise(Nt,dt,sim,num_modes);
if sim.gpu_yes
    At_noise = gpuArray(At_noise);
end

%% Run the step function over each step
run_start = tic;
% -------------------------------------------------------------------------
if sim.gain_model == 2 % rate-equation-gain model
    % For fundamental mode, the computation of the gain amplification
    % factor is faster with CPU if the number of point < ~2^20.
    if sim.gpu_yes && num_modes > 1 % multimode
        gain_rate_eqn.cross_sections        = gpuArray(gain_rate_eqn.cross_sections);
        gain_rate_eqn.overlap_factor.signal = gpuArray(gain_rate_eqn.overlap_factor.signal);
        gain_rate_eqn.N_total               = gpuArray(gain_rate_eqn.N_total);
        gain_rate_eqn.FmFnN                 = gpuArray(gain_rate_eqn.FmFnN);
        gain_rate_eqn.GammaN                = gpuArray(gain_rate_eqn.GammaN);
    end
    
    [A_out,Power,...
     save_z,save_dz,...
     T_delay_out,...
     N] = SteppingCaller_adaptive_rategain(sim,gain_rate_eqn,...
                                           save_z,save_points,...
                                           initial_condition,...
                                           n2_prefactor,...
                                           omegas, D_op,...
                                           SK_info, SRa_info, SRb_info,...
                                           haw, hbw,...
                                           At_noise);
% -------------------------------------------------------------------------
else % No gain, Gaussian gain
    [A_out,...
     save_z,save_dz,...
     T_delay_out] = SteppingCaller_adaptive(sim,...
                                            G, saturation_parameter,...
                                            save_z,save_points,...
                                            initial_condition,...
                                            n2_prefactor,...
                                            D_op,...
                                            SK_info, SRa_info, SRb_info,...
                                            haw, hbw,...
                                            At_noise);
end

% -------------------------------------------------------------------------
% Just to get an accurate timing, wait before recording the time
if sim.gpu_yes
    sim.betas = gather(sim.betas);
    At_noise = gather(At_noise);
    wait(sim.gpuDevice.Device);
end
fulltime = toc(run_start);

%% Save the results in a struct
foutput = struct('z', save_z,...
                 'dz', save_dz,...
                 'fields', A_out,...
                 'dt', initial_condition.dt,...
                 'betas', sim.betas,...
                 'seconds', fulltime,...
                 't_delay', T_delay_out,...
                 'shot_noise',At_noise);
if sim.gain_model == 2
    foutput.Power = Power;
    foutput.population = N;
end

end