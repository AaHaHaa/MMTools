function foutput = GMMNLSE_propagate_without_adaptive_rmc(fiber, initial_condition, sim, gain_rate_eqn)
%GMMNLSE_PROPAGATE_WITHOUT_ADAPTIVE_LMC Propagate an initial multimode pulse through an 
%arbitrary distance of an optical fiber without the use of adaptive-step 
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
%
%       Gain properties (for gain model 1,2,3; for rate-eqn gain model, see "gain_info.m") -->
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
%               saturation_intensity - for Taylor or multimode-Gaussian-gain model, the scale intensity in J/m^2
%                                      This is defined by h*f/(sigma*tau), where f is the center frequency,
%                                                                                sigma is the sum of the emission and absorption cross sections,
%                                                                                tau is the lifetime of the higher energy level for population inversion,
%                   OR
%               saturation_energy - for SM gain model, the scale energy in nJ
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
%           f0 - center frequency, in THz
%           deltaZ - step size, in m
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
%       Algorithms to use -->
%
%           gpu_yes - 1(true) = GPU
%                     0(false) = CPU
%
%                     Whether or not to use the GPU. Using the GPU is HIGHLY recommended, as a speedup of 50-100x should be possible.
%
%           Raman_model - 0 = ignore Raman effect
%                         1 = Raman model approximated analytically by a single vibrational frequency of silica molecules
%                                 (Ch. 2.3, p.42, Nonlinear Fiber Optics (5th), Agrawal)
%                         2 = Raman model including the anisotropic contribution
%                                 ("Ch. 2.3, p.43" and "Ch. 8.5, p.340", Nonlinear Fiber Optics (5th), Agrawal)
%                                 For more details, please read "Raman response function for silica fibers", by Q. Lin and Govind P. Agrawal (2006)
%
%               include_sponRS - true or false; whether to include spontaneous Raman term or not
%
%           gain_model - 0 = no gain
%                        1 = Gaussian-gain model
%                        2 = Gain-rate-equation model: see "gain_info.m" for details
%
%       Others -->
%
%           pulse_centering - 1(true) = center the pulse according to the time window, 0(false) = do not
%                             The time delay will be stored in time_delay after running GMMNLSE_propagate().
%           num_photon_noise_per_bin - a scalar; include photon noise (typically one photon per spectral discretization bin)
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
% foutput.fields - (N, num_modes, num_save_points) matrix with the multimode field at each save point
% foutput.dt - time grid point spacing, to fully identify the field
% foutput.z - the propagation length of the saved points
% foutput.deltaZ - the (small) step size for each saved points
% foutput.betas - the [betas0,betas1] used for the moving frame
% foutput.t_delay - the time delay of each pulse which is centered in the time window during propagation
% foutput.seconds - time spent in the main loop
%
% For gain-rate-equation model:
%   foutput.Power.pump.forward - (1,1,num_save_points); the forward pump power along the fiber
%   foutput.Power.pump.backward - (1,1,num_save_points); the backward pump power along the fiber
%   *If N2 is exported,
%       foutput.N2 - (Nx,Nx,num_save_points); the doped ion density of the upper state

%% Consider only the last fields of the initial condition
initial_condition.fields = initial_condition.fields(:,:,end);

%% Always use MPA stepping algorithm for the random mode coupling
% Because in our computation of random mode coupling, the dispersion term
% needs to be pre-computed for significant performance improvement, the
% varying gain term can't be treated as a dispersion operator anymore;
% instead, it needs to be treated as a nonlinear term, which is done in our
% MPA algorithm.
sim.step_method = 'MPA';
sim.MPA.coeff = permute(MPA_AM_coeff(sim.MPA.M),[1,3,2]);

%% Check the validity of input parameters
if sim.save_period == 0
    sim.save_period = fiber.L0;
end

% MPA performs M "sim.deltaZ/M" in parallel
num_zSteps = fiber.L0/sim.deltaZ;
num_zPoints_persave = sim.save_period/sim.deltaZ;
num_saveSteps = fiber.L0/sim.save_period;

% Because of machine error, "rem" alone to determine divisibility isn't
% enough. eps() is included here.
if rem(num_zSteps,1) && rem(num_zSteps+eps(num_zSteps),1) && rem(num_zSteps-eps(num_zSteps),1)
    error('GMMNLSE_propagate:SizeIncommemsurateError',...
        'The large step size is %f m and the fiber length is %f m, which are not commensurate', sim.deltaZ, fiber.L0)
else
    num_zSteps = round(num_zSteps);
end

if rem(num_zPoints_persave,1) && rem(num_zPoints_persave+eps(num_zPoints_persave),1) && rem(num_zPoints_persave-eps(num_zPoints_persave),1)
    error('GMMNLSE_propagate:SizeIncommemsurateError',...
        'The large step size is %f m and the save period is %f m, which are not commensurate', sim.deltaZ, sim.save_period)
else
    num_zPoints_persave = round(num_zPoints_persave);
end

if rem(num_saveSteps,1) && rem(num_saveSteps+eps(num_saveSteps),1) && rem(num_saveSteps-eps(num_saveSteps),1)
    error('GMMNLSE_propagate:SizeIncommemsurateError',...
        'The save period is %f m and the fiber length is %f m, which are not commensurate', sim.save_period, fiber.L0)
else
    num_saveSteps = round(num_saveSteps);
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
     sim.cuda_sponRS,sim.cuda_num_operations_sponRS,...
     sim.cuda_MPA_psi_update,...
     sim.rmc.cuda_mypagemtimes] = setup_stepping_kernel(sim,Nt,num_modes,sim.rmc.model);
end

%% Pre-calculate the dispersion term
% The "omegas" here is actually (omega - omega0), omega: true angular frequency
%                                                 omega0: central angular frequency (=2*pi*f0)
if sim.gpu_yes
    dt = gpuArray(initial_condition.dt);
else
    dt = initial_condition.dt;
end
omegas = 2*pi*ifftshift(linspace(-floor(Nt/2), floor((Nt-1)/2), Nt))'/(Nt*dt); % in 1/ps, in the order that the fft gives

% The dispersion term in the GMMNLSE, in frequency space
[D_op,sim] = calc_D_op(fiber,sim,Nt,dt,omegas,initial_condition.fields);

%% Pre-calculate the factor used in GMMNLSE (the nonlinear constant)
c = 2.99792458e-4; % speed of ligth m/ps
if ~isfield(fiber,'n2') || isempty(fiber.n2)
    fiber.n2 = 2.3e-20; % m^2/W
end
prefactor = 1i*fiber.n2*(omegas+2*pi*sim.f0)/c; % m/W

%% Deal with deltaZ
% After this line, the code starts to use deltaZ in variables important in simulations.
if sim.gpu_yes
    sim.deltaZ = gpuArray(sim.deltaZ);
end

% Small step size for MPA
% It's the step size between parallelization planes
sim.small_deltaZ = sim.deltaZ/sim.MPA.M;

%% % We can pre-compute exp(D_op*z) and exp(-D_op*z) for all z
% Variable "z" for the computation of the dispersion operator
fftshift_omegas = fftshift(omegas,1);
spectrum = sum(abs(fftshift(ifft(initial_condition.fields),1)).^2,2);
omega0 = sum(fftshift_omegas.*spectrum)/sum(spectrum); % 2*pi*THz; the pulse center frequency (under shifted omega)
idx0 = find(fftshift_omegas>=omega0,1);
if idx0 > floor(Nt/2)
    idx0 = idx0 - floor(Nt/2);
else
    idx0 = idx0 + ceil(Nt/2);
end

z = sim.small_deltaZ*permute((0:sim.MPA.M)',[2,3,1]); % to make "Dz" into (N,num_modes,M+1)
Dz = (D_op-D_op(idx0,:)).*z;
D = struct('pos', exp( Dz),...
           'neg', exp(-Dz));

rmc_idx = (1:num_modes+1:num_modes^2)' + permute((0:num_zSteps-1)*num_modes^2,[1,3,2]);
D_idx = repmat((idx0 + (0:num_modes-1)*Nt)',1,1,num_zSteps);
sim.rmc.matrices(rmc_idx) = D_op(D_idx);
sim.rmc.matrices = permute(sim.rmc.matrices,[1,2,4,3]); % make it (num_modes,num_modes,1,num_zSteps); the 3rd dimension is for the MPA's M-parallelization
sim.rmc.D = struct('pos', myexpm_( sim.rmc.matrices.*z,[],[],true,false,true),...
                   'neg', myexpm_(-sim.rmc.matrices.*z,[],[],true,false,true));
sim.rmc = rmfield(sim.rmc,'matrices'); % remove unused "matrices" to save some memory
sim.rmc.D = cell2struct([squeeze(num2cell(sim.rmc.D.pos,[1,2,3])),...
                         squeeze(num2cell(sim.rmc.D.neg,[1,2,3]))],{'pos','neg'},2);

%% Pre-compute the Raman response in frequency space
if ~isfield(fiber,'material')
    fiber.material = 'silica';
end
[fiber,haw,hbw] = Raman_model( fiber,sim,Nt,dt);

%% Pre-calculate the Gaussian gain term if necessary
[G,saturation_parameter] = Gaussian_gain(fiber,sim,omegas);

%% Work out the overlap tensor details
[SK_info, SRa_info, SRb_info] = calc_SRSK(fiber,sim,num_spatial_modes);

%% Include the shot noise: one photon per mode
initial_condition.fields = include_shot_noise(sim,omegas,Nt,dt,initial_condition.fields);

%% Include spontaneous Raman scattering
sim.include_sponRS = (sim.Raman_model ~= 0 && sim.include_sponRS);
if sim.include_sponRS
    sponRS_prefactor = spontaneous_Raman(Nt,dt,sim);
    haw_sponRS = haw; haw_sponRS = 1i*imag(haw_sponRS);
    hbw_sponRS = haw; hbw_sponRS = 1i*imag(hbw_sponRS);
else
    sponRS_prefactor = 0; % dummy variable
    haw_sponRS = [];
    hbw_sponRS = [];
end

%% Create a damped frequency window to prevent aliasing
sim.damped_freq_window = create_damped_freq_window(Nt);

%% Setup the exact save points
% We will always save the initial condition as well
save_points = int64(num_saveSteps + 1);
num_zPoints = num_zSteps + 1;

%% Run the step function over each step
run_start = tic;
if sim.gain_model == 2 % rate-equation-gain model
    if ispc
        sep_char = '\';
    else % unix
        sep_char = '/';
    end

    % Load the folder of rate-eqn gain including random mode coupling
    current_path = mfilename('fullpath');
    sep_pos = strfind(current_path,sep_char);
    upper_folder = current_path(1:sep_pos(end-1));
    addpath([upper_folder 'Gain_rate_eqn/']);
% -------------------------------------------------------------------------
    % For single mode, the computation of the gain amplification
    % factor is faster with CPU if the number of point < ~2^20.
    if sim.gpu_yes && num_modes > 1 % multimode
        gain_rate_eqn.cross_sections        = structfun(@(c) gpuArray(c),gain_rate_eqn.cross_sections,'UniformOutput',false);
        gain_rate_eqn.overlap_factor.signal = gpuArray(gain_rate_eqn.overlap_factor.signal);
        gain_rate_eqn.N_total               = gpuArray(gain_rate_eqn.N_total);
        gain_rate_eqn.FmFnN                 = gpuArray(gain_rate_eqn.FmFnN);
        gain_rate_eqn.GammaN                = gpuArray(gain_rate_eqn.GammaN);
    end

    % Set up the ASE data for the simulation later
    if gain_rate_eqn.include_ASE
        initial_condition.Power.ASE.forward  = ifftshift(initial_condition.Power.ASE.forward (:,:,end),1); % in the order of fft
        initial_condition.Power.ASE.backward = ifftshift(initial_condition.Power.ASE.backward(:,:,end),1); % in the order of fft
        % Put ASE initial condition into GPU for the GPU computation
        if sim.gpu_yes
            initial_condition.Power.ASE.forward  = gpuArray(initial_condition.Power.ASE.forward);
            initial_condition.Power.ASE.backward = gpuArray(initial_condition.Power.ASE.backward);
        end
    end

    first_iteration  = ~isfield(gain_rate_eqn,'saved_data') || isempty(gain_rate_eqn.saved_data);
    % The condition below is to check whether the current pulse propagation is the first iteration of a linear oscillator.
    % Since for a linear oscillator, it's in general a self-consistent system,
    % it's not necessary to get the result of the first propagation really accurately.
    if gain_rate_eqn.linear_oscillator && first_iteration
        gain_rate_eqn.tol = gain_rate_eqn.tol*10; % the tolerance of the iterations of the gain rate-eqn model is increased
    end
    % =====================================================================
    % Run the correct step function depending on the options chosen.
    % =====================================================================
    % For a linear oscillator, two propagations are coupled by their pump power, ASE power, and pulses; 
    % therefore, the first propagation is run normally in "GMMNLSE_rategain" 
    % and then use "GMMNLSE_rategain_linear_oscillator" for the latter propagations
    % to consider this coupling effect.
    if gain_rate_eqn.linear_oscillator && ~first_iteration
        %function_name = 'SteppingCaller_rategain_linear_oscillator_rmc'; % start the pulse propgation in a linear oscillator
        error('Random mode coupling is not supported for linear oscillators yet.');
    else
        function_name = 'SteppingCaller_rategain_rmc';
    end
    GMMNLSE_rategain_func = str2func(function_name);

    [A_out,T_delay_out,...
     Power,N2,...
     saved_data] = GMMNLSE_rategain_func(sim,gain_rate_eqn,...
                                         num_zPoints,save_points,num_zPoints_persave,...
                                         initial_condition,...
                                         prefactor,...
                                         SRa_info, SRb_info, SK_info,...
                                         omegas, D,...
                                         haw, hbw,...
                                         haw_sponRS, hbw_sponRS, sponRS_prefactor,...
                                         gain_rate_eqn.saved_data);
else % No gain, Gaussian gain
    [A_out,T_delay_out] = SteppingCaller_rmc(sim,...
                                            G, saturation_parameter,...
                                            num_zPoints,save_points,num_zPoints_persave,...
                                            initial_condition,...
                                            prefactor,...
                                            SRa_info, SRb_info, SK_info,...
                                            D,...
                                            haw, hbw,...
                                            haw_sponRS, hbw_sponRS, sponRS_prefactor);
end

% -------------------------------------------------------------------------
% Just to get an accurate timing, wait before recording the time
if sim.gpu_yes
    sim.betas = gather(sim.betas);
    wait(sim.gpuDevice.Device);
end
fulltime = toc(run_start);

%% Save the results in a struct
foutput = struct('z', sim.save_period*(0:num_saveSteps)',...
                 'fields', A_out,...
                 'dt', initial_condition.dt,...
                 'betas', sim.betas,...
                 'seconds', fulltime,...
                 't_delay', T_delay_out);
if sim.gain_model == 2 % rate-equation-gain model
    foutput.Power = Power; % pump and ASE powers
    if gain_rate_eqn.export_N2
        foutput.N2 = N2; % upper-state population
    end
    if gain_rate_eqn.reuse_data
        foutput.saved_data = saved_data; % the saved data for all pulses, ASE and pump powers for the future oscillator roundtrips
    end
end

end