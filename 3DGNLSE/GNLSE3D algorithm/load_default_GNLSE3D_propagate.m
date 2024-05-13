function [fiber,sim] = load_default_GNLSE3D_propagate( input_fiber,input_sim )
%LOAD_DEFAULT_GNLSE3D_PROPAGATE It loads the default settings for "fiber"
%and "sim" for different types of modes used.
%
%   If a user has specified some of the parameters of "fiber" and "sim",
%   user-defined one will be chosen instead of the default ones.
%
%   If you want to use a different one, specify it as fiber.xxx or sim.xxx 
%   and send them into this function. The output will use your parameters 
%   besides other default parameters a user doesn't specify.
%
%% ========================================================================
%   Because some parameters are correlated to each other, sometimes it's 
%   necessary to load input parameters first.
%   However, this makes loading default parameters a complicated function.
%
%   The procedure of loading default parameters is described below.
%   The order matters!
% -------------------------------------------------------------------------
%
%   <-- Uncorrelated parameters are loaded directly -->
%
%       sim.f0 - depend on input f0 or lambda0
%                If no input f0 or lambda0, f0=3e5/1030e-9 (THz)
%
% If there's user-defined one, use user's for the parameters below.
% Below I list the default values -- >
%
%       fiber.fiber_type = 'silica';
%       fiber.n2 = 2.3e-20;
%
%       sim.deltaZ = 1000e-6;
%       sim.save_period = 0;
%
%       sim.scalar = true;
%
%       sim.adaptive_deltaZ.DW_threshold = 1e-6;
%       sim.adaptive_deltaZ.threshold = 1e-6;
%
%       sim.gpu_yes = true;
%       sim.Raman_model = 1;
%
%       sim.pulse_centering = true;
%       sim.include_noise = false;
%       sim.gpuDevice.Index = 1;
%       sim.progress_bar = true;
%       sim.progress_bar_name = '';
%       sim.cuda_dir_path = [folder_of_this_function 'GNLSE3D/cuda'];
%
%
% For both single-mode and multimode -->
%
%   fiber.L0 = (1) input L0
%              (2) 2 (m)
%
%% ========================================================================
%   Details of each parameter
% -------------------------------------------------------------------------
%
% Example Use:
%
%    % User-defined parameters
%    fiber.L0 = 3;
%    
%    % Incorporate default settings
%    [fiber,sim] = load_default_GNLSE3D_propagate(fiber,[]); % single_mode
%
%    % If there are "sim" settings
%    sim.adaptive_deltaZ.model = 0;
%    [fiber,sim] =  load_default_GNLSE3D_propagate(fiber,sim); % single_mode
%
%    % Use only user-defined "sim", not "fiber"
%    [fiber,sim] = load_default_GNLSE3D_propagate([],sim); % single_mode
%
% -------------------------------------------------------------------------
%
%	Additional parameters:
%
%       input_sim.lambda0 - central wavelength, in m
%
% -------------------------------------------------------------------------
%
%   If both "lambda0" and "f0" are set by users, the final value will depend on "f0".
%
% -------------------------------------------------------------------------
%
%   "fiber" is a structure with the field:
%
%       Basic properties -->
%
%           n2 - the nonlinear coefficient (default to 2.3e-20 if not set)
%           L0 - length of fiber, in m
%
% -------------------------------------------------------------------------
%
%   "initial_condition" is a structure with the field:
%
%       dt - time step
%       field - initial field
%
% -------------------------------------------------------------------------
%
%   "sim" is a structure with the field:
%
%       Basic settings -->
%
%           betas - the betas for the slowly varying approximation and the moving frame
%           f0 - center frequency, in THz
%           save_period - spatial period between saves, in m
%                         0 = only save input and output (save_period = fiber.L0)
%
%       Polarization included -->
%
%           scalar - 0(false) = consider polarization mode coupling
%                    1(true) = don't consider polarization mode coupling
%
%       Adaptive method -->
%
%           adaptive_deltaZ.threshold - a scalar;
%                                       the accuracy used to determined whether to increase or decrease the step size.
%           adaptive_deltaZ.DW_threshold - a scalar;
%                                          the accuracy used to determined whether to increase or decrease the step size for the dispersoin + waveguide.
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
%       Others -->
%
%           pulse_centering - 1(true) = center the pulse according to the time window, 0(false) = do not
%                             The time delay will be stored in time_delay after running GNLSE3D_propagate().
%           include_noise - 1(true) = include photon noise (one photon per frequency band), 0(false) = do not
%           gpuDevice.Index - a scalar; the GPU to use
%           gpuDevice.Device - the output of MATLAB "gpuDevice(gpu_index)"
%           cuda_dir_path - path to the cuda directory into which ptx files will be compiled and stored
%           progress_bar - 1(true) = show progress bar, 0(false) = do not
%                          It'll slow down the code slightly. Turn it off for performance.
%           progress_bar_name - the name of the GNLSE3D propagation shown on the progress bar.
%                               If not set (no "sim.progress_bar_name"), it uses a default empty string, ''.
%
% =========================================================================

%% Current path (or the folder where this "load_default_GNLSE3D_propagate.m" is)
if ispc
    sep = '\';
else % unix
    sep = '/';
end
current_path = mfilename('fullpath');
sep_pos = strfind(current_path,sep);
upper_folder = current_path(1:sep_pos(end-1));

%% Default settings below:

% Supperss warnings generated by the function, "catstruct", due to there
% are definitely duplicate elements between default and input.
warning('off','catstruct:DuplicatesFound');

if ~exist('input_fiber','var')
    input_fiber = [];
end
if ~exist('input_sim','var')
    input_sim = [];
end

% -------------------------------------------------------------------------
% Set some default parameters early here because these parameters will be
% used for loading files if multimode.
% If single-mode, only "fiber.lambda0" or "fiber.f0" is important.
% -------------------------------------------------------------------------
c = 2.99792458e-4; % speed of ligth, m/ps

% Get lambda0 from input f0 or lambda0
if isfield(input_sim,'f0')
    default_sim.lambda0 = c/input_sim.f0;
else
    if isfield(input_sim,'lambda0')
        default_sim.lambda0 = input_sim.lambda0;
    else
        default_sim.lambda0 = 1030e-9;
    end
end

% -------------------------------------------------------------------------
% fiber
% -------------------------------------------------------------------------
% Basic properties
default_fiber.fiber = '1060XP';
default_fiber.material = 'fused silica'; % for finding the Raman response in GNLSE3D_propagate()
default_fiber.n2 = 2.3e-20; % m^2/W

% L0 is necessary to be put into "default_fiber" first for "gain_coeff" calculation.
if isfield(input_fiber,'L0')
    default_fiber.L0 = input_fiber.L0;
else
    default_fiber.L0 = 2; % m
end

% -------------------------------------------------------------------------
% sim
% -------------------------------------------------------------------------
% Basic settings
default_sim.f0 = c/default_sim.lambda0; % THz
default_sim.save_period = 0; % m

% Polarization modes
default_sim.scalar = true;

% Adaptive method
% From my test with 1-m-long wave breaking in a 1060XP fiber, the following
% numbers are required. If either of them is larger, the beam will not
% exhibit a LP01 spatial profile. In a 3D-GNLSE, resolving the spatial
% dimension correctly requires a super small step size.
default_sim.adaptive_deltaZ.DW_threshold = 1e-6; % the threshold of the dispersion/waveguide part in the adaptive method
default_sim.adaptive_deltaZ.threshold = 1e-6; % the threshold of the RK4IP part in the adaptive method

% Algorithms to use
default_sim.gpu_yes = true;
default_sim.Raman_model = 1; % isotropic Raman model
default_sim.gain_model = 0;

% Others
default_sim.pulse_centering = true; % center the pulse according to the time window
default_sim.Raman_sponRS = true; % include spontaneous Raman scattering
default_sim.num_photon_noise_per_band = 0; % don't include photon noise
default_sim.gpuDevice.Index = 1; % the gpuDevice to use
default_sim.progress_bar = true;
default_sim.progress_bar_name = '';
default_sim.cuda_dir_path = fullfile(upper_folder,'cuda');

%%
% =========================================================================
% Merge settings with the input, which have higher priorities than the
% default ones.
% =========================================================================
if isempty(input_fiber)
    fiber = default_fiber;
elseif isstruct(input_fiber)
    fiber = catstruct(default_fiber, input_fiber);
else
    error('LoadDefaultGNLSE3DPropagate:InputFiberError',...
            '"input_fiber" should be a "structure".');
end
if isempty(input_sim)
    sim = default_sim;
elseif isstruct(input_sim)
    sim = catstruct(default_sim, input_sim);
else
    error('LoadDefaultGNLSE3DPropagate:InputSimError',...
            '"input_sim" should be a "structure".');
end

end