function gain_rate_eqn = gain_info( fiber,sim,gain_rate_eqn,lambda,mode_profiles )
%GAIN_INFO Computes several information related to the gain medium.
%
% =========================================================================
% =============== Call this function with the following code ==============
% =========================================================================
% f = ifftshift( (-N/2:N/2-1)/N/dt + sim.f0 ); % offset angular frequency; in the order of "Omega" in the GMMNLSE_propagate()
% c = 299792.458; % nm/ps
% lambda = c./f; % nm
%
% gain_rate_eqn = gain_info( fiber,sim,gain_rate_eqn,lambda );
% output_field = GMMNLSE_propagate(fiber,input_field,sim,gain_rate_eqn);
%
% =========================================================================
%
%   gain_rate_eqn:
%
%       multimode mode-profile folder -->
%
%           MM_folder - a string; where the betas.mat and S_tensor_?modes.mat are
%
%       oscillator info -->
%
%           reuse_data - true(1) or false(0);
%                        For a ring or linear cavity, the pulse will enter a steady state eventually.
%                        If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
%
%           linear_oscillator - false: not a linear oscillator
%                               true:  a linear oscillator
%                           For a linear oscillator, there are pulses from both directions simultaneously, which will both contribute to saturating the gain;
%                           therefore , the backward-propagating pulses need to be taken into account.
%
%                           How to use it:
%                               gain_rate_eqn.saved_data = rate_gain_saved_data;
%                               prop_output = GMMNLSE_propagate(fiber, input_field, sim, gain_rate_eqn);
%                               (next)rate_gain_saved_data = prop_output.saved_data;
%
%                               "rate_gain_saved_data" contains
%                               	signal_fields, signal_fields_backward,
%                               	Power_pump_forward, Power_pump_backward,
%                                   Power_ASE_forward, Power_ASE_backward
%
%       fiber info -->
%
%           core_diameter - um; where the doped ion stays
%           cladding_diameter - um
%           core_NA - numerical aperture of the gain fiber
%
%       doped ion info -->
%
%           gain_medium - doped ion name
%           absorption_wavelength_to_get_N_total - nm
%           absorption_to_get_N_total - dB/m
%
%       pump info -->
%
%           pump_wavelength - nm
%           pump_direction - "Don't set this parameter" since it'll be automatically determined from "copump_power" and "counterpump_power" below.
%                            'co', 'counter', or 'bi'
%                            'co' reads only the copump_power, 'counter' reads only the counterpump_power, while 'bi' reads both.
%           copump_power - W
%           counterpump_power - W
%       
%       mode profiles inside the gain core area -->
%
%           downsampling_factor - reduce the size of mode profiles for multimode
%
%       computational info -->
%
%           t_rep - the roundtrip time (1/repetition rate) of the pulse, which is used to calculate the power of the "signal" pulse; in "s"
%
%       ASE info -->
%
%           sponASE_spatial_modes - the number of ASE's supported spatial modes
%                                   In LMA fibers, the number of ASE modes can be larger than one as the signal field, so this factor is used to correctly considered ASE.
%                                   If empty like [], it's length(sim.midx).
%
%       rate equation model algorithm info -->
%
%           ignore_ASE - 1(true) or 0(false)
%           max_iterations - the maximum number of iterations
%           tol - the tolerance of this iteration loop. If the difference between the last two results is smaller than this tolerance, we're done. 
%           verbose - show the information(final pulse energy) during iterations
%           memory_limit: the RAM limit for this computation
%                         (This can be found by default.)
%
%   lambda - the wavelengths of the computational region; in "nm"
%
%   mode_profiles (If not specified, this code will read them from "fiber.MM_folder"):
%
%       mode_profiles - (Nx,Nx,num_spatial_modes); the eigenmode field profiles of modes;
%                       It'll be normalized into the unit of "1/um"
%       mode_profiles_x - the x-position for mode_profiles; in "um"
%
%   =======================================================================
%   Output (in "gain_rate_eqn"):
%
%       cross_sections_pump - um^2
%       cross_sections - um^2
%       overlap_factor - 1/um^2
%       N_total - a scalar for single mode and (Nx,Nx) for multimode; the doped ion density; in "1/um^3"
%       FmFnN - precalculate the integral2(overlap_factor*N_total) for the signal and ASE
%       GammaN - precalculate the integral2(overlap_factor*N_total) for the pump

%% Pick the stepping algorithm to use depending on the number of modes
if isscalar(sim.midx) % single mode
    sim.step_method = 'RK4IP';
else % multimode
    sim.step_method = 'MPA';
end
% Only MPA is implemented in Random mode coupling
if sim.rmc.model
    sim.step_method = 'MPA';
end

%% Add the folder of functions of gain-rate-equation model and its functions
% Besides loading mode coupling folder, this "sep_char" is also used in GPU setup below.
if ispc
    sep_char = '\';
else % unix
    sep_char = '/';
end
current_path = mfilename('fullpath');
sep_pos = strfind(current_path,sep_char);
upper_folder = current_path(1:sep_pos(end-1));
addpath([upper_folder 'Steady_state_Gain_rate_eqn/'],...
        [upper_folder 'Steady_state_Gain_rate_eqn/gain cross sections/'],...
        [upper_folder 'Steady_state_Gain_rate_eqn/Judd-Ofelt theory/Material data/']);

%% Add more-verbose parameters
gain_rate_eqn.include_ASE = ~gain_rate_eqn.ignore_ASE;
gain_rate_eqn.load_profiles = (isfield(fiber,'MM_folder') && ~isempty(fiber.MM_folder));

%% linear-oscillator model
if gain_rate_eqn.linear_oscillator
    addpath([upper_folder 'Steady_state_Gain_rate_eqn/linear oscillator/']);

    gain_rate_eqn.reuse_data = true; % Force to reuse the previously calculated data because it's a linear oscillator
end

%% "lambda" error check
% "lambda" must be a column vector.
if size(lambda,1) == 1
    lambda = lambda.';
end
if any(lambda(:)<0)
    error('gain_info:lambdaError',...
          ['Wavelength, the "lambda" input variable, cannot be negative.\n',...
           'If possible, don''t use the pulse frequency as the center of the frequency window, which helps offset the frequency window from getting close to zero.\n',...
           'Use find_tw_f0() to help do this.']);
end
% MATLAB before 2017 doesn't have "issortedrows()" and 'monotonic' argument in "issorted()"
MATLAB_version = version('-release'); MATLAB_version = str2double(MATLAB_version(1:4));
if MATLAB_version < 2017
    do_ifftshift_on_lambda = (issorted(lambda) || issorted(flipud(lambda)));
else
    do_ifftshift_on_lambda = issortedrows(lambda,'monotonic'); % = issorted(lambda,'monotonic');
end
if do_ifftshift_on_lambda
    lambda = ifftshift(lambda,1);
end

%% Narrowband transformation due to the scaled Fourier transform
% lambda needs to be re-defined
if sim.cs.cs > 1
    lambda0 = lambda; % save the original lambda for extracting cross sections in the original lambda
    lambda = lambda(1:sim.cs.cs:end);
end

%% Function container to read necessary parameters based on the gain medium to use later
gain_func = gain_medium();

% Read necessary parameters and cross sections based on the gain medium to use
gain_rate_eqn = gain_func.load_medium_parameters(gain_rate_eqn);
[gain_rate_eqn,...
 cross_sections,cross_sections_pump,...
 GSA_find_Ntotal] = gain_func.load_cross_sections(gain_rate_eqn,lambda);

% Compute the original cross sections without the narrowband transformation (scaled Fourier transform)
if sim.cs.cs > 1
    [~,cross_sections0] = gain_func.load_cross_sections(gain_rate_eqn,lambda0);
end

%% Overlap factor of the field and the doping area
% Load mode profiles for multimode
if gain_rate_eqn.load_profiles % ready to load mode profiles
    lambda0 = int16(lambda(1));
    load(sprintf('%smode%uwavelength%u.mat',fiber.MM_folder,sim.midx(1),lambda0*10),'phi','x'); % load 1st mode first to get the size and dimension of the mode profile
    mode_profiles.x = x; % x-position vector
    mode_profiles.mode_profiles = zeros(length(x),length(x),length(sim.midx)); % initialization
    mode_profiles.mode_profiles(:,:,1) = phi; % the 1st mode
    for ni = 2:length(sim.midx)
        n = sim.midx(ni);
        load(sprintf('%smode%uwavelength%u.mat',fiber.MM_folder,n,lambda0*10),'phi');
        mode_profiles.mode_profiles(:,:,ni) = phi;
    end
    
    gain_rate_eqn.mode_profile_dx = abs(mode_profiles.x(2)-mode_profiles.x(1)); % unit: um; This variable will be used later in "solve_gain_rate_eqn", so I put it in "gain_rate_eqn".
    % Normalization
    norms = sqrt(sum( abs(mode_profiles.mode_profiles).^2 ,[1,2]))*gain_rate_eqn.mode_profile_dx;
    mode_profiles.mode_profiles = mode_profiles.mode_profiles./norms; % unit: 1/um
    
    % Choose only the core region for the mode profiles to save the memory.
    % Also downsample the spatial profiles.
    chosen_region_left_idx = find(mode_profiles.x > -gain_rate_eqn.core_diameter/2,1) - 1;
    chosen_region_right_idx = find(mode_profiles.x < gain_rate_eqn.core_diameter/2,1,'last') + 1;
    core_region_idx = chosen_region_left_idx:chosen_region_right_idx;
    mode_profiles.x = downsample(mode_profiles.x(core_region_idx),gain_rate_eqn.downsampling_factor);
    old_large_mode_profiles = mode_profiles.mode_profiles;
    mode_profiles.mode_profiles = zeros(length(mode_profiles.x),length(mode_profiles.x),length(sim.midx));
    for ni = 1:length(sim.midx)
        mode_profiles.mode_profiles(:,:,ni) = downsample(downsample(old_large_mode_profiles(core_region_idx,core_region_idx,ni), gain_rate_eqn.downsampling_factor)', gain_rate_eqn.downsampling_factor)'; % downsample;;
    end
    
    gain_rate_eqn.mode_profile_dx = gain_rate_eqn.mode_profile_dx*gain_rate_eqn.downsampling_factor;
    
    % Core region where active ions lie
    [x,y] = meshgrid(mode_profiles.x,mode_profiles.x);
    core_region = (x.^2 + y.^2) <= (gain_rate_eqn.core_diameter/2)^2;
    
    % Consider the core region only
    mode_profiles.mode_profiles = mode_profiles.mode_profiles.*core_region;
end

% pump overlap factor = A_doping/A_cladding
if isequal(sim.step_method,'RK4IP') % single spatial mode
    overlap_factor.pump = 1/(pi*(gain_rate_eqn.cladding_diameter/2)^2);
else % multimode
    overlap_factor.pump = 1/(pi*(gain_rate_eqn.cladding_diameter/2)^2)*core_region;
end

% signal overlap factor = F_m*conj(F_n)
if ~gain_rate_eqn.load_profiles % fundamental mode of a step-index fiber
    V = 2*pi*(gain_rate_eqn.core_diameter/2)/(lambda(1)*1e-3)*gain_rate_eqn.core_NA; % V-number, normalized frequency
    if V < 0.8 || V > 2.8
        warning('For the computation of the fundamental mode, I use "Whitley''s Gaussian-mode approximation". It works only in the range of V=0.8-2.8.\nCurrent V number is %4.2f',V);
    end
    % w_over_a = 0.65+1.619*V^(-3/2)+2.879*V^(-6); % Marcuse et al., "Loss analysis of single-mode fiber splices" (1977)
    w_over_a = 0.616+1.66*V^(-3/2)+0.987*V^(-6); % Whitley et al., "Alternative Gaussian spot size polynomial for use with doped fiber amplifiers" (1993)
    overlap_factor.signal = ( 1-exp(-2/w_over_a^2) )/(pi*(gain_rate_eqn.core_diameter/2)^2); % 1/um^2
else % multimode or higher-order modes
    overlap_factor.signal = mode_profiles.mode_profiles.*permute(conj(mode_profiles.mode_profiles),[1 2 4 3]);
    
    if isequal(sim.step_method,'RK4IP') % single higher-order mode
        overlap_factor.signal = sum(overlap_factor.signal(:))*gain_rate_eqn.mode_profile_dx^2/(pi*(gain_rate_eqn.core_diameter/2)^2); % 1/um^2
    end
end

%% Doped ion density
% For small-signal absorption, N1~0 and N0~N_total.
% pump is proportional to "exp(-integral2(overlap_factor.pump)*N_total*absorption_cross_section*L)", L: propagation length
% absorption_dB/m =10*log10( exp(-integral2(overlap_factor.pump)*N_total*absorption_cross_section*(L=1m)) )
N_total = log(10^(gain_rate_eqn.absorption_to_get_N_total/10))./(((gain_rate_eqn.core_diameter/gain_rate_eqn.cladding_diameter)^2)*GSA_find_Ntotal*1e6); % doped ion density based on absorption at a specific wavelength; in "1/um^3"
if gain_rate_eqn.load_profiles && ... % not fundamental-mode of a step-index fiber
       isequal(sim.step_method,'MPA') % multimode
    N_total = N_total*core_region; % size: (Nx,Nx)
end

%% Read necessary parameters based on the gain medium to use
gain_rate_eqn = gain_func.load_N_related_parameters(gain_rate_eqn,N_total);

%% Check the validity of the code
% Number of ASE spatial modes
% For LMA fibers, there are more than one spatial modes for ASE although
% the signal field mostly stays only within the fundamental mode. In this
% situation, the simulations are mostly run with single mode, so
% consideration of multimode ASE is included with this 
% sponASE_spatial_modes factor.
if gain_rate_eqn.include_ASE
    if isempty(gain_rate_eqn.sponASE_spatial_modes)
        gain_rate_eqn.sponASE_spatial_modes = length(sim.midx);
    else
        if gain_rate_eqn.sponASE_spatial_modes < length(sim.midx)
            error('gain_info:NumASEError',...
                  'The number of ASE spatial modes need to be larger than that of the signal field.');
        end
    end
end

%% Pre-compute the integral of "overlap_factor*N_total"
% For coupled equtions with more than or equal to two levels,
% pre-computations aren't necessary.
if length(gain_rate_eqn.energy_levels) == 2 % two-level system
    if gain_rate_eqn.load_profiles && ... % there are loaded mode profiles for multimode, higher-order modes, or user-defined modes
           isequal(sim.step_method,'MPA') % multimode with MPA
        trapz2 = @(x) trapz(trapz(x,1),2)*gain_rate_eqn.mode_profile_dx^2; % take the integral w.r.t. the x-y plane
        FmFnN = trapz2(overlap_factor.signal.*N_total);
        GammaN = trapz2(overlap_factor.pump.*N_total);
    else
        FmFnN = [];
        GammaN = [];
    end
else % multi-level system
    FmFnN = [];
    GammaN = [];
end


%% Query the gpuDevice
if sim.gpu_yes
    if ~isfield(sim,'gpuDevice') || ~isfield(sim.gpuDevice,'Device')
        sim.gpuDevice.Device = gpuDevice(sim.gpuDevice.Index); % use the specified GPU device
    end
end

%% Find the memory limit
if ~isfield(gain_rate_eqn,'memory_limit')
    if sim.gpu_yes
        gain_rate_eqn.memory_limit = sim.gpuDevice.Device.AvailableMemory/2;
    else
        if ispc % Windows
            userview = memory;
            gain_rate_eqn.memory_limit = userview.MaxPossibleArrayBytes/2; % B
        elseif isunix % Unix, linux
            [~,w] = unix('free -b | grep Mem'); % Display the memory info in Bytes
            stats = str2double(regexp(w, '[0-9]*', 'match'));
            %memsize = stats(1)/1e6;
            freemem = stats(end); % B; availabel memory
            gain_rate_eqn.memory_limit = freemem/2;
        else % iOS
            error('gain_info:OSError',...
                  'iOS is not well supported yet.');
        end
    end
end

%% Saved data
if ~isfield(gain_rate_eqn,'saved_data')
    gain_rate_eqn.saved_data = [];
end

%% Reorganize cross sections into an array for faster numerical computations
% The 8th dimension is to store various cross-section data
% Dimension: (1,1,1,1,Nt,1,1,num_cross_sections),
% which each stands for, in gain computations, (Nx,Nx,num_spatial_modes,num_spatial_modes,Nt,M,num_polarization,num_cross_sections), where M: parallelization in MPA
cross_sections_pump = permute(cell2mat(struct2cell(cross_sections_pump)),[8,2,3,4,5,6,7,1]);
cross_sections = permute(cell2mat(struct2cell(cross_sections)),[8,2,3,4,5,6,7,1]);

if sim.cs.cs > 1
    cross_sections0 = permute(cell2mat(struct2cell(cross_sections0)),[8,2,3,4,5,6,7,1]);
end

%% Put all information into "gain_rate_eqn"
gain_rate_eqn.cross_sections_pump = cross_sections_pump;
gain_rate_eqn.cross_sections      = cross_sections;
gain_rate_eqn.overlap_factor      = overlap_factor;
gain_rate_eqn.N_total             = N_total;
gain_rate_eqn.FmFnN               = FmFnN;
gain_rate_eqn.GammaN              = GammaN;

% Save original cross sections without the narrowband transformation (scaled Fourier transform)
if sim.cs.cs > 1
    gain_rate_eqn.cross_sections0 = cross_sections0;
end

end