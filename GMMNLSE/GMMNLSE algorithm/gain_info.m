function gain_rate_eqn = gain_info( fiber,sim,gain_rate_eqn,lambda,mode_profiles )
%GAIN_INFO Computes several information related to the gain medium.
%
% =========================================================================
% =============== Call this function with the following code ==============
% =========================================================================
% f = ifftshift( (-N/2:N/2-1)/N/dt + sim.f0 ); % in the order of "omegas" in the GMMNLSE_propagate()
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
%           absorption_wavelength_to_get_N_total - nm
%           absorption_to_get_N_total - dB/m
%           cross_section_filename - the filename for the doped ion cross section data
%                                    Currently I have 'Liekki Yb_AV_20160530.txt' for Yb and 'optiwave Er' for Er
%                                    Er data is from optiwave website: https://optiwave.com/resources/applications-resources/optical-system-edfa-basic-concepts/
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
%           tau - the lifetime of the upper states, which is used in spontaneous emission; in "s".
%                 (1) lifetime of Yb in F_(5/2) state = 840 us (Paschotta et al., "Lifetme quenching in Yb-doped fibers")
%                 (2) lifetime of Er in (^4)I_(13/2) state = 8-10 ms (Z. Y. Zhang et al., "Fluorescence decay-time characteristics of erbium-doped optical fiber at elevated temperatures")
%           t_rep - the roundtrip time (1/repetition rate) of the pulse, which is used to calculate the power of the "signal" pulse; in "s"
%
%       rate equation model algorithm info -->
%
%           export_N2 - 1(true) or 0(false); Whether to export N2, the ion density of the upper state, or not
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
if length(sim.midx) == 1 % single mode
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
addpath([upper_folder 'Gain_rate_eqn/'],[upper_folder 'Gain_rate_eqn/gain cross sections/']);

%% Add a more-verbose parameter
gain_rate_eqn.include_ASE = ~gain_rate_eqn.ignore_ASE;
gain_rate_eqn.load_profiles = (isfield(fiber,'MM_folder') && ~isempty(fiber.MM_folder));

%% linear-oscillator model
if gain_rate_eqn.linear_oscillator
    addpath([upper_folder 'Gain_rate_eqn/linear oscillator/']);

    gain_rate_eqn.reuse_data = true; % Force to reuse the previously calculated data because it's a linear oscillator
end

%% Read information based on the gain medium to use
gain_rate_eqn = gain_medium(gain_rate_eqn);

%% Cross sections
% "lambda" must be a column vector.
if size(lambda,1) == 1
    lambda = lambda.';
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
% Read cross sections from the file.
necessary_lambda = [gain_rate_eqn.absorption_wavelength_to_get_N_total*1e-9; ...
                    gain_rate_eqn.pump_wavelength*1e-9; ...
                    lambda*1e-9];
switch gain_rate_eqn.gain_medium
    case 'Yb'
        [absorption,emission] = read_cross_sections_Yb(gain_rate_eqn.cross_section_filename,necessary_lambda); % read the file
        absorption = absorption*1e12; % change the unit to um^2
        emission = emission*1e12;
        
        cross_sections_pump = struct('absorption',absorption(2),'emission',emission(2)); % pump
        cross_sections = struct('absorption',absorption(3:end),'emission',emission(3:end)); % signal, ASE
    case {'Er','Nd'} % consider excited-state absorption
        [absorption,emission,ESA] = read_cross_sections_Er_Nd(gain_rate_eqn.cross_section_filename,necessary_lambda); % read the file
        absorption = absorption*1e12; % change the unit to um^2
        emission = emission*1e12;
        ESA = ESA*1e12;
        
        cross_sections_pump = struct('absorption',absorption(2),'emission',emission(2),'ESA',ESA(2)); % pump
        cross_sections = struct('absorption',absorption(3:end),'emission',emission(3:end),'ESA',ESA(3:end)); % signal, ASE
end

cross_sections = structfun(@(x) permute(x,[2 3 4 5 1]),cross_sections,'UniformOutput',false); % change it to the size (1,1,1,1,N)

%% Overlap factor of the field and the dopping area
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
% For small-signal absorption, N2~0 and N1~N_total.
% pump is proportional to "exp(-integral2(overlap_factor.pump)*N_total*absorption_cross_section*L)", L: propagation length
% absorption_dB/m =10*log10( exp(-integral2(overlap_factor.pump)*N_total*absorption_cross_section*(L=1m)) )
N_total = log(10^(gain_rate_eqn.absorption_to_get_N_total/10))./(((gain_rate_eqn.core_diameter/gain_rate_eqn.cladding_diameter)^2)*absorption(1)*1e6); % doped ion density based on absorption at a specific wavelength; in "1/um^3"
if gain_rate_eqn.load_profiles && ... % not fundamental-mode of a step-index fiber
       isequal(sim.step_method,'MPA') % multimode
    N_total = N_total*core_region; % size: (Nx,Nx)
end

%% Check the validity of the code
% Because I use the approximation, sqrt(1+x)=1+x/2 if x is small, in
% calculating signal fields with MPA, the code will give error here if
% this approximation is bad.
if isequal(sim.step_method,'MPA') && ... % multimode with MPA
        (gain_rate_eqn.include_ASE || gain_rate_eqn.reuse_data || gain_rate_eqn.linear_oscillator) % with a fixed-step method (so the appropriateness of deltaZ needs to be checked here)
    relative_N2 = 0:0.1:0.5;
    relative_gain = relative_N2.*cross_sections.emission - (1-relative_N2).*cross_sections.absorption;
    relative_gain(relative_gain<0) = 0; % I think (?) it only neds to resolve the "gain" correctly
    
    tol_approximation = 1e-4; % I've found that 1e-3 is not enough
    approx_error = @(x)abs((sqrt(1+x)-(1+x/2))./sqrt(1+x));
    if approx_error( 2*(sim.deltaZ/sim.MPA.M*1e6)*max(N_total(:))*max(relative_gain(:)) ) > tol_approximation
        error('gain_info:deltaZError',...
              'The deltaZ is too large for this code to run because of the approximation of "sqrt(1+x)=1+x/2" I use for calculating the gain for multimode cases.');
    end
end

% This code assumes the population inversion reaches the steady state
% because of high repetition rate of pulses.
% I'll check the highest recovery lifetime for the inversion. This should
% be larger than repetition rate for this assumption to be true.
h = 6.626e-34; % J*s
c = 299792458; % m/s
tc = 1/(max(overlap_factor.pump(:))*(gain_rate_eqn.pump_wavelength*1e-9)/(h*c)*(cross_sections_pump.absorption+cross_sections_pump.emission)*max(gain_rate_eqn.copump_power,gain_rate_eqn.counterpump_power)+1/gain_rate_eqn.tau);
if tc < gain_rate_eqn.t_rep
    warning('The repetition rate isn''t high enough for this code to be accurate.');
end

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
if gain_rate_eqn.load_profiles && ... % there are loaded mode profiles for multimode, higher-order modes, or user-defined modes
       isequal(sim.step_method,'MPA') % multimode with MPA
    trapz2 = @(x) trapz(trapz(x,1),2)*gain_rate_eqn.mode_profile_dx^2; % take the integral w.r.t. the x-y plane
    FmFnN = trapz2(overlap_factor.signal.*N_total);
    GammaN = trapz2(overlap_factor.pump.*N_total);
else
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
            error('OSError:GMMNLSE_gain_rate_eqn',...
                  'iOS is not well supported yet.');
        end
    end
end

%% Saved data
if ~isfield(gain_rate_eqn,'saved_data')
    gain_rate_eqn.saved_data = [];
end

%% Put all information into "gain_rate_eqn"
gain_rate_eqn.cross_sections_pump = cross_sections_pump;
gain_rate_eqn.cross_sections      = cross_sections;
gain_rate_eqn.overlap_factor      = overlap_factor;
gain_rate_eqn.N_total             = N_total;
gain_rate_eqn.FmFnN               = FmFnN;
gain_rate_eqn.GammaN              = GammaN;

end