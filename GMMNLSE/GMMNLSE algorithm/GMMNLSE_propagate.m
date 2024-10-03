function foutput = GMMNLSE_propagate(fiber, initial_condition, sim, gain_rate_eqn)
%GMMNLSE_PROPAGATE Propagate an initial multimode pulse through an arbitrary distance of an optical fiber
%   This is a caller function, calling
%   GMMNLSE_propagate_with_adaptive() or
%   GMMNLSE_propagate_without_adaptive()
%   based on whether to use adaptive step-size method or not.
% -------------------------------------------------------------------------
%   It uses adaptive step-size method whenever possible, such as
%       all passive simulations, single-mode or multimode,
%       Gaussian-gain simulations,
%       rate-eqn-gain simulations without ASE, etc.
%   Only when random mode coupling, ASE, re-using rate-eqn data, or running
%   a linear-oscillator scheme, does it not use the adaptive step-size
%   method.

%%
if ispc
    sep_char = '\';
else % unix
    sep_char = '/';
end

% Load the folder
current_path = mfilename('fullpath');
sep_pos = strfind(current_path,sep_char);
upper_folder = current_path(1:sep_pos(end-1));

sim.cuda_dir_path = [upper_folder 'cuda'];

%% Random mode coupling
if sim.rmc.model
    if ~isfield(sim.rmc,'matrices') || ...
            size(sim.rmc.matrices,1) == 1 % There is only one mode, so there are no other modes for the linear mode coupling to happen
        sim.rmc.model = false;
    else
        % Add the folder regarding random mode coupling
        addpath([upper_folder 'Random_Mode_Coupling/'],...
                [upper_folder 'Random_Mode_Coupling/GMMNLSE algorithm/'],...
                [upper_folder 'Random_Mode_Coupling/expm/']);
    end
end

%% Gain model
% Check the gain model to use
if ~ismember(sim.gain_model,[0,1,2])
    error('GMMNLSE_propagate:gain_modelError',...
          '"sim.gain_model" can only be 0 (no gain), 1 (Gaussian gain), or 2 (rate-eqn gain).');
end

%% Determine whether to use adaptive-step-size method based on the gain model and the random mode coupling
adaptive_dz_str = 'with';
if sim.gain_model == 2 % rate-eqn-gain mdoel
    if gain_rate_eqn.include_ASE || gain_rate_eqn.reuse_data || gain_rate_eqn.linear_oscillator
        adaptive_dz_str = 'without';
        
        if sim.rmc.model
            error('GMMNLSE_propagate:rmc_modelError',...
                  'Random mode coupling does not support computations with ASE, reusing data, or of a linear oscillator.');
        end
    end
else % no gain or Gaussian-gain model
    gain_rate_eqn = [];
end
if sim.rmc.model % random mode coupling
    adaptive_dz_str = 'without';
    rmc_str = '_rmc';
else
    rmc_str = '';
end

GMMNLSE_propgation_func = str2func(['GMMNLSE_propagate_', adaptive_dz_str, '_adaptive', rmc_str]);

%% Run the pulse propagation
foutput = GMMNLSE_propgation_func(fiber, initial_condition, sim, gain_rate_eqn);

end