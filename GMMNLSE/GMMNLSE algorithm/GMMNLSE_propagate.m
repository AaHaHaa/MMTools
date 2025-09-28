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
if sim.gain_model == 2 % rate-eqn-gain model
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

%% Apply the narrowband transformation (due to the scaled Fourier transform)
scaledFT_func = narrowband_scaledFT();

if sim.cs.cs > 1
    Nt = size(initial_condition.fields,1);
    num_modes = size(initial_condition.fields,2);
    transformed_At = zeros(round(Nt/sim.cs.cs),num_modes); % I use round() here to prevent error. The error check for "cs" will be done in scaledFT_func.convert later, so we just let it pass here.
    for midx = 1:num_modes
        transformed_At(:,midx) = scaledFT_func.convert(initial_condition.fields(:,midx,end),sim.cs.cs);
    end

    initial_condition.dt = initial_condition.dt*sim.cs.cs;
    initial_condition.fields = transformed_At;
    
    % Narrowband transformation for the ASE powers
    if sim.gain_model == 2
        if gain_rate_eqn.include_ASE
            transformed_Power_ASE_forward = zeros(Nt/sim.cs.cs,num_modes);
            transformed_Power_ASE_backward = zeros(Nt/sim.cs.cs,num_modes);
            for midx = 1:num_modes
                tmp = fftshift(initial_condition.Power.ASE.forward(:,midx,end),1);
                transformed_Power_ASE_forward(:,midx) = ifftshift(tmp(1:sim.cs.cs:end)*sim.cs.cs,1);
    
                tmp = fftshift(initial_condition.Power.ASE.backward(:,midx,end),1);
                transformed_Power_ASE_backward(:,midx) = ifftshift(tmp(1:sim.cs.cs:end)*sim.cs.cs,1);
            end
            initial_condition.Power.ASE.forward = transformed_Power_ASE_forward;
            initial_condition.Power.ASE.backward = transformed_Power_ASE_backward;
        end
    end
end

%% Run the pulse propagation
foutput = GMMNLSE_propgation_func(fiber, initial_condition, sim, gain_rate_eqn);

%% Recover from the narrowband transformation
if sim.cs.cs > 1
    foutput.dt = foutput.dt/sim.cs.cs;
    foutput.fields = scaledFT_recover_field(scaledFT_func,foutput.fields,sim.cs.cs);

    if sim.gain_model == 2
        if gain_rate_eqn.include_ASE
            foutput.Power.ASE.forward = scaledFT_recover_power(foutput.Power.ASE.forward,sim.cs.cs);
            foutput.Power.ASE.backward = scaledFT_recover_power(foutput.Power.ASE.backward,sim.cs.cs);
        end
    end
end

end

%% Helper functions for recovering the field from the narrowband transformation due to the scaled Fourier transform
function rA = scaledFT_recover_field(scaledFT_func,A,cs)

Nt = size(A,1);
num_modes = size(A,2);
Nz = size(A,3);

rA = zeros(Nt*cs,num_modes,Nz);
for zi = 1:Nz
    for midx = 1:num_modes
        rA(:,midx,zi) = scaledFT_func.recover(A(:,midx,zi),cs);
    end
end

end

function rP = scaledFT_recover_power(P,cs)

Nt = size(P,1);
num_modes = size(P,2);
Nz = size(P,3);

recover_idx = linspace(1,Nt+1,Nt*cs+1)';
recover_idx = recover_idx(1:end-1); % remove the last "Nt+1" point
if mod(Nt,2) == 1
    recover_idx = [recover_idx(end-ceil((cs-1)/2)+1:end)-Nt;recover_idx(1:end-ceil((cs-1)/2))];
end

rP = zeros(Nt*cs,num_modes,Nz);
for zi = 1:Nz
    for midx = 1:num_modes
        rP(:,midx,zi) = interp1((1:Nt)',P(:,midx,zi),recover_idx,'linear','extrap')/cs;
    end
end

end