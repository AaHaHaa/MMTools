function foutput = UPPE3D_propagate(fiber, initial_condition, sim)
%UPPE3D_PROPAGATE Propagate an initial full-field pulse through an arbitrary 
% distance of a nonlinear medium, such as an optical fiber
%   This is a caller function, calling
%   UPPE_propagate_xy() or
%   UPPE_propagate_r()
%   based on whether the field is radially symmetric.

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

%% Determine whether to use the radially-symmetric scheme
rxy_str = 'xy';
if isfield(initial_condition,'r')
    rxy_str = 'r';
end

UPPE3D_propgation_func = str2func(['UPPE3D_propagate_', rxy_str]);

%% Run the pulse propagation
foutput = UPPE3D_propgation_func(fiber, initial_condition, sim);

end