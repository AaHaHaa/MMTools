function output = saturable_absorber_action_3d_modified(input, saturation_power, moddepth, mode_profiles, dx, Aeff, sim)
%SATURABLE_ABSORBER_ACTION_3D_MODIFIED Apply an ideal spatial-temporal saturable absorber effect
%
% input - either a struct or a matrix:
%   input.fields - a (N, num_modes, m) matrix with the fields for each mode, in the time domain
%  OR
%   input - a (N, num_mdoes, m) matrix with the fields for each mode, in the time domain
%
% saturation_power - the scale power for the saturable absorber, in W
% moddepth - the modulation depth, 0-1
% mode_profiles - a (Nx, Nx, num_modes) matrix with the spatial profile of each mode
% dx - the spatial step, in m
% Aeff - in m^2
% sim - "sim.gpu = true or false" determines whether to use GPU or not
%
% This function takes the spatial profiles of the modes into account when
% applying saturable absorption, which should be much more accurate for
% most real or effective saturable absorbers.

if isstruct(input)
    input_field = input.fields(:, :, end);
else
    input_field = input(:, :, end);
end

if exist('sim','var')
    if sim.single_yes
        input_field = single(input_field);
        mode_profiles = single(mode_profiles);
    end
    if sim.gpu_yes
        input_field = gpuArray(input_field);
        mode_profiles = gpuArray(mode_profiles);
    end
end

Nt = size(input_field,1);
Nx = size(mode_profiles,1);
num_modes = size(input_field, 2);

% Calculate the normalization constants
if exist('sim','var')
    if sim.single_yes
        if sim.gpu_yes
            trans_field = zeros(Nt,num_modes,'single','gpuArray');
        else
            trans_field = zeros(Nt,num_modes,'single');
        end
    else
        if sim.gpu_yes
            trans_field = zeros(Nt,num_modes,'gpuArray');
        else
            trans_field = zeros(Nt,num_modes);
            %reflect_field = zeros(Nt,num_modes);
        end
    end
else
    trans_field = zeros(Nt,num_modes);
    %reflect_field = zeros(Nt,num_modes);
end
norms = sqrt(sum(sum( abs(mode_profiles).^2 ,1),2))*dx;

% This is rather brutal force, but it works and it only needs to be
% calculated one per round trip. I would suggest downsampling the spatial
% mode profiles, however, or this calculation may take longer than the
% propagation.

normalized_mode_profiles_xym = bsxfun(@rdivide, mode_profiles,norms); % Normalization

clear mode_profiles

% Calculate the result with a few time periods instead of finishing it at once
if exist('sim','var')
    parts = memory_control(Nt,Nx,sim);
else
    parts = memory_control(Nt,Nx);
end

for i = 1:length(parts)
    part_i = parts{i};
    input_field_part = input_field(part_i,:);

    full_field_txy = recompose_into_space(normalized_mode_profiles_xym, input_field_part);

    % Apply 3D ideal saturable absorption
    %SA_trans_field_txy = full_field_txy.*sqrt(1 - moddepth./(1+abs(full_field_txy).^2/(saturation_power/Aeff)));
    SA_trans_field_txy = full_field_txy.*sqrt(1 - moddepth.*exp(-abs(full_field_txy).^2/(saturation_power/Aeff)));
    %SA_reflect_field_txy = full_field_txy.*sqrt(moddepth./(1+abs(full_field_txy).^2/(saturation_power/Aeff)));
    %SA_reflect_field_txy = full_field_txy.*sqrt(moddepth.*exp(-abs(full_field_txy).^2/(saturation_power/Aeff)));

    % Project the field back into the modes
    trans_field(part_i,:) = decompose_into_modes(normalized_mode_profiles_xym, SA_trans_field_txy, dx);
    %reflect_field(split_i,:) = decompose_into_modes(normalized_mode_profiles_mxy, SA_reflect_field_txy, dx);
end

%output = struct('trans_field',trans_field,'reflect_field',reflect_field);
if exist('sim','var') && sim.gpu_yes
    output = gather(trans_field);
else
    output = trans_field;
end

end

%% memory_control
function parts = memory_control(Nt,Nx,sim)
% MEMORY_CONTROL Split the large (t,x,y) matrices w.r.t the 1st dimension(t)
% in the calculation above to prevent memory overload because of a large
% amount of spatial points used in mode profiles of each spatial mode.

if exist('sim','var') && sim.gpu_yes
    gd = gpuDevice;
    memory_limit = gd.AvailableMemory;
    %memory_limit = 5*2^(30); % 5GB
else
    memory_limit = 30*2^(30); % 30GB
end

% "single" precision is assumed here.
%    single: 4 bytes
%    double: 8 bytes
if exist('sim','var') && sim.gpu_yes
    if sim.single_yes
        precision = 4;
    else
        precision = 8;
    end
else
    precision = 4;
end
mem_complex_number = precision*2;

num_parts_matrix = ceil(mem_complex_number*Nt*Nx^2/memory_limit)*4; % Choose 4 because 3 would blow up the GPU memory on my desktop, but there should be a stricter calculation of this.
if num_parts_matrix == 1
    parts = {1:Nt};
else
    num_each_split = floor(Nt/num_parts_matrix);
    num_parts_matrix = ceil(Nt/num_each_split);
    index = [1:Nt zeros(1,num_parts_matrix*num_each_split-Nt)];
    tmp = reshape(index,num_each_split,num_parts_matrix);
    parts = cell(num_parts_matrix,1);
    for i = 1:num_parts_matrix
        parts{i} = tmp(:,i);
    end
    parts{end} = parts{end}(parts{end}~=0);
end

end