function output = saturable_absorber_action_3d(input, saturation_power, moddepth, mode_profiles, dx, Aeff)
% SATURABLE_ABSORBER_ACTION_3D Apply an ideal spatial-temporal saturable absorber effect
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
%
% This function takes the spatial profiles of the modes into account when
% applying saturable absorption, which should be much more accurate for
% most real or effective saturable absorbers.

if isstruct(input)
    input_field = input.fields(:, :, end);
else
    input_field = input(:, :, end);
end

Nt = size(input_field, 1);
num_modes = size(input_field, 2);
Nx = size(mode_profiles, 1);

% Calculate the normalization constants
norms = zeros(num_modes, 1);
for midx = 1:num_modes
     norms(midx) = sqrt(sum(sum(abs(mode_profiles(:, :, midx)).^2)))*dx;
end

% This is rather brute force, but it works and it only needs to be
% calculated one per round trip. I would suggest downsampling the spatial
% mode profiles, however, or this calculation may take longer than the
% propagation.
output_field = zeros(Nt, num_modes);
for ii = 1:Nt
    % First, build the full field at time t_i
    full_field_i = zeros(Nx, Nx);
    for midx = 1:num_modes
        full_field_i = full_field_i + mode_profiles(:, :, midx)/norms(midx)*input_field(ii, midx);
    end
    
    % Then apply 3D ideal saturable absorbtion
    %full_field_i = full_field_i.*sqrt(1 - moddepth./(1+abs(full_field_i).^2*Aeff/saturation_power));
    full_field_i = full_field_i.*sqrt(1 - moddepth.*exp(-abs(full_field_i).^2*Aeff/saturation_power));
    
    % Then project the field back into the modes
    for midx = 1:num_modes
        output_field(ii, midx) = sum(sum(mode_profiles(:, :, midx).*full_field_i/norms(midx)))*dx^2;
    end
end

if isstruct(input)
    output = input;
    output.fields = output_field;
else
    output = output_field;
end

end