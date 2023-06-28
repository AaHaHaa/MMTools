function output = saturable_absorber_action_simple_modified(input, saturation_power, moddepth)
%SATURABLE_ABSORBER_ACTION_SIMPLE_MODIFIED Apply an ideal saturable absorber effect to each mode individually
%
% input - either a struct or a matrix:
%   input.fields - a (N, num_modes, m) matrix with the fields for each mode, in the time domain
%  OR
%   input - a (N, num_mdoes, m) matrix with the fields for each mode, in the time domain
%
% saturation_power - the scale power for the saturable absorber, in W
% moddepth - the modulation depth, 0-1

if isstruct(input)
    input_field = input.fields(:, :, end);
    output = input;
else
    input_field = input(:, :, end);
end

powers = sum(abs( input_field ).^2,2);
%trans_field = input_field.*sqrt(1 - moddepth./(1 + powers/saturation_power));
%reflect_field = input_field.*sqrt(moddepth./(1 + powers/saturation_power));
trans_field = input_field.*sqrt(1 - moddepth*exp(-powers/saturation_power));
reflect_field = input_field.*sqrt(moddepth*exp(-powers/saturation_power));

output.trans_field = trans_field;
output.reflect_field = reflect_field;

end