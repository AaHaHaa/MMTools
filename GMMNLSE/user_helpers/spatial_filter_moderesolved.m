function output = spatial_filter_moderesolved(input, filter_matrix)
% SPATIAL_FILTER_MODERESOLVED Apply the spatial filter using a pre-calculated filter matrix
% input - a (N, num_modes) matrix with the time profile of each mode
% filter_matrix -  a (num_modes, num_modes) matrix that couples the modes through the filter

if isstruct(input)
    output = input;
    output.fields = (filter_matrix*(input.fields(:, :, end)'))';
else
    output = (filter_matrix*(input(:, :, end)'))';
end

end