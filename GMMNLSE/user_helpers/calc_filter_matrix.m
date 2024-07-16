function filter_matrix = calc_filter_matrix(mode_profiles, x, spatial_hwhm, offcenter, filter_type)
%CALC_FILTER_MATRIX Calculate the matrix that applies a spatial filter to the field
% mode_profiles - a (Nx, Nx, num_modes) matrix with the spatial profile of each mode
% x - a Nx vector with the spatial coordinates in m
% spatial_hwhm - the half-width at half max of the filter, i.e. the radial extent at half max, in m
% offcenter - a (2,1) or (1,2) vector indicating offcenter values in x and y directions, respectively
% filter_type - the type of the filter, Gaussian or pinhole(a binary hard filter)
%
% The filter matrix T is the matrix such that A' = TA, where A is the
% vector of mode field envelopes before the filter and A' is the vector
% after the filter.
%
% The matrix is calculated by Tij = integral{Fi*Fj*F dxdy},
% where Fi and Fj are the spatial fields of modes i and j,
% and F is the profile of the filter

num_modes = size(mode_profiles, 3);

% Calculate the normalization constants for both sets of modes
norms_profiles = sqrt(sum(sum(abs(mode_profiles).^2,1),2));
mode_profiles = mode_profiles./norms_profiles;

% Setup the spatial grid
[X, Y] = meshgrid(x, x);

switch filter_type
    case 'Gaussian'
        spatial_w0 = spatial_hwhm/(2*sqrt(log(2))); % 2*sqrt(log(2))=1.665
        filter_profile = exp(-((X-offcenter(1)).^2+(Y-offcenter(2)).^2)/(2*spatial_w0^2));
    case 'pinhole'
        P = sqrt((X-offcenter(1)).^2+(Y-offcenter(2)).^2);
        filter_profile = (P<spatial_hwhm);
end

% Calculate the filter matrix
% Original code (easy to understand):
%   filter_matrix = zeros(num_modes, num_modes);
%   for midx1 = 1:num_modes
%       for midx2 = 1:num_modes
%            filter_matrix(midx1, midx2) = sum(sum(mode_profiles(:, :, midx1).*mode_profiles(:, :, midx2).*filter_profile));
%       end
%   end
% I vectorize it below.
[midx1,midx2] = ind2sub([num_modes,num_modes],1:num_modes^2);
filter_matrix = reshape(sum(sum(mode_profiles(:,:,midx1).*mode_profiles(:,:,midx2).*filter_profile,1),2),[num_modes,num_modes]);

end