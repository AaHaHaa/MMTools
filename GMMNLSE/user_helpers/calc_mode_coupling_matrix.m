function coupling_matrix = calc_mode_coupling_matrix(x1, mode_profiles1, x2, mode_profiles2)
%CALC_MODE_COUPLING_MATRIX Calculate the matrix that couples the fields
% from one fiber's mode basis to another fiber's mode basis.
% x1 - a Nx1 vector with the spatial coordinates for the original fiber
% mode_profiles1 - a (Nx1, Nx1, num_modes1) matrix with the spatial mode profiles for each mode of the original fiber
% x2 - a Nx2 vector with the spatial coordinates for the final fiber
% mode_profiles2 - a (Nx2, Nx2, num_modes2) matrix with the spatial mode profiles for each mode of the final fiber
%
% x1 and x2 should be in the same units, but they do not need to have the
% same extents or the same number of grid points. It is important, however,
% that the 0 point of x1 corresponds to the true center for fiber1, and the
% 0 point of x2 corresponds to the true center for fiber2.
%
% The order of 1 and 2 does matter, the resulting coupling matrix is
% defined such that A2 = T A1, i.e. T is a num_modes2xnum_modes1 matrix
% To go the other way, it will just be the transpose of T
%
% The matrix is calculated by Tij = integral{Fi/Ni*Fj/Nj dxdy}, where Fi
% and Fj are the spatial fields of modes i in the original fiber and j in
% the final fiber, and Ni and Nj are their normalization constants

num_modes1 = size(mode_profiles1, 3);
num_modes2 = size(mode_profiles2, 3);

% The fiber with the larger spatial extents will be used as the default,
% and the one with the smaller spatial extents will be interpolated
if max(x1) > max(x2)
    x = x1;
else
    x = x2;
end

Nx = length(x);

[X, Y] = meshgrid(x, x);
standard_mode_profiles1 = zeros(Nx, Nx, num_modes1);
standard_mode_profiles2 = zeros(Nx, Nx, num_modes2);

% Interpolate the fields with the smaller spatial extents
if max(x1) > max(x2)
    % In this case we keep the original modes the same...
    standard_mode_profiles1 = mode_profiles1;
    
    % And interpolate the final modes to fit
    [X2, Y2] = meshgrid(x2, x2);
    for midx = 1:num_modes2
        standard_mode_profiles2(:, :, midx) = interp2(X2, Y2, mode_profiles2(:, :, midx), X, Y);
    end
    standard_mode_profiles2(isnan(standard_mode_profiles2)) = 0;
else
    % In this case we keep the final modes the same...
    standard_mode_profiles2 = mode_profiles2;
    
    % And interpolate the original modes to fit
    [X1, Y1] = meshgrid(x1, x1);
    for midx = 1:num_modes1
        standard_mode_profiles1(:, :, midx) = interp2(X1, Y1, mode_profiles1(:, :, midx), X, Y);
    end
    standard_mode_profiles1(isnan(standard_mode_profiles1)) = 0;
end

% Calculate the normalization constants for both sets of modes
norms1 = zeros(num_modes1, 1);
for midx = 1:num_modes1
    norms1(midx) = sqrt(sum(sum(abs(standard_mode_profiles1(:, :, midx)).^2)));
end

norms2 = zeros(num_modes2, 1);
for midx = 1:num_modes2
    norms2(midx) = sqrt(sum(sum(abs(standard_mode_profiles2(:, :, midx)).^2)));
end

% Finally, calcaulte the coupling matrix
coupling_matrix = zeros(num_modes2, num_modes1);

for midx1 = 1:num_modes2
    for midx2 = 1:num_modes1
         coupling_matrix(midx1, midx2) = sum(sum(standard_mode_profiles2(:, :, midx1).*standard_mode_profiles1(:, :, midx2)))/ ...
                    (norms2(midx1)*norms1(midx2));
    end
end

end

