function [n, x, dx] = build_step(Nx, spatial_window, core_radius, clad_radius, n)
%BUILD_STEP Build a step refractive index profile
%
%   Nx - the number of spatial points in each dimension
%   spatial_window - the total size of space in each dimension, in um
%   radius - the radius of the step, in um
%   n - the refractive index of the core and cladding; a cell as {n_core,n_clad}

% refractive index
[n_core, n_clad, n_coat] = deal(n{:});

dx = spatial_window/Nx; % um

x = (-Nx/2:Nx/2-1)*dx;
[X, Y] = meshgrid(x, x);

%% Step index profile
% clad index
n = n_clad.*ones(1, Nx, Nx);

% core index
idx = permute(X.^2+Y.^2,[3,1,2]) <= core_radius^2;
n_core = repmat(n_core,1,sum(idx(:)));
n(repmat(idx ,size(n_core,1),1,1)) = n_core(:);

% coat index
idx = permute(X.^2+Y.^2,[3,1,2]) > clad_radius^2;
n_coat = repmat(n_coat,1,sum(idx(:)));
n(repmat(idx ,size(n_coat,1),1,1)) = n_coat(:);

% Make it symmetric
n = sqrt((n.^2+flip(n.^2,2))/2);
n = sqrt((n.^2+flip(n.^2,3))/2);

end