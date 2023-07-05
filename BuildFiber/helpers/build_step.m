function [epsilon, x, dx] = build_step(Nx, spatial_window, core_radius, clad_radius, n)
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

% Step index profile
epsilon = n_clad^2*ones(Nx, Nx);
epsilon(X.^2+Y.^2 <= core_radius^2) = n_core^2;

epsilon(sqrt(X.^2+Y.^2)>clad_radius) = n_coat^2;

% Make it symmetric
epsilon = (epsilon+flipud(epsilon))/2;
epsilon = (epsilon+fliplr(epsilon))/2;

end