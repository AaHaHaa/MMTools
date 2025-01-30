function [n, x, dx] = build_GRIN_xy(Nx, spatial_window, core_radius, clad_radius, n, alpha)
%BUILD_GRIN_XY Build a GRIN's 2D refractive index profile
%
%   Nx - the number of spatial points in each dimension
%   spatial_window - the total size of space in each dimension, in um
%   core_radius - the radius of the core, in um
%   clad_radius - the radius of the cladding, in um
%   n - the refractive index of the core and cladding; a cell as {n_core,n_clad,n_coat}
%   alpha - the shape parameter for the GRIN fiber

% refractive index
[n_core, n_clad, n_coat] = deal(n{:});

dx = spatial_window/Nx; % um

x = (-floor(Nx/2):floor((Nx-1)/2))*dx;
[X, Y] = meshgrid(x, x);

%% GRIN index profile
n = max(n_clad, n_core - (n_core-n_clad)*(sqrt(X.^2+Y.^2)/core_radius).^alpha);
n(sqrt(X.^2+Y.^2)>clad_radius) = n_coat;

% Make it symmetric
n = sqrt((n.^2+flip(n.^2,2))/2);
n = sqrt((n.^2+flip(n.^2,3))/2);

end