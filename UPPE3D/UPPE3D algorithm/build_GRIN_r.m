function n = build_GRIN_r(r, core_radius, clad_radius, n, alpha)
%BUILD_GRIN_R Build a GRIN's 2D radially-symmetric refractive index profile
% of the fiber. Only 1D radial information is generated.
%
%   r - radial sampling positions
%   core_radius - the radius of the core, in um
%   clad_radius - the radius of the cladding, in um
%   n - the refractive index of the core and cladding; a cell as {n_core,n_clad,n_coat}
%   alpha - the shape parameter for the GRIN fiber

% refractive index
[n_core, n_clad, n_coat] = deal(n{:});

%% GRIN index profile
n = max(n_clad, n_core - (n_core-n_clad)*(r/core_radius).^alpha);
n(r>clad_radius) = n_coat;

end