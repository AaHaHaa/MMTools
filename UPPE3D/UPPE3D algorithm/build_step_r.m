function n = build_step_r(r, core_radius, clad_radius, n)
%BUILD_STEP_R Build a step's 2D radially-symmetric refractive index profile
% of the fiber. Only 1D radial information is generated.
%
%   r - radial sampling positions
%   core_radius - the radius of the core, in um
%   clad_radius - the radius of the cladding, in um
%   n - the refractive index of the core and cladding; a cell as {n_core,n_clad,n_coat}

% refractive index
[n_core, n_clad, n_coat] = deal(n{:});

%% Step index profile
Nr = length(r);

% clad index
n = n_clad.*ones(1, Nr);

% core index
idx = r <= core_radius;
n_core = repmat(n_core,1,sum(idx(:)));
n(repmat(idx ,size(n_core,1),1)) = n_core(:);

% coat index
idx = r > clad_radius;
n_coat = repmat(n_coat,1,sum(idx(:)));
n(repmat(idx ,size(n_coat,1),1)) = n_coat(:);

end