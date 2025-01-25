function Ef_out = add_spherical_lens_phase_xy(Ef,dx,dy,wavelength,radius_of_curvature,varargin)
%ADD_SPERICALLENS_PHASE_XY It adds a spherical-lens phase for a 3D electric field
% The lens is assumed to be a plano-convex/concave lens with its curved
% side facing the beam.
%   Ef: electric field in the frequency domain(=ifft(Et,[],1)); size: (Nt,Nx,Ny,Nz)
%   dx: (m)
%   dy: (m)
%   wavelength: wavelength in the same order as Ef (m)
%   radius_of_curvature: focal length of a lens
%
% Optional input argument:
%   lens_material: (default: BK7)

%% Default optional input arguments
% Accept only 1 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 1
    error('add_spherical_lens_phase:TooManyInputs', ...
          'It takes only at most 1 optional input.');
end

% Set defaults for optional inputs
lens_material = 'BK7';
optargs = {lens_material};

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
lens_material = optargs{:};

%% Lens phase calculation
[Nt,Nx,Ny,Nz] = size(Ef);

x = (-Nx/2:Nx/2-1)*dx;
y = permute((-Ny/2:Ny/2-1)*dy,[1,3,2]);

h2 = x.^2+y.^2;
dr = abs(radius_of_curvature) - sqrt(radius_of_curvature^2 - h2);

beyond_lens = (radius_of_curvature^2 < h2);

% Sellmeier coefficients
[a,b] = Sellmeier_coefficients(lens_material); % assume only silica lens in this function
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
n_lens = n_from_Sellmeier(wavelength*1e6);
k0 = 2*pi./wavelength;
k_lens = k0.*n_lens;

lens_phase = (k0-k_lens).*dr*sign(radius_of_curvature);

Ef_out = Ef.*exp(1i*lens_phase);

beyond_lens = repmat(beyond_lens,[Nt,1,1,Nz]);
Ef_out(beyond_lens) = 0; % the beam is clipped

end

