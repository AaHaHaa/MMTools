function output = build_2Dgaussian_xy(MFD, spatial_window, Nx, varargin)
%BUILD_2DGAUSSIAN_XY Build a superGaussian spatial profile using the following parameters:
%
% MFD - mode-field diameter (m)
% spatial_window - width of the spatial domain (um)
% Nx - number of spatial grid points
%
% Optional inputs (varargin):
%   x_center - x spatial position (default is 0)
%   y_center - y spatial position (default is 0)
%   gaussexpo - supergaussian exponent (~exp(-t^(2*gaussexpo))) (default is 1)

%% Default optional input arguments
% Accept only 3 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 3
    error('build_2Dgaussian:TooManyInputs', ...
        'It takes only at most 3 optional inputs');
end

% Set defaults for optional inputs
optargs = {0,0,1};

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[x_center, y_center, gaussexpo] = optargs{:};

%% Gausian spatial profile
dx = spatial_window/Nx;  % m
x = (-Nx/2:Nx/2-1)'*dx; % m
y = x';
MFR = MFD/2; % mode-field radius; m

gexpo = 2*gaussexpo;

spatial_profile = exp(-sqrt((x-x_center).^2+(y-y_center).^2).^gexpo/(MFR^gexpo));
spatial_profile = spatial_profile/sqrt(sum(abs(spatial_profile(:)).^2)*dx^2);

%% Output as a struct
output = struct('field',spatial_profile,'dx',dx);

end